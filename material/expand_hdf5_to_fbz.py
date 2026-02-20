#!/usr/bin/env python3
"""Expand IBZ phono3py kappa HDF5 into an explicit FBZ HDF5.

This script follows a symmetry-based route similar to the provided Python code:
1) Build reciprocal symmetry operations from POSCAR via phonopy.
2) Rotate each IBZ q-point to its star.
3) Fold rotated q-points into the first Brillouin zone (Wigner-Seitz-like nearest-|k| rule).
4) Remove duplicates in each star and assemble explicit FBZ q-list.
5) Expand q-indexed datasets consistently and write a new HDF5 file.

Main target is C++ reader compatibility (frequency/group_velocity/gamma/qpoint/mesh/weight/temperature).
"""

from __future__ import annotations

import argparse
import pathlib
import sys
from typing import Optional

import h5py
import numpy as np
from phonopy import Phonopy
from phonopy.interface.calculator import read_crystal_structure


def move_to_FBZ(q_in: np.ndarray, reciprocal_lattice: np.ndarray) -> np.ndarray:
    """Fold q-points into FBZ by minimizing |k+G| over G in {-1,0,1}^3."""
    shifts_int = np.array(list(np.ndindex(3, 3, 3)), dtype=np.int64) - 1
    recip_T = reciprocal_lattice.T

    q_current = np.asarray(q_in, dtype=np.float64).copy()
    k_current = q_current @ recip_T
    G_shifts_cart = shifts_int @ recip_T

    active = np.arange(q_current.shape[0], dtype=np.int64)
    for _ in range(3):
        if active.size == 0:
            break
        k_active = k_current[active]
        k_candidates = k_active[:, None, :] + G_shifts_cart[None, :, :]
        dists = np.linalg.norm(k_candidates, axis=-1)
        min_idx = np.argmin(dists, axis=1)
        best_shifts = shifts_int[min_idx]
        is_centered = np.all(best_shifts == 0, axis=1)

        q_current[active] += best_shifts
        k_current[active] += G_shifts_cart[min_idx]
        active = active[~is_centered]

    return q_current


def _infer_q_axis(shape: tuple[int, ...], nq_in: int) -> Optional[int]:
    candidates = [i for i, s in enumerate(shape) if s == nq_in]
    if len(candidates) == 1:
        return candidates[0]
    return None


def build_fbz_mapping(
    qpoints_ibz: np.ndarray,
    weight_ibz: np.ndarray,
    rotations: np.ndarray,
    reciprocal_lattice: np.ndarray,
    decimals: int = 6,
    enforce_weight: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (q_fbz, src_q_index, rot_index, star_size_per_ibz)."""
    q_out_parts = []
    src_parts = []
    rot_parts = []
    star_sizes = np.zeros(qpoints_ibz.shape[0], dtype=np.int64)

    for i, q in enumerate(qpoints_ibz):
        star = np.array([r @ q for r in rotations], dtype=np.float64)
        star = move_to_FBZ(star, reciprocal_lattice)
        star = np.around(star, decimals=decimals)

        star_q, unique_idx = np.unique(star, axis=0, return_index=True)
        raw_size = star_q.shape[0]
        target_size = raw_size
        if enforce_weight and weight_ibz is not None and i < weight_ibz.size:
            wi = int(weight_ibz[i])
            if wi > 0:
                target_size = wi
                if raw_size > wi:
                    star_q = star_q[:wi]
                    unique_idx = unique_idx[:wi]
                elif raw_size < wi and raw_size > 0:
                    rep = np.resize(np.arange(raw_size, dtype=np.int64), wi)
                    star_q = star_q[rep]
                    unique_idx = unique_idx[rep]
        star_sizes[i] = star_q.shape[0]

        q_out_parts.append(star_q)
        src_parts.append(np.full(star_q.shape[0], i, dtype=np.int64))
        rot_parts.append(unique_idx.astype(np.int64))

    q_fbz = np.concatenate(q_out_parts, axis=0)
    src_idx = np.concatenate(src_parts, axis=0)
    rot_idx = np.concatenate(rot_parts, axis=0)

    if weight_ibz is not None and weight_ibz.size == star_sizes.size:
        mismatch = np.where(star_sizes != weight_ibz.astype(np.int64))[0]
        if mismatch.size > 0:
            print(
                f"[warn] star-size and weight mismatch at {mismatch.size} IBZ points "
                f"(first index: {int(mismatch[0])})"
            )

    return q_fbz, src_idx, rot_idx, star_sizes


def rotate_group_velocity(gv_ibz: np.ndarray, src_idx: np.ndarray, rot_idx: np.ndarray, Rcart_list: list[np.ndarray]) -> np.ndarray:
    """Expand and rotate group_velocity from IBZ to FBZ.

    gv_ibz shape: (Nq_ibz, Nband, 3)
    output shape: (Nq_fbz, Nband, 3)
    """
    out = np.empty((src_idx.size, gv_ibz.shape[1], 3), dtype=np.float64)
    for j, (si, ri) in enumerate(zip(src_idx, rot_idx, strict=False)):
        v = gv_ibz[int(si)]
        R = Rcart_list[int(ri)]
        out[j] = (R @ v.T).T
    return out


def maybe_plot(
    outdir: pathlib.Path,
    q_fbz: np.ndarray,
    reciprocal_lattice: np.ndarray,
    frequency_fbz: np.ndarray,
    gv_fbz: np.ndarray,
    gamma_fbz: np.ndarray,
    temperatures: np.ndarray,
    T0: float,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"[warn] matplotlib unavailable, skip plotting: {exc}")
        return

    outdir.mkdir(parents=True, exist_ok=True)

    k = q_fbz @ reciprocal_lattice.T
    x, y = k[:, 0], k[:, 1]

    def _save_six(values: np.ndarray, title: str, cbar: str, fn: str) -> None:
        nband = values.shape[1]
        nb = min(6, nband)
        vmin, vmax = np.nanmin(values[:, :nb]), np.nanmax(values[:, :nb])
        fig, axes = plt.subplots(2, 3, figsize=(11, 7), dpi=220, constrained_layout=True)
        axes = axes.ravel()
        sc = None
        for b in range(nb):
            sc = axes[b].scatter(x, y, c=values[:, b], s=2.5, alpha=0.7, cmap="viridis", vmin=vmin, vmax=vmax)
            axes[b].set_aspect("equal", adjustable="box")
            axes[b].set_title(f"Branch {b}")
            axes[b].set_xlabel("$k_x$")
            axes[b].set_ylabel("$k_y$")
        for b in range(nb, 6):
            axes[b].axis("off")
        if sc is not None:
            c = fig.colorbar(sc, ax=axes[:nb], shrink=0.9)
            c.set_label(cbar)
        fig.suptitle(title)
        fp = outdir / fn
        fig.savefig(fp, bbox_inches="tight")
        plt.close(fig)
        print(f"[ok] saved -> {fp}")

    _save_six(frequency_fbz, "Frequency (THz)", "f (THz)", "FBZ_freq.png")
    speed = np.linalg.norm(gv_fbz, axis=-1)
    _save_six(speed, "Group velocity magnitude", "|v|", "FBZ_speed.png")

    t_idx = int(np.argmin(np.abs(temperatures - T0))) if temperatures.size else 0
    gT = gamma_fbz[t_idx]
    tau = np.zeros_like(gT)
    m = gT > 1e-12
    tau[m] = 1.0 / (2.0 * gT[m])
    _save_six(tau, f"Tau at T~{temperatures[t_idx]:.1f} K", "tau (ps)", "FBZ_tau.png")


def maybe_export_csv(out_csv: Optional[pathlib.Path], q_fbz: np.ndarray, frequency_fbz: np.ndarray, gv_fbz: np.ndarray, gamma_fbz: np.ndarray, heat_capacity_fbz: Optional[np.ndarray], temperatures: np.ndarray, T0: float) -> None:
    if out_csv is None:
        return
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    Nq, Nm = frequency_fbz.shape
    t_idx = int(np.argmin(np.abs(temperatures - T0))) if temperatures.size else 0

    gT = gamma_fbz[t_idx]
    if heat_capacity_fbz is not None:
        Cv = heat_capacity_fbz[t_idx]
    else:
        Cv = np.zeros_like(gT)

    q_idx = np.repeat(np.arange(Nq, dtype=np.int64), Nm)
    s_idx = np.tile(np.arange(Nm, dtype=np.int64), Nq)

    rows = np.column_stack([
        q_idx,
        s_idx,
        np.repeat(q_fbz[:, 0], Nm),
        np.repeat(q_fbz[:, 1], Nm),
        np.repeat(q_fbz[:, 2], Nm),
        frequency_fbz.reshape(-1),
        gv_fbz[..., 0].reshape(-1),
        gv_fbz[..., 1].reshape(-1),
        gv_fbz[..., 2].reshape(-1),
        gT.reshape(-1),
        Cv.reshape(-1),
    ])

    header = f"q_idx,s_idx,qx,qy,qz,freq_THz,vx,vy,vz,gamma_THz,Cv; T={float(temperatures[t_idx]):.1f}K"
    np.savetxt(out_csv, rows, delimiter=",", header=header, comments="", fmt="%d,%d,%.6f,%.6f,%.6f,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g")
    print(f"[ok] CSV saved -> {out_csv}")


def run(args: argparse.Namespace) -> int:
    hdf_path = pathlib.Path(args.hdf).expanduser().resolve()
    poscar_path = pathlib.Path(args.poscar).expanduser().resolve() if args.poscar else hdf_path.with_name("POSCAR")
    out_path = pathlib.Path(args.out).expanduser().resolve() if args.out else hdf_path.with_name(f"{hdf_path.stem}-fbz.hdf5")

    if not hdf_path.exists():
        print(f"Error: HDF5 not found: {hdf_path}", file=sys.stderr)
        return 1
    if not poscar_path.exists():
        print(f"Error: POSCAR not found: {poscar_path}", file=sys.stderr)
        return 1

    unitcell, _ = read_crystal_structure(str(poscar_path), interface_mode="vasp")
    phonon = Phonopy(unitcell, supercell_matrix=np.eye(3, dtype=int), primitive_matrix="auto", symprec=1e-5)

    lattice_prim = np.array(phonon.primitive.cell, dtype=float)
    reciprocal_lattice = 2.0 * np.pi * np.linalg.inv(lattice_prim)
    rotations = np.asarray(phonon.primitive_symmetry.reciprocal_operations, dtype=np.int64)

    Rinvt = np.linalg.inv(reciprocal_lattice)
    Rcart_list = [reciprocal_lattice @ r @ Rinvt for r in rotations]

    with h5py.File(hdf_path, "r") as fin:
        if "qpoint" not in fin or "weight" not in fin:
            print("Error: input must contain 'qpoint' and 'weight' datasets", file=sys.stderr)
            return 1
        if "frequency" not in fin or "group_velocity" not in fin or "gamma" not in fin:
            print("Error: input missing one of required datasets: frequency/group_velocity/gamma", file=sys.stderr)
            return 1

        q_ibz = np.asarray(fin["qpoint"][...], dtype=np.float64)
        w_ibz = np.asarray(fin["weight"][...], dtype=np.int64)
        if q_ibz.ndim != 2 or q_ibz.shape[1] != 3:
            print(f"Error: qpoint shape invalid: {q_ibz.shape}", file=sys.stderr)
            return 1
        if w_ibz.ndim != 1 or w_ibz.shape[0] != q_ibz.shape[0]:
            print(f"Error: weight shape invalid: {w_ibz.shape}", file=sys.stderr)
            return 1

        q_fbz, src_idx, rot_idx, star_sizes = build_fbz_mapping(
            q_ibz,
            w_ibz,
            rotations,
            reciprocal_lattice,
            decimals=args.decimals,
            enforce_weight=args.enforce_weight,
        )

        nq_in = q_ibz.shape[0]
        nq_out = q_fbz.shape[0]
        weight_sum = int(np.sum(w_ibz))
        print(f"[info] nq_in={nq_in}, nq_out={nq_out}, sum(weight)={weight_sum}")
        if nq_out != weight_sum:
            print("[warn] nq_out != sum(weight). This can happen at BZ boundaries after folding/dedup.")

        out_path.parent.mkdir(parents=True, exist_ok=True)
        freq_fbz = None
        gv_fbz = None
        gamma_fbz = None
        hc_fbz = None
        with h5py.File(out_path, "w") as fout:
            for k, v in fin.attrs.items():
                fout.attrs[k] = v
            fout.attrs["expanded_from_ibz"] = np.int8(1)
            fout.attrs["expanded_source"] = str(hdf_path)
            fout.attrs["nq_in"] = np.int64(nq_in)
            fout.attrs["nq_out"] = np.int64(nq_out)
            fout.attrs["weight_sum_in"] = np.int64(weight_sum)
            fout.attrs["decimals"] = np.int64(args.decimals)
            fout.attrs["symmetry_ops"] = np.int64(rotations.shape[0])

            for name in fin.keys():
                arr = np.asarray(fin[name][...])
                if name == "qpoint":
                    out_arr = q_fbz
                elif name == "weight":
                    out_arr = np.ones(nq_out, dtype=arr.dtype)
                elif name == "group_velocity":
                    out_arr = rotate_group_velocity(arr, src_idx, rot_idx, Rcart_list)
                else:
                    axis = _infer_q_axis(arr.shape, nq_in) if arr.ndim > 0 else None
                    if axis is None:
                        out_arr = arr
                    else:
                        out_arr = np.take(arr, src_idx, axis=axis)

                ds_in = fin[name]
                kwargs = {}
                if ds_in.compression is not None:
                    kwargs["compression"] = ds_in.compression
                    if ds_in.compression_opts is not None:
                        kwargs["compression_opts"] = ds_in.compression_opts
                if ds_in.shuffle:
                    kwargs["shuffle"] = True
                if ds_in.fletcher32:
                    kwargs["fletcher32"] = True

                ds_out = fout.create_dataset(name, data=out_arr, **kwargs)
                for ak, av in ds_in.attrs.items():
                    ds_out.attrs[ak] = av
                if name == "frequency":
                    freq_fbz = out_arr
                elif name == "group_velocity":
                    gv_fbz = out_arr
                elif name == "gamma":
                    gamma_fbz = out_arr
                elif name == "heat_capacity":
                    hc_fbz = out_arr

            # write helper mapping datasets for debugging/repro
            fout.create_dataset("_fbz_src_q_index", data=src_idx)
            fout.create_dataset("_fbz_rot_index", data=rot_idx)
            fout.create_dataset("_fbz_star_size", data=star_sizes)

        print(f"[ok] FBZ HDF5 saved -> {out_path}")

        # optional CSV export
        temperatures = np.asarray(fin["temperature"][...]) if "temperature" in fin else np.array([args.T0], dtype=float)
        if freq_fbz is None or gv_fbz is None or gamma_fbz is None:
            print("Error: failed to build required expanded datasets.", file=sys.stderr)
            return 1

    outdir = pathlib.Path(args.outdir).expanduser().resolve() if args.outdir else out_path.parent

    if args.csv:
        csv_path = pathlib.Path(args.csv)
        if not csv_path.is_absolute():
            csv_path = outdir / csv_path
        maybe_export_csv(csv_path, q_fbz, freq_fbz, gv_fbz, gamma_fbz, hc_fbz, temperatures, args.T0)

    if args.plot:
        maybe_plot(outdir, q_fbz, reciprocal_lattice, freq_fbz, gv_fbz, gamma_fbz, temperatures, args.T0)

    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Expand IBZ kappa HDF5 to FBZ using phonopy symmetry.")
    p.add_argument("--hdf", required=True, help="Input HDF5 path (e.g. kappa.hdf5)")
    p.add_argument("--poscar", default="", help="POSCAR path (default: sibling POSCAR of HDF)")
    p.add_argument("--out", default="", help="Output HDF5 path (default: <input>-fbz.hdf5)")
    p.add_argument("--decimals", type=int, default=6, help="Rounding decimals before unique")
    p.add_argument("--no-enforce-weight", dest="enforce_weight", action="store_false", help="Do not force each star size to input weight")
    p.set_defaults(enforce_weight=True)
    p.add_argument("--T0", type=float, default=300.0, help="Reference temperature for tau/CSV")
    p.add_argument("--csv", default="", help="Optional CSV output path (relative path uses --outdir)")
    p.add_argument("--plot", action="store_true", help="Generate FBZ_freq/speed/tau png files")
    p.add_argument("--outdir", default="", help="Output directory for CSV/plots (default: output hdf folder)")
    return p


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
