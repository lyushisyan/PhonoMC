#!/usr/bin/env python3
"""Convert scalar, full-grid phono3py collisions to a PhonoMC energy generator.

The supported raw input is phono3py's ``--write-collision`` output from
``--reducible-colmat --no-kappa-stars``, with every grid point and every band.
See ``material/phono3py_collision.md`` for units, preparation and limitations.
Only numpy and h5py are required. No scattering data are synthesized.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
import math
from pathlib import Path
import sys
import tempfile

import h5py
import numpy as np


# Exact SI constants; THz is cycles/ps, not angular frequency.
PLANCK_OVER_BOLTZMANN_THZ_K = 6.62607015e-34 * 1.0e12 / 1.380649e-23
RATE_FACTOR = 4.0 * np.pi


def _array(group, name: str, shape=None, *, integer=False) -> np.ndarray:
    if name not in group:
        raise ValueError(f"{group.filename}: missing '{name}' dataset")
    raw = np.asarray(group[name][...])
    if raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain real numbers")
    if not np.all(np.isfinite(raw)):
        raise ValueError(f"{name} contains nonfinite values")
    if shape is not None and raw.shape != shape:
        raise ValueError(f"{name} shape {raw.shape} must be {shape}")
    if integer and (np.any(raw != np.rint(raw)) or np.any(np.abs(raw) > 2**53)):
        raise ValueError(f"{name} must contain exactly representable integers")
    return np.asarray(raw, dtype=np.int64 if integer else np.float64)


def _temperature_index(group, temperature: float) -> int:
    temperatures = _array(group, "temperature")
    if temperatures.ndim != 1 or not temperatures.size:
        raise ValueError("temperature must be a nonempty vector")
    matched = np.flatnonzero(np.isclose(temperatures, temperature, rtol=0, atol=1e-6))
    if len(matched) != 1:
        raise ValueError(f"{group.filename}: temperature {temperature:g} K must occur exactly once")
    return int(matched[0])


def _mesh(group) -> np.ndarray:
    mesh = _array(group, "mesh", (3,), integer=True)
    if np.any(mesh <= 0):
        raise ValueError("mesh must contain three positive integers")
    # Generalized regular grids need a different GRG address transformation.
    if "Q_matrix" in group and not np.array_equal(_array(group, "Q_matrix"), np.eye(3)):
        raise ValueError("Generalized grids (nonidentity Q_matrix) are not supported")
    if "grid_matrix" in group:
        matrix = _array(group, "grid_matrix", (3, 3), integer=True)
        if not np.array_equal(matrix, np.diag(mesh)):
            raise ValueError("Generalized or nondiagonal grid_matrix is not supported")
    return mesh


def _grid_ids(qpoints: np.ndarray, mesh: np.ndarray, tolerance: float) -> np.ndarray:
    """Identify regular-grid points modulo reciprocal-lattice vectors."""
    nq = math.prod(int(x) for x in mesh)
    if qpoints.shape != (nq, 3) or not np.all(np.isfinite(qpoints)):
        raise ValueError("qpoint must contain exactly the full mesh, with shape (Nq,3)")
    nearest = np.rint(qpoints * mesh)
    if np.any(np.abs(qpoints - nearest / mesh) > tolerance):
        raise ValueError("qpoint lies outside the declared Gamma-centered regular mesh")
    addresses = nearest.astype(np.int64) % mesh
    ids = addresses[:, 0] + mesh[0] * (addresses[:, 1] + mesh[1] * addresses[:, 2])
    if np.unique(ids).size != nq:
        raise ValueError("qpoint contains duplicate points modulo reciprocal-lattice vectors")
    return ids


def energy_scale(frequency: np.ndarray, temperature: float) -> np.ndarray:
    """Return D up to a common factor: x/(2 sinh(x/2)), x=h*nu/(kB*T).

    D is proportional to h*nu*sqrt(n0*(n0+1)), equivalently sqrt(mode Cv).
    The exponential form is stable both for acoustic modes and large x.
    """
    x = np.asarray(frequency) * PLANCK_OVER_BOLTZMANN_THZ_K / temperature
    with np.errstate(over="ignore", under="ignore", invalid="ignore", divide="ignore"):
        scale = x * np.exp(-x / 2.0) / (-np.expm1(-x))
    if np.any(~np.isfinite(scale)) or np.any(scale <= np.finfo(float).tiny):
        raise ValueError("Energy transformation is numerically singular at this temperature")
    return scale


def _reject_diagonal_extras(group) -> None:
    for name in ("gamma_isotope", "gamma_elph"):
        if name in group and np.any(_array(group, name) != 0):
            raise ValueError(
                f"Nonzero {name} provides only a diagonal loss, not a conserving full "
                "scattering kernel. Generate pure phonon-phonon collisions."
            )
    if "boundary_mfp" in group and np.any(_array(group, "boundary_mfp") != 0):
        raise ValueError("boundary_mfp diagonal scattering is not supported; use PhonoMC boundaries")


def _check_collision_input(group, metadata) -> None:
    if "collision_matrix" not in group or group["collision_matrix"].dtype.kind not in "iuf":
        raise ValueError("collision_matrix must contain real numbers")
    if "band_index" in group:
        raise ValueError("Per-band collision fragments are unsupported; provide all bands per grid point")
    for key in ("sigma", "sigma_cutoff_width"):
        if (key in group) != (key in metadata):
            raise ValueError(f"Collision and collision-kappa {key} metadata differ")
        if key in group and not np.array_equal(_array(group, key), _array(metadata, key)):
            raise ValueError(f"Collision and collision-kappa {key} differ")


def _read_collision(paths, metadata, temperature: float, nq: int, nb: int,
                    input_format: str) -> tuple[np.ndarray, list[str]]:
    """Read raw whole tensors or all per-grid-point pieces without guessing axes."""
    size = nq * nb
    mt = _temperature_index(metadata, temperature)
    nt = len(metadata["temperature"])
    reference_gamma = None
    if input_format == "phono3py-raw":
        reference_gamma = _array(metadata, "gamma", (nt, nq, nb))[mt]
        if np.any(reference_gamma < 0):
            raise ValueError("gamma must be nonnegative")
        _reject_diagonal_extras(metadata)
    with h5py.File(paths[0], "r") as first:
        is_grid_piece = "grid_point" in first and first["grid_point"].shape == ()
        whole = len(paths) == 1 and not is_grid_piece
    versions = []
    if whole:
        with h5py.File(paths[0], "r") as group:
            _check_collision_input(group, metadata)
            ti = _temperature_index(group, temperature)
            expected = (len(group["temperature"]), nq, nb, nq, nb)
            if group["collision_matrix"].shape != expected:
                raise ValueError(
                    f"Scalar full-BZ collision_matrix must have shape {expected}; "
                    "irreducible Cartesian tensors, flattened arrays and partial bands are unsupported"
                )
            collision = np.asarray(group["collision_matrix"][ti], dtype=float).reshape(size, size)
            if input_format == "phono3py-raw":
                _reject_diagonal_extras(group)
                gamma = _array(group, "gamma", expected[:3])[ti]
                if not np.allclose(gamma, reference_gamma, rtol=1e-7, atol=1e-12):
                    raise ValueError("Raw gamma disagrees with collision-kappa: incomplete rows or different run/order")
                collision.flat[::size + 1] += gamma.reshape(-1)
            if "version" in group:
                versions.append(str(group["version"][()]))
    else:
        if input_format != "phono3py-raw":
            raise ValueError("Per-grid collision pieces require --input-format phono3py-raw")
        source_gp = _array(metadata, "grid_point", (nq,), integer=True)
        if np.unique(source_gp).size != nq:
            raise ValueError("collision-kappa grid_point identifiers must be unique")
        rows = {int(gp): i for i, gp in enumerate(source_gp)}
        seen = set()
        collision = np.zeros((nq, nb, nq, nb), dtype=float)
        for path in paths:
            with h5py.File(path, "r") as group:
                _check_collision_input(group, metadata)
                gp = int(_array(group, "grid_point", (), integer=True))
                if gp not in rows or gp in seen:
                    raise ValueError(f"Unknown or duplicate collision grid_point {gp}")
                seen.add(gp)
                row = rows[gp]
                ti = _temperature_index(group, temperature)
                expected = (len(group["temperature"]), nb, nq, nb)
                if group["collision_matrix"].shape != expected:
                    raise ValueError(f"Per-grid scalar collision_matrix shape must be {expected}; IR Cartesian data are unsupported")
                _reject_diagonal_extras(group)
                gamma = _array(group, "gamma", expected[:2])[ti]
                if not np.allclose(gamma, reference_gamma[row], rtol=1e-7, atol=1e-12):
                    raise ValueError(f"Raw gamma at grid_point {gp} disagrees with collision-kappa")
                collision[row] = group["collision_matrix"][ti]
                for branch in range(nb):
                    collision[row, branch, row, branch] += gamma[branch]
                if "version" in group:
                    versions.append(str(group["version"][()]))
        if len(seen) != nq:
            raise ValueError(f"Incomplete collision pieces: received {len(seen)} of {nq} grid points")
        collision = collision.reshape(size, size)
    if not np.all(np.isfinite(collision)):
        raise ValueError("collision_matrix contains nonfinite values")
    if len(set(versions)) > 1:
        raise ValueError("Collision pieces were generated by different phono3py versions")
    return collision, sorted(set(versions))


def _relative_column_residual(matrix: np.ndarray) -> float:
    norm = np.abs(matrix).sum(axis=0)
    # Exactly empty columns have zero residual; tiny but nonzero columns still count.
    residual = np.divide(np.abs(matrix.sum(axis=0)), norm,
                         out=np.zeros_like(norm), where=norm != 0)
    return float(np.max(residual, initial=0))


def make_generator(collision: np.ndarray, frequency: np.ndarray, temperature: float,
                   *, project_energy: bool = False, tolerance: float = 1e-10):
    """Return A=-4*pi*D*Omega/D, with A acting on column energy deviations.

    Symmetrization matches phono3py's full-grid LBTE assembly. The optional
    explicit P*Omega*P projection imposes the energy null vector D, preserving
    symmetry and positive semidefiniteness when the input has that property.
    """
    collision = (collision + collision.T) * 0.5
    scale = energy_scale(frequency, temperature)
    initial = -RATE_FACTOR * collision * scale[:, None] / scale[None, :]
    original_residual = _relative_column_residual(initial)
    projection_relative_norm = 0.0
    if project_energy:
        unit = scale / np.linalg.norm(scale)
        action = collision @ unit
        correction = (np.outer(unit, action) + np.outer(action, unit)
                      - float(unit @ action) * np.outer(unit, unit))
        norm = np.linalg.norm(collision)
        projection_relative_norm = float(np.linalg.norm(correction) / norm) if norm else 0.0
        collision -= correction
        initial = -RATE_FACTOR * collision * scale[:, None] / scale[None, :]
    if not np.all(np.isfinite(initial)):
        raise ValueError("Nonfinite values after energy-variable transformation")
    residual = _relative_column_residual(initial)
    if residual > tolerance:
        raise ValueError(
            f"Energy-conservation relative column residual {residual:.6g} exceeds {tolerance:g}. "
            "Check data completeness, temperature and units. To explicitly impose the energy "
            "sum rule, use --project-energy-conservation and inspect its reported correction."
        )
    diagonal = np.diag(initial)
    if np.any(diagonal > 1e-13 * np.abs(initial).sum(axis=0)):
        raise ValueError("Collision generator contains a positive growth diagonal")
    return initial, {"energy_residual_before_projection": original_residual,
                     "energy_residual_after_projection": residual,
                     "projection_relative_frobenius_norm": projection_relative_norm}


def convert(collision_paths, material_path, output_path, temperature: float, *,
            collision_kappa=None, input_format="phono3py-raw", project_energy=False,
            tolerance=1e-10, q_tolerance=1e-6, frequency_tolerance=1e-10,
            max_memory_gb=2.0, overwrite=False) -> dict:
    """Validate, assemble and atomically write a portable CSC collision file."""
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("Reference temperature must be finite and positive")
    for label, value in [("tolerance", tolerance), ("q_tolerance", q_tolerance),
                         ("frequency_tolerance", frequency_tolerance), ("max_memory_gb", max_memory_gb)]:
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{label} must be finite and positive")
    if tolerance > 1e-10 or frequency_tolerance > 1e-10:
        raise ValueError("Conservation and relative frequency tolerances may not exceed 1e-10")
    if input_format not in ("phono3py-raw", "phono3py-assembled"):
        raise ValueError("Unknown input_format")
    collision_paths = [Path(path).expanduser().resolve() for path in collision_paths]
    if not collision_paths or len(set(collision_paths)) != len(collision_paths):
        raise ValueError("Provide nonempty, distinct collision input paths")
    material_path, output_path = Path(material_path).resolve(), Path(output_path).resolve()
    metadata_path = Path(collision_kappa).resolve() if collision_kappa else collision_paths[0]
    if output_path in [*collision_paths, material_path, metadata_path]:
        raise ValueError("Output must not overwrite an input file")
    if output_path.exists() and not overwrite:
        raise ValueError("Output already exists; use --overwrite to replace it")
    with ExitStack() as stack:
        material = stack.enter_context(h5py.File(material_path, "r"))
        metadata = stack.enter_context(h5py.File(metadata_path, "r"))
        mesh = _mesh(material)
        if not np.array_equal(mesh, _mesh(metadata)):
            raise ValueError("Material and collision-kappa meshes differ")
        nq = math.prod(int(x) for x in mesh)
        target_frequency = _array(material, "frequency")
        if target_frequency.ndim != 2 or target_frequency.shape[0] != nq or not target_frequency.shape[1]:
            raise ValueError("Material frequency must have shape (Nq,Nband)")
        nb = target_frequency.shape[1]
        # Includes several dense working arrays and worst-case CSC buffers.
        estimated_gb = (nq * nb)**2 * 8 * 8 / 1024**3
        if estimated_gb > max_memory_gb:
            raise ValueError(f"Dense conversion needs approximately {estimated_gb:.2f} GiB; "
                             f"limit is {max_memory_gb:g} GiB. Use a coarser mesh or increase --max-memory-gb")
        for group in (material, metadata):
            if "weight" not in group or not np.array_equal(_array(group, "weight", (nq,)), np.ones(nq)):
                raise ValueError("Both material and collision-kappa require full-grid unit weights")
        source_ids = _grid_ids(_array(metadata, "qpoint"), mesh, q_tolerance)
        target_ids = _grid_ids(_array(material, "qpoint"), mesh, q_tolerance)
        if input_format == "phono3py-raw" and not np.array_equal(source_ids, np.arange(nq)):
            raise ValueError("Raw collision-kappa must use phono3py GRG ordering (first mesh index fastest)")
        source_frequency = _array(metadata, "frequency", (nq, nb))
        target_by_id = np.argsort(target_ids)
        source_to_target = target_by_id[source_ids]
        if not np.allclose(source_frequency, target_frequency[source_to_target], rtol=frequency_tolerance, atol=0):
            raise ValueError("Frequencies disagree after q-point mapping; material, band ordering or calculation differs")
        _temperature_index(material, temperature)
        collision, versions = _read_collision(collision_paths, metadata, temperature, nq, nb, input_format)

    source_flat_frequency = source_frequency.reshape(-1)
    positive = source_flat_frequency > 0
    if not np.any(positive):
        raise ValueError("No positive-frequency modes")
    if np.any(~positive):
        scale = np.max(np.abs(collision), initial=0)
        inactive_coupling = max(np.max(np.abs(collision[~positive, :]), initial=0),
                                np.max(np.abs(collision[:, ~positive]), initial=0))
        if inactive_coupling > tolerance * scale:
            raise ValueError("Nonpositive-frequency modes have nonzero collision rows/columns; cannot discard them")
        collision = collision[np.ix_(positive, positive)]
    frequency = source_flat_frequency[positive]
    generator, diagnostics = make_generator(collision, frequency, temperature,
                                             project_energy=project_energy, tolerance=tolerance)
    mode_indices = np.column_stack((np.repeat(source_to_target, nb), np.tile(np.arange(nb), nq)))[positive]
    # Canonical output follows material q-major/branch-minor ordering.
    permutation = np.lexsort((mode_indices[:, 1], mode_indices[:, 0]))
    mode_indices, frequency = mode_indices[permutation], frequency[permutation]
    generator = generator[np.ix_(permutation, permutation)]
    size = len(frequency)
    offsets = np.empty(size + 1, dtype=np.int64)
    offsets[0] = 0
    rows, values = [], []
    for column in range(size):
        indices = np.flatnonzero(generator[:, column] != 0)
        rows.append(indices.astype(np.int64))
        values.append(generator[indices, column])
        offsets[column + 1] = offsets[column] + len(indices)
    diagnostics.update(modes=size, nonzeros=int(offsets[-1]), reference_temperature=temperature,
                       omitted_nonpositive_modes=int(np.count_nonzero(~positive)))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix=output_path.name + ".", suffix=".tmp",
                                     dir=output_path.parent, delete=False) as temporary:
        temporary_path = Path(temporary.name)
    try:
        with h5py.File(temporary_path, "w") as output:
            output.create_dataset("schema_version", data=np.int64(1))
            output.create_dataset("reference_temperature", data=float(temperature))
            output.create_dataset("mode_indices", data=mode_indices.astype(np.int64))
            output.create_dataset("frequency", data=frequency)
            output.create_dataset("column_offsets", data=offsets)
            output.create_dataset("row_indices", data=np.concatenate(rows))
            output.create_dataset("values", data=np.concatenate(values))
            output.attrs["representation"] = "energy_deviation_generator"
            output.attrs["equation"] = "df/dt = A f; column source, row destination"
            output.attrs["units"] = "ps^-1"
            output.attrs["source_format"] = input_format
            output.attrs["source_files"] = "\n".join(str(path) for path in collision_paths)
            output.attrs["collision_metadata"] = str(metadata_path)
            output.attrs["material_file"] = str(material_path)
            output.attrs["phono3py_versions"] = ", ".join(versions)
            output.attrs["energy_projection_applied"] = np.int8(project_energy)
            for key, value in diagnostics.items():
                output.attrs[key] = value
        # The second check protects against an output created during conversion.
        if output_path.exists() and not overwrite:
            raise ValueError("Output was created during conversion; refusing to overwrite")
        temporary_path.replace(output_path)
    finally:
        temporary_path.unlink(missing_ok=True)
    return diagnostics


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--collision", nargs="+", required=True, help="Whole collision HDF5 or all per-grid-point files")
    parser.add_argument("--material", required=True, help="PhonoMC kappa-fbz.hdf5")
    parser.add_argument("--collision-kappa", help="Same-run all-grid kappa file supplying GRG qpoint, frequency, gamma and grid_point")
    parser.add_argument("--temperature", type=float, required=True, help="Exact stored temperature in K; no interpolation")
    parser.add_argument("--out", required=True, help="Portable collision HDF5 output")
    parser.add_argument("--input-format", choices=["phono3py-raw", "phono3py-assembled"], required=True,
                        help="raw: add stored gamma diagonal; assembled: already-combined scalar Omega, before diagonalization")
    parser.add_argument("--project-energy-conservation", action="store_true",
                        help="Explicitly apply P*Omega*P energy sum-rule projection and report its correction")
    parser.add_argument("--conservation-tolerance", type=float, default=1e-10)
    parser.add_argument("--qpoint-tolerance", type=float, default=1e-6, help="Fractional-coordinate tolerance")
    parser.add_argument("--frequency-tolerance", type=float, default=1e-10, help="Relative frequency tolerance, at most 1e-10")
    parser.add_argument("--max-memory-gb", type=float, default=2.0, help="Maximum estimated dense working memory in GiB")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = convert(args.collision, args.material, args.out, args.temperature,
                         collision_kappa=args.collision_kappa, input_format=args.input_format,
                         project_energy=args.project_energy_conservation,
                         tolerance=args.conservation_tolerance, q_tolerance=args.qpoint_tolerance,
                         frequency_tolerance=args.frequency_tolerance,
                         max_memory_gb=args.max_memory_gb, overwrite=args.overwrite)
    except (ValueError, KeyError, OSError, MemoryError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    print(f"Wrote {args.out}: {result['modes']} modes, {result['nonzeros']} entries at {args.temperature:g} K")
    print(f"Energy residual: {result['energy_residual_before_projection']:.6g} -> "
          f"{result['energy_residual_after_projection']:.6g}; "
          f"projection relative norm: {result['projection_relative_frobenius_norm']:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
