#!/usr/bin/env python3
"""Plot 3D temperature distribution from NTMC convergence output."""

from __future__ import annotations

import argparse
import math
import re
import struct
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm, colors

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError as exc:  # pragma: no cover
    raise SystemExit("Python 3.11+ is required (tomllib missing).") from exc

Vec3 = Tuple[float, float, float]
Tri = Tuple[int, int, int]


def _load_toml(path: Path) -> dict:
    with path.open("rb") as f:
        return tomllib.load(f)


def _resolve_results_dir(input_path: Path, cfg: dict, explicit: str | None) -> Path:
    if explicit:
        p = Path(explicit)
        if not p.is_absolute():
            p = (input_path.parent / p).resolve()
        return p

    io_cfg = cfg.get("io", {})
    base = io_cfg.get("output_folder")
    if not base:
        raise ValueError("`[io].output_folder` not found in input TOML, and --results was not provided.")
    base_path = Path(base)
    if not base_path.is_absolute():
        base_path = (input_path.parent / base_path).resolve()

    if (base_path / "convergence.txt").exists():
        return base_path

    parent = base_path.parent
    stem = base_path.name
    pat = re.compile(rf"^{re.escape(stem)}_(\d+)$")
    best_idx = -1
    best_dir: Path | None = None
    if parent.exists():
        for p in parent.iterdir():
            if not p.is_dir():
                continue
            m = pat.match(p.name)
            if not m:
                continue
            idx = int(m.group(1))
            if idx > best_idx and (p / "convergence.txt").exists():
                best_idx = idx
                best_dir = p
    if best_dir is None:
        raise FileNotFoundError(f"Cannot find results folder from base: {base_path}")
    return best_dir


def _load_convergence(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing convergence file: {path}")

    lines = path.read_text(encoding="utf-8").splitlines()
    header_line = None
    for ln in lines:
        if ln.startswith("#"):
            header_line = ln
            break
    if header_line is None:
        raise ValueError("convergence.txt header not found.")

    cols = header_line.lstrip("#").strip().split()
    if "heatflux" not in cols:
        raise ValueError("Invalid convergence header: missing 'heatflux'.")
    heatflux_idx = cols.index("heatflux")
    temp_cols = cols[1:heatflux_idx]

    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    timesteps = data[:, 0]
    temps = data[:, 1:heatflux_idx]
    rest = data[:, heatflux_idx:]
    return timesteps, temps, rest, temp_cols


def _triangulate_polygon(indices: Sequence[int]) -> List[Tri]:
    if len(indices) < 3:
        return []
    tris: List[Tri] = []
    for i in range(1, len(indices) - 1):
        tris.append((indices[0], indices[i], indices[i + 1]))
    return tris


def _load_obj(path: Path) -> Tuple[List[Vec3], List[Tri]]:
    vertices: List[Vec3] = []
    faces: List[Tri] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if parts[0] == "v" and len(parts) >= 4:
            vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
        elif parts[0] == "f" and len(parts) >= 4:
            idx: List[int] = []
            for tk in parts[1:]:
                v = tk.split("/")[0]
                ii = int(v)
                if ii < 0:
                    ii = len(vertices) + ii
                else:
                    ii = ii - 1
                idx.append(ii)
            faces.extend(_triangulate_polygon(idx))
    if not vertices or not faces:
        raise ValueError(f"OBJ mesh has no vertices/faces: {path}")
    return vertices, faces


def _is_binary_stl(path: Path) -> bool:
    raw = path.read_bytes()
    if len(raw) < 84:
        return False
    ntri = struct.unpack_from("<I", raw, 80)[0]
    expected = 84 + 50 * ntri
    if expected == len(raw):
        return True
    head = raw[:5].decode("ascii", errors="ignore").lower()
    return head != "solid"


def _load_binary_stl(path: Path) -> Tuple[List[Vec3], List[Tri]]:
    raw = path.read_bytes()
    ntri = struct.unpack_from("<I", raw, 80)[0]
    vertices: List[Vec3] = []
    faces: List[Tri] = []
    vmap: dict[Tuple[int, int, int], int] = {}

    def add_vertex(v: Vec3) -> int:
        key = (round(v[0] * 1e9), round(v[1] * 1e9), round(v[2] * 1e9))
        if key in vmap:
            return vmap[key]
        idx = len(vertices)
        vmap[key] = idx
        vertices.append(v)
        return idx

    off = 84
    for _ in range(ntri):
        vals = struct.unpack_from("<12fH", raw, off)
        off += 50
        p0 = (float(vals[3]), float(vals[4]), float(vals[5]))
        p1 = (float(vals[6]), float(vals[7]), float(vals[8]))
        p2 = (float(vals[9]), float(vals[10]), float(vals[11]))
        i0 = add_vertex(p0)
        i1 = add_vertex(p1)
        i2 = add_vertex(p2)
        faces.append((i0, i1, i2))
    if not vertices or not faces:
        raise ValueError(f"Binary STL mesh has no vertices/faces: {path}")
    return vertices, faces


def _load_ascii_stl(path: Path) -> Tuple[List[Vec3], List[Tri]]:
    vertices: List[Vec3] = []
    faces: List[Tri] = []
    vmap: dict[Tuple[int, int, int], int] = {}
    tri_tmp: List[int] = []

    def add_vertex(v: Vec3) -> int:
        key = (round(v[0] * 1e9), round(v[1] * 1e9), round(v[2] * 1e9))
        if key in vmap:
            return vmap[key]
        idx = len(vertices)
        vmap[key] = idx
        vertices.append(v)
        return idx

    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        s = raw.strip()
        if not s.lower().startswith("vertex"):
            continue
        parts = s.split()
        if len(parts) < 4:
            continue
        idx = add_vertex((float(parts[1]), float(parts[2]), float(parts[3])))
        tri_tmp.append(idx)
        if len(tri_tmp) == 3:
            faces.append((tri_tmp[0], tri_tmp[1], tri_tmp[2]))
            tri_tmp.clear()
    if not vertices or not faces:
        raise ValueError(f"ASCII STL mesh has no vertices/faces: {path}")
    return vertices, faces


def _load_mesh(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    ext = path.suffix.lower()
    if ext == ".obj":
        verts, faces = _load_obj(path)
    elif ext == ".stl":
        if _is_binary_stl(path):
            verts, faces = _load_binary_stl(path)
        else:
            verts, faces = _load_ascii_stl(path)
    else:
        raise ValueError(f"Unsupported mesh file: {path}")
    v = np.array(verts, dtype=float)
    f = np.array(faces, dtype=int)
    # Match C++ behavior: shift mesh to origin.
    v = v - v.min(axis=0, keepdims=True)
    return v, f


def _ray_intersects_triangle(orig: np.ndarray, direction: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float | None:
    eps = 1e-12
    e1 = b - a
    e2 = c - a
    pvec = np.cross(direction, e2)
    det = float(np.dot(e1, pvec))
    if abs(det) < eps:
        return None
    inv_det = 1.0 / det
    tvec = orig - a
    u = float(np.dot(tvec, pvec) * inv_det)
    if u < -eps or u > 1.0 + eps:
        return None
    qvec = np.cross(tvec, e1)
    v = float(np.dot(direction, qvec) * inv_det)
    if v < -eps or (u + v) > 1.0 + eps:
        return None
    t = float(np.dot(e2, qvec) * inv_det)
    if t <= eps:
        return None
    return t


def _contains_point_mesh(point: np.ndarray, vertices: np.ndarray, faces: np.ndarray) -> bool:
    # Use +x ray like C++ and deduplicate equal-t hits from triangle pairs on one facet.
    direction = np.array([1.0, 0.0, 0.0], dtype=float)
    ts: List[float] = []
    for i0, i1, i2 in faces:
        t = _ray_intersects_triangle(point, direction, vertices[i0], vertices[i1], vertices[i2])
        if t is not None:
            ts.append(t)
    if not ts:
        return False
    ts.sort()
    unique_hits = 1
    tol = 1e-9
    for i in range(1, len(ts)):
        if abs(ts[i] - ts[i - 1]) > tol:
            unique_hits += 1
    return (unique_hits % 2) == 1


def _build_centers_from_grid(cfg: dict, input_path: Path) -> np.ndarray:
    sim = cfg.get("simulation", {})
    g = sim.get("grid_xyz")
    if not isinstance(g, list) or len(g) != 3:
        raise ValueError("`grid_xyz` must be [nx, ny, nz].")
    nx, ny, nz = int(g[0]), int(g[1]), int(g[2])
    if nx <= 0 or ny <= 0 or nz <= 0:
        raise ValueError("grid_xyz entries must be positive.")

    geo = cfg.get("geometry", {})
    model = geo.get("model")
    if not isinstance(model, str):
        raise ValueError("Missing `[geometry].model` in input TOML.")

    if model == "box":
        dims = geo.get("sizes")
        if not isinstance(dims, list) or len(dims) != 3:
            raise ValueError("Box model requires `sizes = [Lx, Ly, Lz]`.")
        bmin = np.array([0.0, 0.0, 0.0], dtype=float)
        # Input sizes are in nm; solver converts to Angstrom internally.
        bmax = np.array([float(dims[0]), float(dims[1]), float(dims[2])], dtype=float) * 10.0
        ext = bmax - bmin
        centers: List[Vec3] = []
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    p = (
                        float(bmin[0] + (ix + 0.5) * ext[0] / nx),
                        float(bmin[1] + (iy + 0.5) * ext[1] / ny),
                        float(bmin[2] + (iz + 0.5) * ext[2] / nz),
                    )
                    centers.append(p)
        return np.array(centers, dtype=float)

    mesh_path = Path(model)
    if not mesh_path.is_absolute():
        p1 = (Path.cwd() / mesh_path).resolve()
        p2 = (input_path.parent / mesh_path).resolve()
        mesh_path = p1 if p1.exists() else p2
    vertices, faces = _load_mesh(mesh_path)
    if mesh_path.suffix.lower() == ".stl":
        # STL geometry is in nm; solver converts to Angstrom internally.
        vertices = vertices * 10.0
    bmin = vertices.min(axis=0)
    bmax = vertices.max(axis=0)
    ext = bmax - bmin
    centers: List[Vec3] = []

    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                p = np.array(
                    [
                        bmin[0] + (ix + 0.5) * ext[0] / nx,
                        bmin[1] + (iy + 0.5) * ext[1] / ny,
                        bmin[2] + (iz + 0.5) * ext[2] / nz,
                    ],
                    dtype=float,
                )
                if _contains_point_mesh(p, vertices, faces):
                    centers.append((float(p[0]), float(p[1]), float(p[2])))
    return np.array(centers, dtype=float)


def _load_centers_csv(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(path)
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 4:
        raise ValueError(f"Invalid centers csv format: {path}")
    return data[:, 1:4]


def _set_equal_axes(ax, xyz: np.ndarray) -> None:
    xmin, ymin, zmin = xyz.min(axis=0)
    xmax, ymax, zmax = xyz.max(axis=0)
    xr, yr, zr = xmax - xmin, ymax - ymin, zmax - zmin
    r = max(xr, yr, zr)
    cx, cy, cz = (xmin + xmax) * 0.5, (ymin + ymax) * 0.5, (zmin + zmax) * 0.5
    ax.set_xlim(cx - 0.5 * r, cx + 0.5 * r)
    ax.set_ylim(cy - 0.5 * r, cy + 0.5 * r)
    ax.set_zlim(cz - 0.5 * r, cz + 0.5 * r)


def _infer_axis_spacing(values: np.ndarray) -> float:
    uniq = np.unique(np.round(values.astype(float), 12))
    if uniq.size <= 1:
        return 0.0
    diffs = np.diff(np.sort(uniq))
    diffs = diffs[diffs > 1e-15]
    if diffs.size == 0:
        return 0.0
    return float(np.median(diffs))


def _infer_block_size(centers: np.ndarray, scale: float) -> np.ndarray:
    spacing = np.array(
        [
            _infer_axis_spacing(centers[:, 0]),
            _infer_axis_spacing(centers[:, 1]),
            _infer_axis_spacing(centers[:, 2]),
        ],
        dtype=float,
    )
    span = np.ptp(centers, axis=0)
    fallback = np.where(span > 1e-15, span * 0.05, 1.0)
    size = np.where(spacing > 1e-15, spacing, fallback)
    return np.maximum(size * float(scale), 1e-15)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot 3D grid temperature from NTMC convergence output.")
    parser.add_argument("--input", required=True, help="Input TOML path")
    parser.add_argument("--results", default="", help="Results folder path (optional, auto-detect latest if omitted)")
    parser.add_argument("--tail", type=int, default=50, help="Average over last N rows")
    parser.add_argument("--out", default="", help="Output PNG path")
    parser.add_argument("--csv-out", default="", help="Output CSV path")
    parser.add_argument("--cmap", default="turbo", help="Matplotlib colormap")
    parser.add_argument("--block-scale", type=float, default=0.88, help="Block size scale relative to center spacing")
    parser.add_argument("--block-alpha", type=float, default=1.0, help="Block opacity in [0,1]")
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    cfg = _load_toml(input_path)
    results_dir = _resolve_results_dir(input_path, cfg, args.results or None)

    conv_path = results_dir / "convergence.txt"
    ts, temps, _rest, temp_cols = _load_convergence(conv_path)
    if temps.shape[1] == 0:
        raise ValueError("No grid temperature columns found in convergence.txt.")

    tail_n = max(1, min(args.tail, temps.shape[0]))
    temps_tail = temps[-tail_n:, :]
    t_mean = temps_tail.mean(axis=0)
    t_std = temps_tail.std(axis=0, ddof=1 if tail_n > 1 else 0)

    centers_csv = results_dir / "grid_centers.csv"
    if centers_csv.exists():
        centers = _load_centers_csv(centers_csv)
    else:
        centers = _build_centers_from_grid(cfg, input_path)
    if centers.shape[0] != t_mean.shape[0]:
        raise RuntimeError(
            f"Center count mismatch: centers={centers.shape[0]}, temperatures={t_mean.shape[0]}. "
            f"Expected centers file: {centers_csv}"
        )

    out_png = Path(args.out).resolve() if args.out else (results_dir / "temperature_3d_tail_mean.png")
    out_csv = Path(args.csv_out).resolve() if args.csv_out else (results_dir / "temperature_3d_tail_mean.csv")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    csv_data = np.column_stack([centers, t_mean, t_std])
    np.savetxt(
        out_csv,
        csv_data,
        delimiter=",",
        header="x,y,z,temp_mean,temp_std",
        comments="",
        fmt="%.10g",
    )

    fig = plt.figure(figsize=(10, 8), dpi=220)
    ax = fig.add_subplot(111, projection="3d")

    block_size = _infer_block_size(centers, args.block_scale)
    dx, dy, dz = float(block_size[0]), float(block_size[1]), float(block_size[2])
    x0 = centers[:, 0] - 0.5 * dx
    y0 = centers[:, 1] - 0.5 * dy
    z0 = centers[:, 2] - 0.5 * dz

    norm = colors.Normalize(vmin=float(np.min(t_mean)), vmax=float(np.max(t_mean)))
    cmap = plt.get_cmap(args.cmap)
    facecolors = cmap(norm(t_mean))
    facecolors[:, 3] = np.clip(float(args.block_alpha), 0.0, 1.0)

    ax.bar3d(
        x0,
        y0,
        z0,
        np.full_like(x0, dx),
        np.full_like(y0, dy),
        np.full_like(z0, dz),
        color=facecolors,
        shade=False,
        linewidth=0.0,
    )

    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.78, pad=0.08)
    cbar.set_label("Temperature (K)")

    _set_equal_axes(ax, centers)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(
        f"3D Temperature (tail mean, N={tail_n})\n"
        f"timestep: {int(ts[-tail_n])} ~ {int(ts[-1])}"
    )
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)

    print(f"[ok] results_dir: {results_dir}")
    print(f"[ok] convergence: {conv_path}")
    print(f"[ok] tail rows: {tail_n}")
    print(f"[ok] grids: {len(temp_cols)}")
    print(f"[ok] csv: {out_csv}")
    print(f"[ok] png: {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
