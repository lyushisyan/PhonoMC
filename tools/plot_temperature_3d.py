#!/usr/bin/env python3
"""Plot 3D temperature distribution from EPMC convergence output."""

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


def _load_convergence(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, List[str]]:
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
    temp_idx = [i for i, c in enumerate(cols[:heatflux_idx]) if c.startswith("T_")]
    if not temp_idx:
        # Backward compatibility: legacy files may not label T_* explicitly.
        start_idx = 2 if (len(cols) > 1 and cols[1] == "time_ps") else 1
        temp_idx = list(range(start_idx, heatflux_idx))
    temp_cols = [cols[i] for i in temp_idx]

    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    timesteps = data[:, 0]
    if "time_ps" in cols:
        tps_idx = cols.index("time_ps")
        time_ps = data[:, tps_idx]
    else:
        time_ps = np.full_like(timesteps, np.nan, dtype=float)
    temps = data[:, temp_idx] if temp_idx else np.zeros((data.shape[0], 0), dtype=float)
    rest = data[:, heatflux_idx:]
    return timesteps, time_ps, temps, rest, temp_cols


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


def _resolve_model_path(model: str, input_path: Path) -> Path:
    mesh_path = Path(model)
    if mesh_path.is_absolute():
        return mesh_path
    p1 = (Path.cwd() / mesh_path).resolve()
    if p1.exists():
        return p1
    return (input_path.parent / mesh_path).resolve()


def _build_box_surface_mesh(cfg: dict) -> Tuple[np.ndarray, np.ndarray]:
    geo = cfg.get("geometry", {})
    dims = geo.get("sizes")
    if not isinstance(dims, list) or len(dims) != 3:
        raise ValueError("Box model requires `sizes = [Lx, Ly, Lz]`.")
    lx = float(dims[0]) * 10.0
    ly = float(dims[1]) * 10.0
    lz = float(dims[2]) * 10.0
    vertices = np.array(
        [
            [0.0, 0.0, 0.0],
            [lx, 0.0, 0.0],
            [lx, ly, 0.0],
            [0.0, ly, 0.0],
            [0.0, 0.0, lz],
            [lx, 0.0, lz],
            [lx, ly, lz],
            [0.0, ly, lz],
        ],
        dtype=float,
    )
    # 12 triangles (2 per face)
    faces = np.array(
        [
            [0, 1, 2], [0, 2, 3],  # bottom z=0
            [4, 6, 5], [4, 7, 6],  # top z=lz
            [0, 4, 5], [0, 5, 1],  # y=0
            [1, 5, 6], [1, 6, 2],  # x=lx
            [2, 6, 7], [2, 7, 3],  # y=ly
            [3, 7, 4], [3, 4, 0],  # x=0
        ],
        dtype=int,
    )
    return vertices, faces


def _load_surface_mesh_from_config(cfg: dict, input_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    geo = cfg.get("geometry", {})
    model = geo.get("model")
    if not isinstance(model, str):
        raise ValueError("Missing `[geometry].model` in input TOML.")
    if model == "box":
        return _build_box_surface_mesh(cfg)

    mesh_path = _resolve_model_path(model, input_path)
    vertices, faces = _load_mesh(mesh_path)
    if mesh_path.suffix.lower() == ".stl":
        # STL geometry is in nm; solver converts to Angstrom internally.
        vertices = vertices * 10.0
    return vertices, faces


def _idw_interpolate_knn(
    sample_points: np.ndarray,
    sample_values: np.ndarray,
    query_points: np.ndarray,
    k: int,
    power: float,
) -> np.ndarray:
    if sample_points.shape[0] == 0:
        raise ValueError("No sample points for interpolation.")
    n = query_points.shape[0]
    out = np.empty(n, dtype=float)
    k = max(1, min(int(k), sample_points.shape[0]))
    p = max(0.1, float(power))
    eps = 1e-24
    chunk = 256
    for beg in range(0, n, chunk):
        end = min(n, beg + chunk)
        q = query_points[beg:end]  # (cq,3)
        d2 = np.sum((q[:, None, :] - sample_points[None, :, :]) ** 2, axis=2)  # (cq,ns)
        idx = np.argpartition(d2, kth=k - 1, axis=1)[:, :k]
        d2k = np.take_along_axis(d2, idx, axis=1)
        vk = sample_values[idx]
        near = np.argmin(d2k, axis=1)
        exact = d2k[np.arange(d2k.shape[0]), near] <= eps
        w = 1.0 / np.maximum(d2k, eps) ** (0.5 * p)
        val = np.sum(w * vk, axis=1) / np.sum(w, axis=1)
        if np.any(exact):
            val[exact] = vk[np.arange(vk.shape[0]), near][exact]
        out[beg:end] = val
    return out


def _geometry_bounds_from_config(cfg: dict, input_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    geo = cfg.get("geometry", {})
    model = geo.get("model")
    if model == "box":
        dims = geo.get("sizes")
        if not isinstance(dims, list) or len(dims) != 3:
            raise ValueError("Box model requires `sizes = [Lx, Ly, Lz]`.")
        bmin = np.array([0.0, 0.0, 0.0], dtype=float)
        bmax = np.array([float(dims[0]), float(dims[1]), float(dims[2])], dtype=float) * 10.0
        return bmin, bmax
    vertices, _faces = _load_surface_mesh_from_config(cfg, input_path)
    return vertices.min(axis=0), vertices.max(axis=0)


def _compute_heat_source_grid_power_density(
    cfg: dict,
    input_path: Path,
    centers: np.ndarray,
) -> Tuple[np.ndarray, dict]:
    out = np.zeros(centers.shape[0], dtype=float)
    hs = cfg.get("heat_source", {})
    if not isinstance(hs, dict):
        return out, {"enabled": False, "reason": "no_heat_source_section"}

    def _vec3(v: object) -> np.ndarray | None:
        if not isinstance(v, list) or len(v) != 3:
            return None
        try:
            arr = np.array([float(v[0]), float(v[1]), float(v[2])], dtype=float)
        except Exception:
            return None
        if not np.all(np.isfinite(arr)):
            return None
        return arr

    profile = str(hs.get("profile", "uniform")).strip().lower()
    if profile not in ("uniform", "gaussian"):
        profile = "uniform"

    enabled = bool(hs.get("enabled", False))
    min_rel = _vec3(hs.get("min"))
    max_rel = _vec3(hs.get("max"))
    center_rel = _vec3(hs.get("center"))
    sigma_rel = _vec3(hs.get("sigma"))
    try:
        power_density = float(hs.get("power_density", 0.0))
    except Exception:
        power_density = 0.0
    if not enabled and abs(power_density) > 0.0:
        if profile == "uniform":
            enabled = (min_rel is not None) and (max_rel is not None)
        else:
            enabled = (center_rel is not None) and (sigma_rel is not None)
    if not enabled:
        return out, {"enabled": False, "reason": "heat_source_disabled"}
    if not np.isfinite(power_density) or abs(power_density) <= 0.0:
        return out, {"enabled": False, "reason": "heat_source_power_density_zero"}

    bmin, bmax = _geometry_bounds_from_config(cfg, input_path)
    ext = bmax - bmin

    power_density_role = "region_value"
    if profile == "uniform":
        if min_rel is None or max_rel is None:
            return out, {"enabled": False, "reason": "heat_source_min_max_invalid"}
        hs_min = bmin + min_rel * ext
        hs_max = bmin + max_rel * ext
        hs_lo = np.minimum(hs_min, hs_max)
        hs_hi = np.maximum(hs_min, hs_max)
        inside = (
            (centers[:, 0] >= hs_lo[0]) & (centers[:, 0] <= hs_hi[0]) &
            (centers[:, 1] >= hs_lo[1]) & (centers[:, 1] <= hs_hi[1]) &
            (centers[:, 2] >= hs_lo[2]) & (centers[:, 2] <= hs_hi[2])
        )
        selected = int(np.count_nonzero(inside))
        if selected <= 0:
            return out, {"enabled": False, "reason": "heat_source_selects_no_grid"}
        out[inside] = power_density
    else:
        if center_rel is None or sigma_rel is None:
            return out, {"enabled": False, "reason": "heat_source_center_sigma_invalid"}
        if centers.shape[0] <= 0:
            return out, {"enabled": False, "reason": "heat_source_selects_no_grid"}
        center = bmin + center_rel * ext
        sigma = np.ones(3, dtype=float)
        gaussian_axis_enabled = np.array([True, True, True], dtype=bool)
        for k in range(3):
            sv = float(sigma_rel[k])
            if np.isfinite(sv) and sv > 0.0:
                sigma[k] = max(1e-12, sv * ext[k])
                gaussian_axis_enabled[k] = True
            else:
                sigma[k] = 1.0
                gaussian_axis_enabled[k] = False
        d = (centers - center[None, :]) / sigma[None, :]
        for k in range(3):
            if not gaussian_axis_enabled[k]:
                d[:, k] = 0.0
        r2 = np.sum(d * d, axis=1)
        out = power_density * np.exp(-0.5 * r2)
        selected = int(centers.shape[0])
        power_density_role = "peak_value"

    info = {
        "enabled": True,
        "profile": profile,
        "power_density": power_density,
        "power_density_role": power_density_role,
        "selected_grids": selected,
    }
    return out, info


def _render_3d_scalar_field(
    cfg: dict,
    input_path: Path,
    centers: np.ndarray,
    values: np.ndarray,
    out_png: Path,
    cmap_name: str,
    cbar_label: str,
    title_template: str,
    view3d_render: str,
    block_scale: float,
    block_alpha: float,
    surface_knn: int,
    surface_idw_power: float,
    surface_alpha: float,
    block_positive_only: bool = False,
) -> str:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(10, 8), dpi=220)
    ax = fig.add_subplot(111, projection="3d")

    vmin = float(np.nanmin(values)) if values.size else 0.0
    vmax = float(np.nanmax(values)) if values.size else 0.0
    if abs(vmax - vmin) <= 1e-12:
        vmax = vmin + 1e-12
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)

    _ = (cfg, input_path, view3d_render, surface_knn, surface_idw_power, surface_alpha)
    if block_positive_only:
        mask = values > 0.0
        if np.any(mask):
            pcenters = centers[mask]
            pvals = values[mask]
        else:
            pcenters = centers
            pvals = values
    else:
        pcenters = centers
        pvals = values

    block_size = _infer_block_size(centers, block_scale)
    dx, dy, dz = float(block_size[0]), float(block_size[1]), float(block_size[2])
    x0 = pcenters[:, 0] - 0.5 * dx
    y0 = pcenters[:, 1] - 0.5 * dy
    z0 = pcenters[:, 2] - 0.5 * dz
    facecolors = cmap(norm(pvals))
    facecolors[:, 3] = np.clip(float(block_alpha), 0.0, 1.0)

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
    _set_equal_axes(ax, centers)
    render_mode_3d = "block"

    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.78, pad=0.08)
    cbar.set_label(cbar_label)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(title_template.format(mode=render_mode_3d))
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    return render_mode_3d


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

    mesh_path = _resolve_model_path(model, input_path)
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
    header = path.read_text(encoding="utf-8").splitlines()[0].strip().lower()
    is_nm = ("x_nm" in header) and ("y_nm" in header) and ("z_nm" in header)
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 4:
        raise ValueError(f"Invalid centers csv format: {path}")
    centers = data[:, 1:4]
    if is_nm:
        # Solver internals and mesh overlays use Angstrom.
        centers = centers * 10.0
    return centers


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


def _nearest_slice_mask(values: np.ndarray, target: float) -> Tuple[np.ndarray, float]:
    vals = values.astype(float)
    if vals.size == 0:
        raise ValueError("Cannot build slice on empty coordinates.")
    uniq = np.unique(np.round(vals, 12))
    nearest = float(uniq[int(np.argmin(np.abs(uniq - target)))])
    spacing = _infer_axis_spacing(vals)
    tol = max(1e-12, 0.51 * spacing) if spacing > 0.0 else 1e-12
    mask = np.abs(vals - nearest) <= tol
    if not np.any(mask):
        idx = int(np.argmin(np.abs(vals - nearest)))
        mask = np.zeros(vals.shape, dtype=bool)
        mask[idx] = True
    return mask, nearest


def _plot_slice_distribution(
    points2d: np.ndarray,
    temps: np.ndarray,
    out_png: Path,
    out_csv: Path,
    title: str,
    xlabel: str,
    ylabel: str,
    cmap_name: str,
    render_mode: str,
    smooth_levels: int,
    smooth_grid_scale: float,
) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    csv_data = np.column_stack([points2d, temps])
    np.savetxt(out_csv, csv_data, delimiter=",", header=f"{xlabel},{ylabel},temp_mean", comments="", fmt="%.10g")

    # Build a regular 2D grid for block mode (no visual gaps between neighboring cells).
    x_centers = np.unique(np.round(points2d[:, 0].astype(float), 12))
    y_centers = np.unique(np.round(points2d[:, 1].astype(float), 12))
    x_centers.sort()
    y_centers.sort()
    x_to_i = {float(v): i for i, v in enumerate(x_centers)}
    y_to_i = {float(v): i for i, v in enumerate(y_centers)}

    grid = np.full((y_centers.size, x_centers.size), np.nan, dtype=float)
    for (xv, yv), tv in zip(points2d, temps):
        ix = x_to_i[float(round(float(xv), 12))]
        iy = y_to_i[float(round(float(yv), 12))]
        grid[iy, ix] = float(tv)

    def _centers_to_edges(c: np.ndarray) -> np.ndarray:
        if c.size == 1:
            return np.array([c[0] - 0.5, c[0] + 0.5], dtype=float)
        mids = 0.5 * (c[:-1] + c[1:])
        left = c[0] - 0.5 * (c[1] - c[0])
        right = c[-1] + 0.5 * (c[-1] - c[-2])
        return np.concatenate(([left], mids, [right]))

    x_edges = _centers_to_edges(x_centers)
    y_edges = _centers_to_edges(y_centers)

    fig, ax = plt.subplots(figsize=(7.4, 6.0), dpi=220)
    vmin = float(np.nanmin(temps))
    vmax = float(np.nanmax(temps))
    same_level = abs(vmax - vmin) <= 1e-12
    norm = colors.Normalize(vmin=vmin, vmax=vmax if not same_level else vmin + 1e-12)
    cmap = plt.get_cmap(cmap_name)
    grid_masked = np.ma.masked_invalid(grid)
    smooth_ok = (render_mode == "smooth") and (not same_level)
    if smooth_ok:
        valid_nodes = int(np.count_nonzero(~grid_masked.mask))
        smooth_ok = valid_nodes >= 3
    if smooth_ok:
        # Smooth contour on regular grid with mask kept: void/outside cells stay unpainted.
        Xc, Yc = np.meshgrid(x_centers, y_centers)
        levels = max(16, int(round(float(smooth_levels) * max(0.5, min(2.0, float(smooth_grid_scale) / 3.0)))))
        try:
            pm = ax.contourf(
                Xc,
                Yc,
                grid_masked,
                levels=levels,
                cmap=cmap,
                norm=norm,
                corner_mask=False,
                antialiased=True,
            )
        except Exception:
            smooth_ok = False
    if not smooth_ok:
        pm = ax.pcolormesh(
            x_edges,
            y_edges,
            grid_masked,
            cmap=cmap,
            norm=norm,
            shading="flat",
            edgecolors="none",
            linewidth=0.0,
            antialiased=False,
        )

    cbar = fig.colorbar(pm, ax=ax, shrink=0.86, pad=0.02)
    cbar.set_label("Temperature (K)")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def _plot_topk_convergence(
    timesteps: np.ndarray,
    time_ps: np.ndarray,
    temps: np.ndarray,
    topk: int,
    out_png: Path,
    out_csv: Path,
) -> int:
    if temps.size == 0:
        raise ValueError("No temperature data available for top-k convergence plot.")
    k = max(1, min(int(topk), temps.shape[1]))
    topk_vals = np.partition(temps, temps.shape[1] - k, axis=1)[:, -k:]
    topk_vals.sort(axis=1)
    topk_vals = topk_vals[:, ::-1]  # top1 >= top2 >= ...

    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    if np.all(np.isfinite(time_ps)):
        x = time_ps
        xlabel = "time (ps)"
    else:
        x = timesteps
        xlabel = "timestep"

    fig, ax = plt.subplots(figsize=(8.2, 5.2), dpi=220)
    for i in range(k):
        ax.plot(x, topk_vals[:, i], lw=1.5, label=f"Top{i + 1}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Temperature (K)")
    ax.set_title(f"Top-{k} Highest Grid Temperatures vs {xlabel}")
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.legend(frameon=False, ncol=1)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)

    cols = [timesteps, time_ps] + [topk_vals[:, i] for i in range(k)]
    header = ["timestep", "time_ps"] + [f"top{i + 1}" for i in range(k)]
    np.savetxt(
        out_csv,
        np.column_stack(cols),
        delimiter=",",
        header=",".join(header),
        comments="",
        fmt="%.10g",
    )
    return k


def _plot_xavg_profile_vs_z(
    xz_points: np.ndarray,
    temps: np.ndarray,
    out_png: Path,
    out_csv: Path,
    title: str,
) -> None:
    if xz_points.shape[0] == 0 or temps.size == 0:
        raise ValueError("No slice points available for x-averaged z profile.")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    z = xz_points[:, 1].astype(float)
    z_key = np.round(z, 12)
    z_unique = np.unique(z_key)
    z_unique.sort()

    z_vals = np.empty(z_unique.size, dtype=float)
    t_mean = np.empty(z_unique.size, dtype=float)
    nx_count = np.empty(z_unique.size, dtype=int)
    for i, zk in enumerate(z_unique):
        mask = z_key == zk
        z_vals[i] = float(np.mean(z[mask]))
        t_mean[i] = float(np.mean(temps[mask]))
        nx_count[i] = int(np.count_nonzero(mask))

    order = np.argsort(z_vals)
    z_vals = z_vals[order]
    t_mean = t_mean[order]
    nx_count = nx_count[order]

    np.savetxt(
        out_csv,
        np.column_stack([z_vals, t_mean, nx_count]),
        delimiter=",",
        header="z_nm,temp_xavg,n_x_points",
        comments="",
        fmt=["%.10g", "%.10g", "%d"],
    )

    fig, ax = plt.subplots(figsize=(7.4, 5.2), dpi=220)
    ax.plot(z_vals, t_mean, lw=1.8, marker="o", markersize=2.8)
    ax.set_xlabel("z_nm")
    ax.set_ylabel("Temperature (K)")
    ax.set_title(title)
    ax.grid(alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot 3D grid temperature from EPMC convergence output.")
    parser.add_argument("--input", required=True, help="Input TOML path")
    parser.add_argument("--results", default="", help="Results folder path (optional, auto-detect latest if omitted)")
    parser.add_argument("--tail", type=int, default=50, help="Average over last N rows")
    parser.add_argument("--out", default="", help="Output PNG path")
    parser.add_argument("--csv-out", default="", help="Output CSV path")
    parser.add_argument("--cmap", default="turbo", help="Matplotlib colormap")
    parser.add_argument(
        "--view3d-render",
        choices=["surface", "block"],
        default="block",
        help="3D render mode (currently forced to block voxels)",
    )
    parser.add_argument("--surface-knn", type=int, default=12, help="KNN count for surface temperature interpolation")
    parser.add_argument("--surface-idw-power", type=float, default=2.0, help="IDW power for surface interpolation")
    parser.add_argument("--surface-alpha", type=float, default=1.0, help="Surface opacity in [0,1]")
    parser.add_argument("--block-scale", type=float, default=0.88, help="Block size scale relative to center spacing")
    parser.add_argument("--block-alpha", type=float, default=1.0, help="Block opacity in [0,1]")
    parser.add_argument(
        "--slice-render",
        choices=["smooth", "block"],
        default="smooth",
        help="2D slice render mode: smooth contour or block map",
    )
    parser.add_argument("--slice-smooth-levels", type=int, default=96, help="Contour levels for smooth slice rendering")
    parser.add_argument(
        "--slice-smooth-grid-scale",
        type=float,
        default=3.0,
        help="Smooth contour density scale (multiplies contour levels)",
    )
    parser.add_argument("--x-slice-rel", type=float, default=0.5, help="Relative x location for YZ slice (0~1)")
    parser.add_argument("--y-slice-rel", type=float, default=0.6, help="Relative y location for XZ slice (0~1)")
    parser.add_argument("--topk", type=int, default=5, help="Number of highest-temperature convergence lines")
    args = parser.parse_args()
    if args.view3d_render != "block":
        print("[info] --view3d-render=surface is ignored; using block voxel rendering.")

    input_path = Path(args.input).resolve()
    cfg = _load_toml(input_path)
    results_dir = _resolve_results_dir(input_path, cfg, args.results or None)

    conv_path = results_dir / "convergence.txt"
    ts, tps, temps, _rest, temp_cols = _load_convergence(conv_path)
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

    render_mode_3d = _render_3d_scalar_field(
        cfg=cfg,
        input_path=input_path,
        centers=centers,
        values=t_mean,
        out_png=out_png,
        cmap_name=args.cmap,
        cbar_label="Temperature (K)",
        title_template=(
            f"3D Temperature ({{mode}}, tail mean, N={tail_n})\n"
            f"timestep: {int(ts[-tail_n])} ~ {int(ts[-1])}"
        ),
        view3d_render=args.view3d_render,
        block_scale=args.block_scale,
        block_alpha=args.block_alpha,
        surface_knn=args.surface_knn,
        surface_idw_power=args.surface_idw_power,
        surface_alpha=args.surface_alpha,
    )

    source_wm3, source_info = _compute_heat_source_grid_power_density(cfg, input_path, centers)
    out_source_csv = results_dir / "heat_source_distribution.csv"
    out_source_png = results_dir / "heat_source_distribution_3d.png"
    source_render_mode = "disabled"
    if source_info.get("enabled", False):
        source_csv = np.column_stack([centers, source_wm3])
        np.savetxt(
            out_source_csv,
            source_csv,
            delimiter=",",
            header="x,y,z,source_wm3",
            comments="",
            fmt="%.10g",
        )
        profile = str(source_info.get("profile", "uniform"))
        selected = int(source_info.get("selected_grids", 0))
        power0 = float(source_info.get("power_density", 0.0))
        power_role = str(source_info.get("power_density_role", "value"))
        source_render_mode = _render_3d_scalar_field(
            cfg=cfg,
            input_path=input_path,
            centers=centers,
            values=source_wm3,
            out_png=out_source_png,
            cmap_name=args.cmap,
            cbar_label="Heat source (W/m^3)",
            title_template=(
                f"Heat Source Distribution ({{mode}})\n"
                f"profile={profile}, power_density({power_role})={power0:.6g} W/m^3, selected_grids={selected}"
            ),
            view3d_render=args.view3d_render,
            block_scale=args.block_scale,
            block_alpha=args.block_alpha,
            surface_knn=args.surface_knn,
            surface_idw_power=args.surface_idw_power,
            surface_alpha=args.surface_alpha,
            block_positive_only=True,
        )
    else:
        reason = source_info.get("reason", "unknown")
        print(f"[info] heat source plot skipped: {reason}")

    centers_nm = centers / 10.0
    xmin, ymin, _zmin = centers_nm.min(axis=0)
    xmax, ymax, _zmax = centers_nm.max(axis=0)

    x_rel = float(np.clip(args.x_slice_rel, 0.0, 1.0))
    y_rel = float(np.clip(args.y_slice_rel, 0.0, 1.0))
    x_target = xmin + x_rel * (xmax - xmin)
    y_target = ymin + y_rel * (ymax - ymin)

    x_mask, x_plane = _nearest_slice_mask(centers_nm[:, 0], x_target)
    x_slice_points = centers_nm[x_mask][:, [1, 2]]  # y,z
    x_slice_t = t_mean[x_mask]
    out_xslice_png = results_dir / f"temperature_slice_xrel{x_rel:.3f}_yz.png"
    out_xslice_csv = results_dir / f"temperature_slice_xrel{x_rel:.3f}_yz.csv"
    _plot_slice_distribution(
        points2d=x_slice_points,
        temps=x_slice_t,
        out_png=out_xslice_png,
        out_csv=out_xslice_csv,
        title=f"YZ Slice @ x={x_plane:.3f} nm (x_rel={x_rel:.3f})",
        xlabel="y_nm",
        ylabel="z_nm",
        cmap_name=args.cmap,
        render_mode=args.slice_render,
        smooth_levels=args.slice_smooth_levels,
        smooth_grid_scale=args.slice_smooth_grid_scale,
    )

    y_mask, y_plane = _nearest_slice_mask(centers_nm[:, 1], y_target)
    y_slice_points = centers_nm[y_mask][:, [0, 2]]  # x,z
    y_slice_t = t_mean[y_mask]
    out_yslice_png = results_dir / f"temperature_slice_yrel{y_rel:.3f}_xz.png"
    out_yslice_csv = results_dir / f"temperature_slice_yrel{y_rel:.3f}_xz.csv"
    _plot_slice_distribution(
        points2d=y_slice_points,
        temps=y_slice_t,
        out_png=out_yslice_png,
        out_csv=out_yslice_csv,
        title=f"XZ Slice @ y={y_plane:.3f} nm (y_rel={y_rel:.3f})",
        xlabel="x_nm",
        ylabel="z_nm",
        cmap_name=args.cmap,
        render_mode=args.slice_render,
        smooth_levels=args.slice_smooth_levels,
        smooth_grid_scale=args.slice_smooth_grid_scale,
    )
    out_yzprof_png = results_dir / f"temperature_profile_yrel{y_rel:.3f}_xavg_vs_z.png"
    out_yzprof_csv = results_dir / f"temperature_profile_yrel{y_rel:.3f}_xavg_vs_z.csv"
    _plot_xavg_profile_vs_z(
        xz_points=y_slice_points,
        temps=y_slice_t,
        out_png=out_yzprof_png,
        out_csv=out_yzprof_csv,
        title=f"X-averaged Temperature vs z @ y={y_plane:.3f} nm (y_rel={y_rel:.3f})",
    )

    out_topk_png = results_dir / "temperature_topk_convergence.png"
    out_topk_csv = results_dir / "temperature_topk_convergence.csv"
    k_used = _plot_topk_convergence(
        timesteps=ts,
        time_ps=tps,
        temps=temps,
        topk=args.topk,
        out_png=out_topk_png,
        out_csv=out_topk_csv,
    )

    print(f"[ok] results_dir: {results_dir}")
    print(f"[ok] convergence: {conv_path}")
    print(f"[ok] tail rows: {tail_n}")
    print(f"[ok] grids: {len(temp_cols)}")
    print(f"[ok] 3d_render: {render_mode_3d}")
    print(f"[ok] csv: {out_csv}")
    print(f"[ok] png: {out_png}")
    if source_info.get("enabled", False):
        print(f"[ok] heat_source_3d_render: {source_render_mode}")
        print(f"[ok] heat_source csv: {out_source_csv}")
        print(f"[ok] heat_source png: {out_source_png}")
    print(f"[ok] x-slice png: {out_xslice_png}")
    print(f"[ok] x-slice csv: {out_xslice_csv}")
    print(f"[ok] y-slice png: {out_yslice_png}")
    print(f"[ok] y-slice csv: {out_yslice_csv}")
    print(f"[ok] y-slice xavg-vs-z png: {out_yzprof_png}")
    print(f"[ok] y-slice xavg-vs-z csv: {out_yzprof_csv}")
    print(f"[ok] x-slice plane: x={x_plane:.6g} nm, points={int(np.count_nonzero(x_mask))}")
    print(f"[ok] y-slice plane: y={y_plane:.6g} nm, points={int(np.count_nonzero(y_mask))}")
    print(f"[ok] top{k_used} png: {out_topk_png}")
    print(f"[ok] top{k_used} csv: {out_topk_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
