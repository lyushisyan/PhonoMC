#!/usr/bin/env python3
"""Generate a simple FinFET outline STL (no top split)."""

from __future__ import annotations

import argparse
import math
import struct
from pathlib import Path
from typing import List, Sequence, Tuple

Vec3 = Tuple[float, float, float]
Face = Tuple[int, int, int]


def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(v: Vec3) -> float:
    return math.sqrt(_dot(v, v))


def _normalize(v: Vec3) -> Vec3:
    n = _norm(v)
    if n <= 0.0:
        return (0.0, 0.0, 0.0)
    return (v[0] / n, v[1] / n, v[2] / n)


def _triangle_normal(v0: Vec3, v1: Vec3, v2: Vec3) -> Vec3:
    return _normalize(_cross(_sub(v1, v0), _sub(v2, v0)))


def _centroid(points: Sequence[Vec3]) -> Vec3:
    n = float(len(points))
    return (
        sum(p[0] for p in points) / n,
        sum(p[1] for p in points) / n,
        sum(p[2] for p in points) / n,
    )


def _fix_face_orientation(vertices: Sequence[Vec3], faces: Sequence[Face]) -> List[Face]:
    c = _centroid(vertices)
    fixed: List[Face] = []
    for i0, i1, i2 in faces:
        v0 = vertices[i0]
        v1 = vertices[i1]
        v2 = vertices[i2]
        n = _triangle_normal(v0, v1, v2)
        fc = ((v0[0] + v1[0] + v2[0]) / 3.0, (v0[1] + v1[1] + v2[1]) / 3.0, (v0[2] + v1[2] + v2[2]) / 3.0)
        # Outward if normal points away from global center.
        if _dot(n, _sub(fc, c)) < 0.0:
            fixed.append((i0, i2, i1))
        else:
            fixed.append((i0, i1, i2))
    return fixed


def _write_binary_stl(path: Path, vertices: Sequence[Vec3], faces: Sequence[Face], name: str) -> None:
    header = name.encode("ascii", errors="ignore")[:80].ljust(80, b" ")
    with path.open("wb") as f:
        f.write(header)
        f.write(struct.pack("<I", len(faces)))
        for i0, i1, i2 in faces:
            v0 = vertices[i0]
            v1 = vertices[i1]
            v2 = vertices[i2]
            n = _triangle_normal(v0, v1, v2)
            f.write(
                struct.pack(
                    "<12fH",
                    float(n[0]),
                    float(n[1]),
                    float(n[2]),
                    float(v0[0]),
                    float(v0[1]),
                    float(v0[2]),
                    float(v1[0]),
                    float(v1[1]),
                    float(v1[2]),
                    float(v2[0]),
                    float(v2[1]),
                    float(v2[2]),
                    0,
                )
            )


def _build_simple_finfet(
    base_width: float,
    base_height: float,
    stem_width: float,
    stem_height: float,
    thickness_y: float,
) -> Tuple[List[Vec3], List[Face]]:
    if base_width <= 0.0 or base_height <= 0.0 or stem_width <= 0.0 or stem_height <= 0.0 or thickness_y <= 0.0:
        raise ValueError("All dimensions must be positive.")
    if stem_width > base_width:
        raise ValueError("stem_width must be <= base_width.")

    x_b_half = 0.5 * base_width
    x_s_half = 0.5 * stem_width
    z_base = base_height
    z_top = base_height + stem_height

    # Front face (y=0), 8 points.
    vf: List[Vec3] = [
        (-x_b_half, 0.0, 0.0),      # 0
        (x_b_half, 0.0, 0.0),       # 1
        (x_b_half, 0.0, z_base),    # 2
        (x_s_half, 0.0, z_base),    # 3
        (x_s_half, 0.0, z_top),     # 4
        (-x_s_half, 0.0, z_top),    # 5
        (-x_s_half, 0.0, z_base),   # 6
        (-x_b_half, 0.0, z_base),   # 7
    ]

    vb = [(x, thickness_y, z) for (x, _, z) in vf]  # 8..15
    vertices: List[Vec3] = vf + vb

    faces: List[Face] = []

    # Front / back
    faces.extend([(0, 1, 2), (0, 2, 7)])
    faces.extend([(6, 3, 4), (6, 4, 5)])
    faces.extend([(8, 15, 9), (9, 15, 10)])
    faces.extend([(14, 13, 11), (11, 13, 12)])

    # Side walls / horizontal surfaces (no top split).
    faces.extend([(0, 8, 1), (1, 8, 9)])          # bottom
    faces.extend([(0, 7, 8), (8, 7, 15)])         # left base wall
    faces.extend([(1, 9, 2), (9, 10, 2)])         # right base wall
    faces.extend([(7, 6, 15), (15, 6, 14)])       # base top left flat
    faces.extend([(2, 10, 3), (10, 11, 3)])       # base top right flat
    faces.extend([(6, 5, 14), (14, 5, 13)])       # stem left wall
    faces.extend([(3, 11, 4), (11, 12, 4)])       # stem right wall
    faces.extend([(5, 4, 12), (5, 12, 13)])       # stem top (single flat top)

    faces = _fix_face_orientation(vertices, faces)
    return vertices, faces


def _save_preview_png(vertices: Sequence[Vec3], faces: Sequence[Face], out_png: Path) -> None:
    try:
        import numpy as np
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    except Exception as exc:  # pragma: no cover
        print(f"[warn] Preview skipped ({exc})")
        return

    v = np.array(vertices, dtype=float)
    tri = np.array(faces, dtype=int)
    polys = v[tri]

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    coll = Poly3DCollection(polys, alpha=0.55)
    coll.set_facecolor([0.2, 0.65, 1.0])
    coll.set_edgecolor("k")
    ax.add_collection3d(coll)

    xmin, xmax = float(v[:, 0].min()), float(v[:, 0].max())
    ymin, ymax = float(v[:, 1].min()), float(v[:, 1].max())
    zmin, zmax = float(v[:, 2].min()), float(v[:, 2].max())
    xr, yr, zr = xmax - xmin, ymax - ymin, zmax - zmin
    r = max(xr, yr, zr)
    cx, cy, cz = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax)
    ax.set_xlim(cx - 0.5 * r, cx + 0.5 * r)
    ax.set_ylim(cy - 0.5 * r, cy + 0.5 * r)
    ax.set_zlim(cz - 0.5 * r, cz + 0.5 * r)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Simple FinFET (no top split)")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=250)
    plt.close(fig)
    print(f"[ok] Preview saved: {out_png.resolve()}")


def _make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Generate a simple FinFET STL with unsplit top surface.")
    p.add_argument("--output", default="Model/finfet.stl", help="Output STL path")
    p.add_argument("--name", default="finfet_simple", help="STL header name")

    # Keep the same meaning as your reference code (default values in meters).
    p.add_argument("--base-width", type=float, default=22e-9, help="Base width")
    p.add_argument("--base-height", type=float, default=10e-9, help="Base height")
    p.add_argument("--stem-width", type=float, default=8e-9, help="Stem width")
    p.add_argument("--stem-height", type=float, default=20e-9, help="Stem height")
    p.add_argument("--thickness-y", type=float, default=36e-9, help="Extrusion thickness along Y")

    p.add_argument("--scale", type=float, default=1.0, help="Uniform scale factor applied to all dimensions")
    p.add_argument("--preview", action="store_true", help="Also save a PNG preview")
    p.add_argument("--preview-output", default="Model/finfet_preview.png", help="Preview PNG path")
    return p


def main() -> int:
    args = _make_parser().parse_args()

    s = float(args.scale)
    vertices, faces = _build_simple_finfet(
        base_width=args.base_width * s,
        base_height=args.base_height * s,
        stem_width=args.stem_width * s,
        stem_height=args.stem_height * s,
        thickness_y=args.thickness_y * s,
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    _write_binary_stl(out, vertices, faces, args.name)

    print(f"[ok] STL written: {out.resolve()}")
    print(f"[ok] Total vertices: {len(vertices)}")
    print(f"[ok] Total faces: {len(faces)}")
    print("[ok] Top surface: unsplit (single flat stem top)")

    if args.preview:
        _save_preview_png(vertices, faces, Path(args.preview_output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

