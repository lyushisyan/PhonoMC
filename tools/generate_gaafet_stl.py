#!/usr/bin/env python3
"""Generate a simple double-channel GAAFET-like STL.

Default geometry (nm):
- Left pad:     10 x 20 x 30 (L x W x H)
- Channels: 2 x 32 x 10 x 6  (L x W x H), uniformly distributed along Z
- Right pad:    10 x 20 x 30 (L x W x H)

Coordinates:
- X: length direction
- Y: width direction (centered at 0)
- Z: height direction (bottom at 0)

For the default total height 30 nm and two 6 nm channels, the vertical gap is:
    (30 - 2*6) / (2 + 1) = 6 nm
so the two channel Z-ranges are [6, 12] nm and [18, 24] nm.
"""

from __future__ import annotations

import argparse
import math
import struct
import zlib
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

Vec3 = Tuple[float, float, float]
Face = Tuple[int, int, int]
Box = Tuple[float, float, float, float, float, float]  # x0, x1, y0, y1, z0, z1


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


def _project_iso(p: Vec3) -> Tuple[float, float]:
    x, y, z = p
    u = x - 0.6 * y
    v = z + 0.35 * y
    return u, v


def _draw_line_rgb(img: bytearray, w: int, h: int, x0: int, y0: int, x1: int, y1: int, rgb: Tuple[int, int, int]) -> None:
    dx = abs(x1 - x0)
    dy = -abs(y1 - y0)
    sx = 1 if x0 < x1 else -1
    sy = 1 if y0 < y1 else -1
    err = dx + dy
    x, y = x0, y0

    while True:
        if 0 <= x < w and 0 <= y < h:
            idx = (y * w + x) * 3
            img[idx] = rgb[0]
            img[idx + 1] = rgb[1]
            img[idx + 2] = rgb[2]
        if x == x1 and y == y1:
            break
        e2 = 2 * err
        if e2 >= dy:
            err += dy
            x += sx
        if e2 <= dx:
            err += dx
            y += sy


def _png_chunk(tag: bytes, payload: bytes) -> bytes:
    return struct.pack(">I", len(payload)) + tag + payload + struct.pack(">I", zlib.crc32(tag + payload) & 0xFFFFFFFF)


def _write_png_rgb(path: Path, w: int, h: int, rgb_bytes: bytes) -> None:
    sig = b"\x89PNG\r\n\x1a\n"
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)
    # Filter type 0 for each row.
    raw = bytearray()
    stride = w * 3
    for row in range(h):
        raw.append(0)
        s = row * stride
        raw.extend(rgb_bytes[s : s + stride])
    idat = zlib.compress(bytes(raw), level=6)

    data = bytearray()
    data.extend(sig)
    data.extend(_png_chunk(b"IHDR", ihdr))
    data.extend(_png_chunk(b"IDAT", idat))
    data.extend(_png_chunk(b"IEND", b""))

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(bytes(data))


def _save_preview_png_fallback(vertices: Sequence[Vec3], faces: Sequence[Face], out_png: Path) -> None:
    w, h = 1200, 900
    margin = 60.0
    bg = (248, 251, 255)
    fg = (20, 26, 35)

    uv = [_project_iso(p) for p in vertices]
    umin = min(u for u, _ in uv)
    umax = max(u for u, _ in uv)
    vmin = min(v for _, v in uv)
    vmax = max(v for _, v in uv)
    ur = max(umax - umin, 1e-9)
    vr = max(vmax - vmin, 1e-9)
    scale = min((w - 2.0 * margin) / ur, (h - 2.0 * margin) / vr)

    projected: List[Tuple[int, int]] = []
    for u, v in uv:
        x = int(round((u - umin) * scale + margin))
        y = int(round(h - ((v - vmin) * scale + margin)))
        projected.append((x, y))

    img = bytearray([0] * (w * h * 3))
    for i in range(0, len(img), 3):
        img[i] = bg[0]
        img[i + 1] = bg[1]
        img[i + 2] = bg[2]

    edges: Set[Tuple[int, int]] = set()
    for a, b, c in faces:
        e1 = (a, b) if a < b else (b, a)
        e2 = (b, c) if b < c else (c, b)
        e3 = (c, a) if c < a else (a, c)
        edges.add(e1)
        edges.add(e2)
        edges.add(e3)

    for i0, i1 in edges:
        x0, y0 = projected[i0]
        x1, y1 = projected[i1]
        _draw_line_rgb(img, w, h, x0, y0, x1, y1, fg)

    _write_png_rgb(out_png, w, h, bytes(img))
    print(f"[ok] Preview saved (fallback): {out_png.resolve()}")


def _save_preview_png(vertices: Sequence[Vec3], faces: Sequence[Face], out_png: Path) -> None:
    try:
        import numpy as np
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    except Exception as exc:  # pragma: no cover
        print(f"[warn] matplotlib preview unavailable ({exc}), using fallback renderer")
        _save_preview_png_fallback(vertices, faces, out_png)
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
    ax.set_title("Double-channel GAAFET-like geometry")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=250)
    plt.close(fig)
    print(f"[ok] Preview saved: {out_png.resolve()}")


def _point_in_box(p: Vec3, b: Box) -> bool:
    x, y, z = p
    x0, x1, y0, y1, z0, z1 = b
    return x0 <= x <= x1 and y0 <= y <= y1 and z0 <= z <= z1


def _add_vertex(vmap: Dict[Vec3, int], vertices: List[Vec3], p: Vec3) -> int:
    idx = vmap.get(p)
    if idx is None:
        idx = len(vertices)
        vmap[p] = idx
        vertices.append(p)
    return idx


def _add_oriented_quad(
    vertices: List[Vec3],
    faces: List[Face],
    vmap: Dict[Vec3, int],
    p0: Vec3,
    p1: Vec3,
    p2: Vec3,
    p3: Vec3,
    expected_normal: Vec3,
) -> None:
    i0 = _add_vertex(vmap, vertices, p0)
    i1 = _add_vertex(vmap, vertices, p1)
    i2 = _add_vertex(vmap, vertices, p2)
    i3 = _add_vertex(vmap, vertices, p3)

    n = _triangle_normal(vertices[i0], vertices[i1], vertices[i2])
    if _dot(n, expected_normal) < 0.0:
        faces.append((i0, i2, i1))
        faces.append((i0, i3, i2))
    else:
        faces.append((i0, i1, i2))
        faces.append((i0, i2, i3))


def _merge_rectangles(mask: Sequence[Sequence[bool]]) -> List[Tuple[int, int, int, int]]:
    rows = len(mask)
    if rows == 0:
        return []
    cols = len(mask[0])
    used = [[False for _ in range(cols)] for _ in range(rows)]
    rects: List[Tuple[int, int, int, int]] = []

    for r0 in range(rows):
        for c0 in range(cols):
            if not mask[r0][c0] or used[r0][c0]:
                continue

            c1 = c0
            while c1 + 1 < cols and mask[r0][c1 + 1] and not used[r0][c1 + 1]:
                c1 += 1

            r1 = r0
            grow = True
            while grow and r1 + 1 < rows:
                rr = r1 + 1
                for cc in range(c0, c1 + 1):
                    if not mask[rr][cc] or used[rr][cc]:
                        grow = False
                        break
                if grow:
                    r1 = rr

            for rr in range(r0, r1 + 1):
                for cc in range(c0, c1 + 1):
                    used[rr][cc] = True
            rects.append((r0, r1, c0, c1))

    return rects


def _build_union_surface(boxes: Sequence[Box], merge_faces: bool = True) -> Tuple[List[Vec3], List[Face]]:
    xs = sorted({b[0] for b in boxes} | {b[1] for b in boxes})
    ys = sorted({b[2] for b in boxes} | {b[3] for b in boxes})
    zs = sorted({b[4] for b in boxes} | {b[5] for b in boxes})

    nx, ny, nz = len(xs) - 1, len(ys) - 1, len(zs) - 1
    occ = [[[False for _ in range(nz)] for _ in range(ny)] for _ in range(nx)]

    for i in range(nx):
        cx = 0.5 * (xs[i] + xs[i + 1])
        for j in range(ny):
            cy = 0.5 * (ys[j] + ys[j + 1])
            for k in range(nz):
                cz = 0.5 * (zs[k] + zs[k + 1])
                p = (cx, cy, cz)
                occ[i][j][k] = any(_point_in_box(p, b) for b in boxes)

    vertices: List[Vec3] = []
    faces: List[Face] = []
    vmap: Dict[Vec3, int] = {}

    if not merge_faces:
        for i in range(nx):
            x0, x1 = xs[i], xs[i + 1]
            for j in range(ny):
                y0, y1 = ys[j], ys[j + 1]
                for k in range(nz):
                    z0, z1 = zs[k], zs[k + 1]
                    if not occ[i][j][k]:
                        continue

                    # -X
                    if i == 0 or not occ[i - 1][j][k]:
                        _add_oriented_quad(
                            vertices,
                            faces,
                            vmap,
                            (x0, y0, z0),
                            (x0, y1, z0),
                            (x0, y1, z1),
                            (x0, y0, z1),
                            (-1.0, 0.0, 0.0),
                        )
                    # +X
                    if i == nx - 1 or not occ[i + 1][j][k]:
                        _add_oriented_quad(
                            vertices,
                            faces,
                            vmap,
                            (x1, y0, z0),
                            (x1, y0, z1),
                            (x1, y1, z1),
                            (x1, y1, z0),
                            (1.0, 0.0, 0.0),
                        )
                    # -Y
                    if j == 0 or not occ[i][j - 1][k]:
                        _add_oriented_quad(
                            vertices,
                            faces,
                            vmap,
                            (x0, y0, z0),
                            (x0, y0, z1),
                            (x1, y0, z1),
                            (x1, y0, z0),
                            (0.0, -1.0, 0.0),
                        )
                    # +Y
                    if j == ny - 1 or not occ[i][j + 1][k]:
                        _add_oriented_quad(
                            vertices,
                            faces,
                            vmap,
                            (x0, y1, z0),
                            (x1, y1, z0),
                            (x1, y1, z1),
                            (x0, y1, z1),
                            (0.0, 1.0, 0.0),
                        )
                    # -Z
                    if k == 0 or not occ[i][j][k - 1]:
                        _add_oriented_quad(
                            vertices,
                            faces,
                            vmap,
                            (x0, y0, z0),
                            (x1, y0, z0),
                            (x1, y1, z0),
                            (x0, y1, z0),
                            (0.0, 0.0, -1.0),
                        )
                    # +Z
                    if k == nz - 1 or not occ[i][j][k + 1]:
                        _add_oriented_quad(
                            vertices,
                            faces,
                            vmap,
                            (x0, y0, z1),
                            (x0, y1, z1),
                            (x1, y1, z1),
                            (x1, y0, z1),
                            (0.0, 0.0, 1.0),
                        )
        return vertices, faces

    xminus_masks = [[[False for _ in range(nz)] for _ in range(ny)] for _ in range(len(xs))]
    xplus_masks = [[[False for _ in range(nz)] for _ in range(ny)] for _ in range(len(xs))]
    yminus_masks = [[[False for _ in range(nz)] for _ in range(nx)] for _ in range(len(ys))]
    yplus_masks = [[[False for _ in range(nz)] for _ in range(nx)] for _ in range(len(ys))]
    zminus_masks = [[[False for _ in range(ny)] for _ in range(nx)] for _ in range(len(zs))]
    zplus_masks = [[[False for _ in range(ny)] for _ in range(nx)] for _ in range(len(zs))]

    for i in range(nx):
        x0, x1 = xs[i], xs[i + 1]
        for j in range(ny):
            y0, y1 = ys[j], ys[j + 1]
            for k in range(nz):
                z0, z1 = zs[k], zs[k + 1]
                if not occ[i][j][k]:
                    continue

                if i == 0 or not occ[i - 1][j][k]:
                    xminus_masks[i][j][k] = True
                if i == nx - 1 or not occ[i + 1][j][k]:
                    xplus_masks[i + 1][j][k] = True
                if j == 0 or not occ[i][j - 1][k]:
                    yminus_masks[j][i][k] = True
                if j == ny - 1 or not occ[i][j + 1][k]:
                    yplus_masks[j + 1][i][k] = True
                if k == 0 or not occ[i][j][k - 1]:
                    zminus_masks[k][i][j] = True
                if k == nz - 1 or not occ[i][j][k + 1]:
                    zplus_masks[k + 1][i][j] = True

    for xi, mask in enumerate(xminus_masks):
        x = xs[xi]
        for j0, j1, k0, k1 in _merge_rectangles(mask):
            _add_oriented_quad(
                vertices,
                faces,
                vmap,
                (x, ys[j0], zs[k0]),
                (x, ys[j1 + 1], zs[k0]),
                (x, ys[j1 + 1], zs[k1 + 1]),
                (x, ys[j0], zs[k1 + 1]),
                (-1.0, 0.0, 0.0),
            )

    for xi, mask in enumerate(xplus_masks):
        x = xs[xi]
        for j0, j1, k0, k1 in _merge_rectangles(mask):
            _add_oriented_quad(
                vertices,
                faces,
                vmap,
                (x, ys[j0], zs[k0]),
                (x, ys[j0], zs[k1 + 1]),
                (x, ys[j1 + 1], zs[k1 + 1]),
                (x, ys[j1 + 1], zs[k0]),
                (1.0, 0.0, 0.0),
            )

    for yj, mask in enumerate(yminus_masks):
        y = ys[yj]
        for i0, i1, k0, k1 in _merge_rectangles(mask):
            _add_oriented_quad(
                vertices,
                faces,
                vmap,
                (xs[i0], y, zs[k0]),
                (xs[i0], y, zs[k1 + 1]),
                (xs[i1 + 1], y, zs[k1 + 1]),
                (xs[i1 + 1], y, zs[k0]),
                (0.0, -1.0, 0.0),
            )

    for yj, mask in enumerate(yplus_masks):
        y = ys[yj]
        for i0, i1, k0, k1 in _merge_rectangles(mask):
            _add_oriented_quad(
                vertices,
                faces,
                vmap,
                (xs[i0], y, zs[k0]),
                (xs[i1 + 1], y, zs[k0]),
                (xs[i1 + 1], y, zs[k1 + 1]),
                (xs[i0], y, zs[k1 + 1]),
                (0.0, 1.0, 0.0),
            )

    for zk, mask in enumerate(zminus_masks):
        z = zs[zk]
        for i0, i1, j0, j1 in _merge_rectangles(mask):
            _add_oriented_quad(
                vertices,
                faces,
                vmap,
                (xs[i0], ys[j0], z),
                (xs[i1 + 1], ys[j0], z),
                (xs[i1 + 1], ys[j1 + 1], z),
                (xs[i0], ys[j1 + 1], z),
                (0.0, 0.0, -1.0),
            )

    for zk, mask in enumerate(zplus_masks):
        z = zs[zk]
        for i0, i1, j0, j1 in _merge_rectangles(mask):
            _add_oriented_quad(
                vertices,
                faces,
                vmap,
                (xs[i0], ys[j0], z),
                (xs[i0], ys[j1 + 1], z),
                (xs[i1 + 1], ys[j1 + 1], z),
                (xs[i1 + 1], ys[j0], z),
                (0.0, 0.0, 1.0),
            )

    return vertices, faces


def _build_gaafet(
    pad_length: float,
    pad_width: float,
    pad_height: float,
    channel_length: float,
    channel_width: float,
    channel_height: float,
    channel_count: int = 2,
    merge_faces: bool = True,
) -> Tuple[List[Vec3], List[Face]]:
    """Build two full-height source/drain pads connected by uniformly spaced channels.

    The channels are centered in Y and distributed uniformly along Z.
    For channel_count = 2, this gives equal gaps below, between, and above
    the channels.
    """
    if min(pad_length, pad_width, pad_height, channel_length, channel_width, channel_height) <= 0.0:
        raise ValueError("All dimensions must be positive.")
    if channel_count < 1:
        raise ValueError("channel_count must be >= 1.")
    if channel_width > pad_width:
        raise ValueError("channel_width must be <= pad_width.")
    if channel_count * channel_height >= pad_height:
        raise ValueError("channel_count * channel_height must be < pad_height so gaps can be formed.")

    x0 = 0.0
    x1 = pad_length
    x2 = pad_length + channel_length
    x3 = pad_length + channel_length + pad_length

    py0, py1 = -0.5 * pad_width, 0.5 * pad_width
    cy0, cy1 = -0.5 * channel_width, 0.5 * channel_width

    pz0, pz1 = 0.0, pad_height
    left_pad: Box = (x0, x1, py0, py1, pz0, pz1)
    right_pad: Box = (x2, x3, py0, py1, pz0, pz1)

    # Uniform Z distribution: bottom gap = middle gap = top gap.
    # Default: pad_height=30, channel_count=2, channel_height=6 -> gap=6 nm.
    z_gap = (pad_height - channel_count * channel_height) / (channel_count + 1)
    channels: List[Box] = []
    for n in range(channel_count):
        cz0 = z_gap + n * (channel_height + z_gap)
        cz1 = cz0 + channel_height
        channels.append((x1, x2, cy0, cy1, cz0, cz1))

    return _build_union_surface([left_pad, *channels, right_pad], merge_faces=merge_faces)


def _make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Generate double-channel GAAFET-like STL using two pads plus uniformly distributed channels.")
    p.add_argument("--output", default="model/gaafet.stl", help="Output STL path")
    p.add_argument("--name", default="gaafet_double_channel", help="STL header name")
    p.add_argument("--unit", choices=["nm", "m", "angstrom"], default="nm", help="Unit of input dimensions")
    p.add_argument("--scale", type=float, default=1.0, help="Uniform scale factor")

    p.add_argument("--pad-length", type=float, default=10.0, help="Pad length")
    p.add_argument("--pad-width", type=float, default=20.0, help="Pad width")
    p.add_argument("--pad-height", type=float, default=30.0, help="Overall height / pad height")

    p.add_argument("--channel-length", type=float, default=32.0, help="Channel length")
    p.add_argument("--channel-width", type=float, default=10.0, help="Channel width")
    p.add_argument("--channel-height", type=float, default=6.0, help="Channel height")
    p.add_argument("--channel-count", type=int, default=2, help="Number of uniformly distributed channels")
    p.add_argument("--no-merge-faces", action="store_true", help="Disable coplanar face merge")
    p.add_argument("--no-preview", action="store_true", help="Disable preview PNG generation")
    p.add_argument("--preview-output", default="model/gaafet_preview.png", help="Preview PNG path")
    return p


def main() -> int:
    args = _make_parser().parse_args()

    if args.unit == "nm":
        unit_to_nm = 1.0
    elif args.unit == "m":
        unit_to_nm = 1e9
    else:
        unit_to_nm = 0.1

    s = float(args.scale) * unit_to_nm
    vertices, faces = _build_gaafet(
        pad_length=args.pad_length * s,
        pad_width=args.pad_width * s,
        pad_height=args.pad_height * s,
        channel_length=args.channel_length * s,
        channel_width=args.channel_width * s,
        channel_height=args.channel_height * s,
        channel_count=args.channel_count,
        merge_faces=not args.no_merge_faces,
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    _write_binary_stl(out, vertices, faces, args.name)

    print(f"[ok] STL written: {out.resolve()}")
    z_gap = (args.pad_height - args.channel_count * args.channel_height) / (args.channel_count + 1)
    print(
        "[ok] Geometry: two pads + "
        f"{args.channel_count} channels; "
        f"overall W={args.pad_width:g} nm, H={args.pad_height:g} nm; "
        f"channel W={args.channel_width:g} nm, H={args.channel_height:g} nm"
    )
    print(f"[ok] Uniform vertical gap: {z_gap:g} nm")
    for n in range(args.channel_count):
        cz0 = z_gap + n * (args.channel_height + z_gap)
        cz1 = cz0 + args.channel_height
        print(f"[ok] Channel {n + 1} Z-range: {cz0:g} to {cz1:g} nm")
    print(f"[ok] Coplanar merge: {'off' if args.no_merge_faces else 'on'}")
    print(f"[ok] Total vertices: {len(vertices)}")
    print(f"[ok] Total faces: {len(faces)}")
    if not args.no_preview:
        _save_preview_png(vertices, faces, Path(args.preview_output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
