#!/usr/bin/env python3
"""Draw FinFET variants with a radial 3-sigma source cutoff; no simulation.

Geometry and source center follow the old case. The radial cutoff is a new
plotting specification, not the old solver's per-axis 2-sigma cutoff.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "model/finfet_reference"
META = ROOT / "calculation/finfet_b40_sigma2_rta_mmm_20260915/DEPLOYMENT.json"
INK = "#222222"
BLUE = ["#e5eaf0", "#cdd7e2", "#f2f4f7"]
GREEN = ["#e1eae3", "#c7d8ca", "#f0f4f0"]
BAS = ["#f0e8d9", "#dfd0b6", "#f7f2e9"]
HOT = "#be342b"
SIGMA_NM = 2.0
CUTOFF_SIGMA = 3.0


def source_density(y, z):
    radius2 = ((np.asarray(y)-25)/SIGMA_NM)**2 + ((np.asarray(z)-60)/SIGMA_NM)**2
    return np.where(radius2 <= CUTOFF_SIGMA**2, np.exp(-.5*radius2), 0.)


def project(p):
    p = np.asarray(p)
    x, y, z = np.moveaxis(p, -1, 0)
    return np.stack((x + .64*y, z + .36*y - .20*x), axis=-1)


def face(ax, points, color):
    ax.add_patch(Polygon(project(points), closed=True, facecolor=color,
                         edgecolor=INK, linewidth=.75, joinstyle="miter"))


def line(ax, points, **kwargs):
    p = project(points)
    ax.plot(p[:, 0], p[:, 1], **kwargs)


def dimension(ax, a, b, label, offset=(0, 0), text_offset=(0, 0)):
    a, b = project([a, b])
    d = np.array(offset)
    for p in (a, b):
        ax.plot([p[0], p[0]+d[0]], [p[1], p[1]+d[1]], color="#777777", lw=.5)
    ax.annotate("", xy=b+d, xytext=a+d,
                arrowprops=dict(arrowstyle="<->", lw=.65, color=INK, shrinkA=0, shrinkB=0, mutation_scale=7))
    mid = (a+b)/2+d+text_offset
    ax.text(*mid, label, ha="center", va="center", fontsize=8.5, color=INK,
            bbox=dict(fc="white", ec="none", pad=.6))


def heat_surface(ax, top):
    # Uniform x. Radial cutoff in y-z gives a half-cylinder inside the fin.
    cutoff = CUTOFF_SIGMA * SIGMA_NM
    if top:
        Y, X = np.meshgrid(np.linspace(25-cutoff,25+cutoff,181),np.linspace(7,15,33))
        Z = np.full_like(Y,60)
    else:
        Y, Z = np.meshgrid(np.linspace(25-cutoff,25+cutoff,241),np.linspace(60-cutoff,60,121))
        X = np.full_like(Y,15)
    uv = project(np.stack([X, Y, Z], axis=-1))
    q = np.ma.masked_where(source_density(Y,Z)==0,source_density(Y,Z))
    # Discrete Gaussian contours keep the source legible in vector print output.
    ax.contourf(uv[...,0], uv[...,1], q, levels=[0,.15,.4,.7,1.01],
                colors=["#f6e7e5", "#efc4be", "#dd8b81", HOT], zorder=2)
    if not top:
        theta = np.linspace(0,np.pi,181)
        arc = np.column_stack([np.full_like(theta,15),25+cutoff*np.cos(theta),60-cutoff*np.sin(theta)])
        line(ax,arc,color="#92544f",lw=.55,ls=(0,(2,2)),zorder=3)


def structure(ax, substrate, dimensions=False):
    composite = substrate != "Si"
    lower = GREEN if substrate == "SiC" else BAS
    layers = [(0, 30, lower), (30, 40, BLUE)] if composite else [(0, 40, BLUE)]
    for low, high, colors in layers:
        face(ax, [(0,0,low),(22,0,low),(22,0,high),(0,0,high)], colors[0])
        face(ax, [(22,0,low),(22,36,low),(22,36,high),(22,0,high)], colors[1])
    for lo, hi in [(0,7),(15,22)]:
        face(ax, [(lo,0,40),(hi,0,40),(hi,36,40),(lo,36,40)], BLUE[2])
    face(ax, [(7,0,40),(15,0,40),(15,0,60),(7,0,60)], BLUE[0])
    face(ax, [(15,0,40),(15,36,40),(15,36,60),(15,0,60)], BLUE[1])
    face(ax, [(7,0,60),(15,0,60),(15,36,60),(7,36,60)], BLUE[2])
    heat_surface(ax, False)
    heat_surface(ax, True)
    line(ax, [(7,0,60),(15,0,60),(15,36,60),(7,36,60),(7,0,60)], color=INK, lw=.75, zorder=3)
    if composite:
        line(ax, [(0,0,30),(22,0,30),(22,36,30)], color=INK, lw=.85, ls=(0,(3,2)))
        ax.text(*project([11,0,15]), substrate+"\n30 nm", ha="center", va="center", fontsize=10, color=INK)
        ax.text(*project([11,0,35]), "Si\n10 nm", ha="center", va="center", fontsize=8.5, color=INK, linespacing=1.1)
    else:
        ax.text(*project([11,0,20]), "Si\n40 nm", ha="center", va="center", fontsize=10, color=INK)
    ax.text(*project([11,0,49]), "Si", ha="center", va="center", fontsize=9.5, color=INK)
    line(ax, [(0,0,0),(22,0,0),(22,36,0)], color="#387a96", lw=2)
    ax.text(20,-17,"300 K",fontsize=10,color=INK,ha="center")
    ax.annotate("Heat source ($3\\sigma$ cutoff)", xy=project([15,25,59]), xytext=(27,77),
                fontsize=9, color=INK, ha="center",
                arrowprops=dict(arrowstyle="->", color=INK, lw=.7, mutation_scale=8))
    if dimensions:
        dimension(ax, [0,0,0], [22,0,0], "22 nm", (0,-7))
        dimension(ax, [22,0,0], [22,36,0], "36 nm", (7,-5), (2,0))
        dimension(ax, [0,0,0], [0,0,40], "40 nm", (-8,0))
        dimension(ax, [7,0,40], [7,0,60], "20 nm", (-15,0))
        dimension(ax, [7,0,60], [15,0,60], "8 nm", (0,7))
    # Coordinate triad uses the same projection as the structure.
    if dimensions:
        origin = np.array([51., 19.])
        for vec, label in [([6,0,0],"$x$"),([0,8,0],"$y$"),([0,0,7],"$z$")]:
            end = origin + project(vec)
            ax.annotate("", xy=end, xytext=origin, arrowprops=dict(arrowstyle="->", lw=.65, color=INK, mutation_scale=7))
            ax.text(*(end+[1,.5]), label, fontsize=9, color=INK)
    ax.set(xlim=(-16,65), ylim=(-21,83), aspect="equal")
    ax.axis("off")


def main():
    meta = json.loads(META.read_text())
    assert meta["bounds_nm"] == [22,36,60]
    assert meta["fin_width_nm"] == 8 and meta["substrate_thickness_nm"] == 40
    assert meta["source_position_nm"] == [11,25,60] and meta["source_sigma_nm"] == [0,2,2]
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family":"Arial", "font.size":10,
                         "mathtext.fontset":"dejavusans", "svg.fonttype":"none", "pdf.fonttype":42})
    fig = plt.figure(figsize=(11.7,5.1), facecolor="white")
    for i, (material, title) in enumerate([("Si","(a) Si"),("SiC","(b) Si/SiC"),("BAs","(c) Si/BAs")]):
        ax = fig.add_axes([.015+i*.328,.115,.32,.80])
        structure(ax, material, dimensions=(i==0))
        ax.text(.06,1.02,title,transform=ax.transAxes,fontsize=11,color=INK)
    fig.text(.5,.052,r"$\sigma_y=\sigma_z=\sigma=2$ nm; $q=0$ for $r>3\sigma=6$ nm, where $r^2=(y-y_0)^2+(z-z_0)^2$.",
             ha="center",fontsize=9,color=INK)
    fig.text(.5,.016,r"Uniform across the fin width.   Periodic boundaries in $x$ and $y$.",
             ha="center",fontsize=9,color=INK)
    for ext in ("png","svg","pdf"):
        fig.savefig(OUT/f"finfet_structure.{ext}",dpi=240,facecolor="white")
    plt.close(fig)
    print(OUT / "finfet_structure.png")


if __name__ == "__main__":
    main()
