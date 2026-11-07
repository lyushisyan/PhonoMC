#!/usr/bin/env python3
"""Generate matched FinFET/GAAFET surface meshes with a common substrate (nm).

Only geometry is generated; no phonon transport calculation is performed.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import struct

from generate_gaafet_stl import _build_union_surface, _write_binary_stl


def validate_stl(path, expected_volume):
    """Validate the serialized STL, including edge winding and connectivity."""
    data = path.read_bytes()
    count = struct.unpack_from('<I', data, 80)[0]
    assert len(data) == 84 + 50 * count
    edges = defaultdict(list)
    vertices = set()
    volume = 0.0
    for i in range(count):
        a, b, c = [struct.unpack_from('<3f', data, 84 + i * 50 + 12 + j * 12)
                   for j in range(3)]
        vertices.update((a, b, c))
        u, v = [b[k] - a[k] for k in range(3)], [c[k] - a[k] for k in range(3)]
        cross = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
        assert sum(x*x for x in cross) > 0, 'Degenerate triangle'
        volume += (a[0]*(b[1]*c[2]-b[2]*c[1]) +
                   a[1]*(b[2]*c[0]-b[0]*c[2]) +
                   a[2]*(b[0]*c[1]-b[1]*c[0])) / 6
        for p, q in ((a,b), (b,c), (c,a)):
            edges[tuple(sorted((p,q)))].append((i, p, q))
    adjacency = defaultdict(set)
    for uses in edges.values():
        assert len(uses) == 2, 'Non-manifold or open edge'
        first, second = uses
        assert first[1:] == second[:0:-1], 'Inconsistent winding'
        adjacency[first[0]].add(second[0])
        adjacency[second[0]].add(first[0])
    seen, pending = set(), [0]
    while pending:
        i = pending.pop()
        if i not in seen:
            seen.add(i)
            pending.extend(adjacency[i] - seen)
    assert len(seen) == count, 'Disconnected surface'
    assert volume > 0 and math.isclose(volume, expected_volume, rel_tol=2e-6), (volume, expected_volume)
    return dict(triangles=count, vertices=len(vertices), closed_manifold=True,
                consistent_outward_winding=True, connected=True,
                volume_nm3=volume, expected_volume_nm3=expected_volume,
                bounds_nm=[[min(v[k] for v in vertices), max(v[k] for v in vertices)] for k in range(3)])


def box_volume(box):
    return math.prod(box[k+1]-box[k] for k in (0,2,4))


def validate_periodic_faces(path):
    """Compare opposite exterior face areas and projected boundary edges in X/Y.

    This is a geometry check, not a solver periodic-facet pairing test.
    Removing internal triangulation edges allows opposite diagonal choices.
    """
    data = path.read_bytes()
    count = struct.unpack_from('<I', data, 80)[0]
    triangles = [[struct.unpack_from('<3f', data, 84+i*50+12+j*12)
                  for j in range(3)] for i in range(count)]
    result = {}
    for axis, name in enumerate(('x','y')):
        low = min(p[axis] for tri in triangles for p in tri)
        high = max(p[axis] for tri in triangles for p in tri)
        sides = []
        for plane, sign in ((low,-1),(high,1)):
            edges = Counter()
            area = 0.0
            ntri = 0
            for tri in triangles:
                if not all(p[axis] == plane for p in tri):
                    continue
                a,b,c = tri
                u,v = [b[k]-a[k] for k in range(3)], [c[k]-a[k] for k in range(3)]
                normal = (u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])
                assert sign*normal[axis] > 0, 'Incorrect periodic face normal'
                area += sign*normal[axis]/2
                ntri += 1
                projected = [tuple(p[k] for k in range(3) if k != axis) for p in tri]
                for j in range(3):
                    edges[tuple(sorted((projected[j],projected[(j+1)%3])))]+=1
            assert ntri and all(n in (1,2) for n in edges.values())
            sides.append((area,{edge for edge,n in edges.items() if n == 1},ntri))
        assert math.isclose(sides[0][0],sides[1][0],rel_tol=1e-7), 'Periodic face areas differ'
        assert sides[0][1] == sides[1][1], 'Periodic face outlines differ'
        shift = [0.,0.,0.]
        shift[axis] = high-low
        result[name] = dict(translation_nm=shift,area_each_nm2=sides[0][0],
                            triangle_counts=[s[2] for s in sides],
                            projected_boundaries_identical=True,opposite_normals=True,
                            solver_pairing_verified=False)
    return result


def write_mesh(path, boxes):
    # Keep the shared coordinate lattice: rectangle merging can create T-junctions.
    vertices, faces = _build_union_surface(boxes, merge_faces=False)
    _write_binary_stl(path, vertices, faces, path.stem)
    # Components have disjoint interiors; touching faces disappear in the union.
    return validate_stl(path, sum(map(box_volume, boxes)))


def preview(out, devices, params, heat_length=None, heat_profile='uniform'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb
    from matplotlib.cm import ScalarMappable
    import numpy as np
    from matplotlib.patches import Patch, Rectangle
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    colors = dict(substrate='#b5bdc9', source='#e5aa57', drain='#e5aa57', channel='#368ac1',heat='#d74737')
    fig = plt.figure(figsize=(13,9), facecolor='white')
    fig.suptitle('Proposed drain-side heat sources' if heat_length else
                 'FinFET and GAAFET on a common substrate', fontsize=20, y=.97)
    sx, sy, sz = params['substrate_size_nm']
    total_height = sz + params['pad_height']
    gaussian = heat_length is not None and heat_profile == 'gaussian'
    cmap = plt.get_cmap('YlOrRd')
    if gaussian:
        fig.suptitle('Drain-side Gaussian heat sources',fontsize=20,y=.97)
    for col, (name, parts) in enumerate(devices.items()):
        if heat_length:
            channels, heated = [], []
            for x0,x1,y0,y1,z0,z1 in parts['channel']:
                if x1-heat_length > x0:
                    channels.append((x0,x1-heat_length,y0,y1,z0,z1))
                heated.append((x1-heat_length,x1,y0,y1,z0,z1))
            parts = dict(parts,channel=channels,heat=heated)
        ax = fig.add_subplot(2,2,col+1, projection='3d')
        polygons, facecolors = [], []
        boxes = [box for group in parts.values() for box in group]
        if gaussian:
            # Add fine X coordinates for color sampling without changing the STL.
            for x0,x1,y0,y1,z0,z1 in parts['heat']:
                xs = np.linspace(x0,x1,49)
                boxes.extend((a,b,y0,y1,z0,z1) for a,b in zip(xs[:-1],xs[1:]))
            center_x = sum(parts['heat'][0][:2])/2
            sigma_x = heat_length/4
        def surface_color(part, x):
            if part == 'heat' and gaussian:
                return cmap(math.exp(-.5*((x-center_x)/sigma_x)**2))[:3]
            return to_rgb(colors[part])
        verts, faces = _build_union_surface(boxes, merge_faces=False)
        azimuth_degrees = -85 if heat_length else -62
        elevation, azimuth = math.radians(23), math.radians(azimuth_degrees)
        camera = (math.cos(elevation)*math.cos(azimuth),
                  math.cos(elevation)*math.sin(azimuth),math.sin(elevation))
        # Render the actual union exterior, with no internal contact faces.
        # Its small lattice quads avoid incorrect whole-box depth ordering in
        # mplot3d. Cull rear faces for this orthographic view before sorting.
        for first, second in zip(faces[::2],faces[1::2]):
            edges = [(tri[j],tri[(j+1)%3]) for tri in (first,second) for j in range(3)]
            perimeter = dict((a,b) for a,b in edges if (b,a) not in edges)
            indices = [first[0]]
            for _ in range(3):
                indices.append(perimeter[indices[-1]])
            poly = [verts[i] for i in indices]
            a,b,c = poly[:3]
            u,v = [b[k]-a[k] for k in range(3)], [c[k]-a[k] for k in range(3)]
            normal = (u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])
            if sum(n*d for n,d in zip(normal,camera)) <= 0:
                continue
            center = [sum(p[k] for p in poly)/4 for k in range(3)]
            part = next(part for part,group in parts.items() for box in group
                        if all(box[2*k]-1e-8 <= center[k] <= box[2*k+1]+1e-8 for k in range(3)))
            axis = max(range(3),key=lambda k:abs(normal[k]))
            shade = (.72,.84,1.)[axis]
            polygons.append(poly)
            facecolors.append(tuple(c*shade for c in surface_color(part,center[0])))
        ax.add_collection3d(Poly3DCollection(polygons,facecolors=facecolors,
                                             edgecolors=facecolors,linewidths=.25,
                                             antialiased=False))
        ax.set_proj_type('ortho')
        ax.set(xlim=(0,sx), ylim=(0,sy), zlim=(0,total_height),
               xlabel='X / transport (nm)', ylabel='Y (nm)', zlabel='Z (nm)', title=name)
        ax.set_box_aspect((sx,sy,total_height))
        ax.set_yticks([0,sy/2,sy])
        ax.view_init(elev=23, azim=azimuth_degrees)
        section = fig.add_subplot(2,2,col+3)
        if heat_length:
            for part, group in parts.items():
                for x0,x1,_,_,z0,z1 in group:
                    section.add_patch(Rectangle((x0,z0),x1-x0,z1-z0,
                                               facecolor=colors[part],edgecolor='#374151',lw=.7))
                    if part == 'heat':
                        if gaussian:
                            x = np.linspace(x0,x1,300)
                            intensity = np.exp(-.5*((x-center_x)/sigma_x)**2)
                            section.imshow(np.tile(intensity,(2,1)),extent=(x0,x1,z0,z1),
                                           origin='lower',cmap=cmap,vmin=0,vmax=1,aspect='auto',zorder=2)
                        label = 'P' if name == 'FinFET' else f'P/{len(parts["heat"])}'
                        section.text((x0+x1)/2,(z0+z1)/2,label,ha='center',va='center',color='white',fontsize=9,zorder=3)
                    elif part in ('source','drain'):
                        section.text((x0+x1)/2,(z0+z1)/2,part.title(),rotation=90,ha='center',va='center',fontsize=10)
            start,end = parts['heat'][0][:2]
            section.annotate('',xy=(start,total_height+3),xytext=(end,total_height+3),
                             arrowprops=dict(arrowstyle='<->',color=colors['heat']))
            section.text((start+end)/2,total_height+4,f'{heat_length:g} nm',ha='center',color=colors['heat'],fontsize=10)
            section.set(xlim=(0,sx),ylim=(0,total_height+9),xlabel='X (nm)',ylabel='Z (nm)',
                        title=f'Center slice at Y = {sy/2:g} nm; heat: X = {start:g}–{end:g} nm')
            section.set_aspect('equal')
            section.set_xticks(sorted(set([0,parts['source'][0][1],start,end,sx])))
            section.text(sx/2,sz/2,'Substrate',ha='center',color='#374151')
            continue
        for part in ('substrate','channel'):
            for _,_,y0,y1,z0,z1 in parts[part]:
                section.add_patch(Rectangle((y0,z0),y1-y0,z1-z0,facecolor=colors[part],edgecolor='#374151'))
                if part == 'channel':
                    section.text((y0+y1)/2,(z0+z1)/2,f'{y1-y0:g} x {z1-z0:g}',ha='center',va='center',color='white',fontsize=8)
        section.axhline(sz,color='#64748b',ls=':',lw=1)
        section.set(xlim=(0,sy),ylim=(0,total_height),xlabel='Y (nm)',ylabel='Z (nm)',
                    title='Channel midpoint cross-section (dimensions in nm)')
        section.set_aspect('equal')
        section.text(sy/2,sz/2,'Substrate',ha='center',color='#374151')
    legend = [('substrate','Substrate'),('source','Source / drain'),('channel','Fin / channels')]
    if heat_length:
        legend.append(('heat','Volumetric heat source'))
    fig.legend(handles=[Patch(color=colors[k],label=v) for k,v in legend],
               loc='lower center',ncol=len(legend),frameon=False,bbox_to_anchor=(.5,.018))
    note = (f'Same total power P per FET; GAAFET: P/{len(devices["GAAFET"]["channel"])} per layer. Proposed placement only; no simulation run.'
            if heat_length else 'Geometry only; gate and oxide are not included.')
    fig.text(.5,.008,note,ha='center',fontsize=10,color='#64748b')
    fig.tight_layout(rect=(0,.06,1,.94), h_pad=4.0)
    if gaussian:
        cax = fig.add_axes([.91,.22,.012,.18])
        cb = fig.colorbar(ScalarMappable(cmap=cmap,norm=plt.Normalize(0,1)),cax=cax)
        cb.set_label('Relative heating q / q_peak',fontsize=9)
        cb.ax.tick_params(labelsize=8)
    stem = 'heat_sources_gaussian' if gaussian else ('heat_sources' if heat_length else 'comparison')
    fig.savefig(out/(stem+'.png'),dpi=180)
    fig.savefig(out/(stem+'.svg'))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=Path('model/matched_fets'))
    for name, default in dict(channel_length=32.,pad_length=10.,pad_width=24.,pad_height=30.,
                              fin_width=6.,fin_height=20.,sheet_width=10.,sheet_height=6.,
                              substrate_thickness=20.,margin_x=0.,margin_y=0.).items():
        parser.add_argument('--'+name.replace('_','-'),type=float,default=default)
    parser.add_argument('--sheet-count',type=int,default=2)
    parser.add_argument('--no-preview',action='store_true')
    args = parser.parse_args()
    numeric = [v for k,v in vars(args).items() if isinstance(v,(float,int))
               and not isinstance(v,bool) and k not in ('margin_x','margin_y')]
    if not all(math.isfinite(v) and v > 0 for v in numeric):
        parser.error('All dimensions and sheet count must be finite and positive')
    if not all(math.isfinite(v) and v >= 0 for v in (args.margin_x,args.margin_y)):
        parser.error('Substrate margins must be finite and non-negative')
    if args.fin_width > args.pad_width or args.sheet_width > args.pad_width:
        parser.error('Fin/sheet width must fit inside the source/drain')
    if args.fin_height > args.pad_height or args.sheet_count*args.sheet_height >= args.pad_height:
        parser.error('Fin/sheets must fit within the pad height with positive sheet gaps')
    out = args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    x0 = args.margin_x
    x1 = x0 + args.pad_length
    x2 = x1 + args.channel_length
    x3 = x2 + args.pad_length
    sx,sy,sz = x3 + args.margin_x, args.pad_width + 2*args.margin_y, args.substrate_thickness
    yc = sy/2
    common = dict(substrate=[(0,sx,0,sy,0,sz)],
                  source=[(x0,x1,args.margin_y,sy-args.margin_y,sz,sz+args.pad_height)],
                  drain=[(x2,x3,args.margin_y,sy-args.margin_y,sz,sz+args.pad_height)])
    fin = (x1,x2,yc-args.fin_width/2,yc+args.fin_width/2,sz,sz+args.fin_height)
    gap = (args.pad_height-args.sheet_count*args.sheet_height)/(args.sheet_count+1)
    sheets = [(x1,x2,yc-args.sheet_width/2,yc+args.sheet_width/2,
               sz+gap+i*(args.sheet_height+gap),sz+gap+i*(args.sheet_height+gap)+args.sheet_height)
              for i in range(args.sheet_count)]
    devices = {'FinFET':dict(common,channel=[fin]), 'GAAFET':dict(common,channel=sheets)}
    params = {k:v for k,v in vars(args).items() if k not in ('output_dir','no_preview')}
    params['substrate_size_nm'] = [sx,sy,sz]
    report = dict(units='nm',box_order=['xmin','xmax','ymin','ymax','zmin','zmax'],
                  parameters=params,sheet_gap_nm=gap,devices={},transport_simulated=False,
                  material_assignment=None,gate_and_oxide_included=False)
    for name, parts in devices.items():
        boxes = [b for group in parts.values() for b in group]
        slug = name.lower()
        checks = {slug+'_with_substrate.stl':write_mesh(out/(slug+'_with_substrate.stl'),boxes)}
        comp = out/(slug+'_parts')
        comp.mkdir(exist_ok=True)
        for part, group in parts.items():
            for i, box in enumerate(group):
                path = comp/f'{part}_{i+1}.stl'
                checks[str(path.relative_to(out))] = write_mesh(path,[box])
        report['devices'][name] = dict(boxes=parts,validation=checks,
                periodic_geometry=validate_periodic_faces(out/(slug+'_with_substrate.stl')),
                total_channel_area_nm2=sum((b[3]-b[2])*(b[5]-b[4]) for b in parts['channel']),
                channel_volume_nm3=sum(map(box_volume,parts['channel'])))
    a,b = [d['validation'][name.lower()+'_with_substrate.stl']['bounds_nm'] for name,d in report['devices'].items()]
    assert a == b, 'Device envelopes differ'
    (out/'geometry.json').write_text(json.dumps(report,indent=2)+'\n')
    if not args.no_preview:
        preview(out,devices,params)
        preview(out,devices,params,heat_length=min(6.,args.channel_length))
    print(json.dumps({name: {'area_nm2':d['total_channel_area_nm2'],
                     'mesh':d['validation'][name.lower()+'_with_substrate.stl']}
                     for name,d in report['devices'].items()},indent=2))


if __name__ == '__main__':
    main()
