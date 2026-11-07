#!/usr/bin/env python3
"""Build a three-sheet GAAFET comparison from the existing matched FinFET.

Geometry, input preparation and static audits only; never runs transport.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tomllib

from generate_matched_fets import box_volume, validate_stl, write_mesh
from generate_gaafet_stl import _build_union_surface

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def boundary_kind(center, normal, size, substrate_height):
    x, y, z = center
    sx, sy, sz = size
    if abs(z) < 1e-7 or abs(z-sz) < 1e-7:
        return 'thermal'
    if z < substrate_height and (abs(normal[0])+abs(normal[1]) > .5):
        if min(abs(x),abs(x-sx),abs(y),abs(y-sy)) < 1e-7:
            return 'periodic'
    return 'rough'


def plot_model(out, devices, size, power):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb
    import numpy as np
    from matplotlib.patches import Patch, Rectangle
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    colors = dict(thermal='#179d9b', periodic='#8979b5', rough='#b7c6d4', heat='#c53d37')
    fig = plt.figure(figsize=(12,9), facecolor='white')
    fig.suptitle('Matched FinFET and three-sheet GAAFET', fontsize=19)
    for col,(name,parts) in enumerate(devices.items()):
        boxes = [b for group in parts.values() for b in group]
        # Partition surfaces at the source edge for drawing only.
        heats = [(36,42,b[2],b[3],b[4],b[5]) for b in parts['channel']]
        verts, faces = _build_union_surface(boxes+heats, merge_faces=False)
        ax = fig.add_subplot(2,2,col+1,projection='3d')
        polygons, fc = [], []
        el,az = math.radians(23),math.radians(-85)
        camera = np.array([math.cos(el)*math.cos(az),math.cos(el)*math.sin(az),math.sin(el)])
        for face in faces:
            poly=np.array([verts[i] for i in face]); normal=np.cross(poly[1]-poly[0],poly[2]-poly[0])
            if normal@camera <= 0: continue
            normal=normal/np.linalg.norm(normal); center=poly.mean(axis=0)
            kind=boundary_kind(center,normal,size,20)
            if any(all(b[2*k]-1e-7<=center[k]<=b[2*k+1]+1e-7 for k in range(3)) for b in heats):
                kind='heat'
            shade=(.75,.86,1.)[int(np.argmax(np.abs(normal)))]
            polygons.append(poly);fc.append(tuple(v*shade for v in to_rgb(colors[kind])))
        ax.add_collection3d(Poly3DCollection(polygons,facecolors=fc,edgecolors=fc,linewidths=.3))
        ax.set(xlim=(0,size[0]),ylim=(0,size[1]),zlim=(0,size[2]),xlabel='X (nm)',ylabel='Y (nm)',zlabel='Z (nm)',title=name)
        ax.set_box_aspect(size);ax.set_proj_type('ortho');ax.view_init(elev=23,azim=-85)
        sl=fig.add_subplot(2,2,col+3)
        for part,group in parts.items():
            for x0,x1,y0,y1,z0,z1 in group:
                sl.add_patch(Rectangle((x0,z0),x1-x0,z1-z0,facecolor=colors['rough'],edgecolor='white',lw=.5))
                if part=='channel':
                    xs=np.linspace(36,42,150)
                    sl.imshow(np.tile(np.exp(-.5*((xs-39)/1.5)**2),(2,1)),extent=(36,42,z0,z1),origin='lower',cmap='Reds',vmin=0,vmax=1,aspect='auto',zorder=2)
                    sl.text(25,(z0+z1)/2,f'{z1-z0:g} nm',ha='center',va='center',fontsize=9)
        for x0,x1 in [(0,10),(42,52)]:
            sl.plot([x0,x1],[50,50],color=colors['thermal'],lw=4)
            sl.text((x0+x1)/2,52,'300 K',ha='center',fontsize=10,color='#087c79')
        sl.plot([0,52],[0,0],color=colors['thermal'],lw=4)
        sl.text(26,3,'Substrate; bottom bath: 300 K',ha='center',fontsize=9)
        for x in [0,52]:sl.plot([x,x],[0,20],color=colors['periodic'],lw=4)
        sl.annotate('',xy=(50,13),xytext=(2,13),arrowprops=dict(arrowstyle='<->',color=colors['periodic'],lw=1.5))
        sl.text(26,15,'Periodic X; periodic Y out of plane',ha='center',fontsize=8,color='#665197')
        sl.text(5,35,'Source',rotation=90,ha='center',va='center',fontsize=10)
        sl.text(47,35,'Drain',rotation=90,ha='center',va='center',fontsize=10)
        sl.set(xlim=(-2,54),ylim=(-2,57),xlabel='X (nm)',ylabel='Z (nm)',title=f'Y = 12 nm; {power*1e6:g} µW total; {power*1e6/len(parts["channel"]):.3g} µW / channel')
        sl.set_aspect('equal')
    fig.legend(handles=[Patch(color=colors[k],label=v) for k,v in [('thermal','300 K thermal bath'),('periodic','Substrate sides: periodic X / Y'),('rough','Other surfaces: R, η = 1 nm'),('heat','Gaussian volume heating')]],loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.5,.025))
    fig.text(.5,.012,'Same Si material, contact areas, channel area and total power. Geometry and boundary setup; no transport results.',ha='center',fontsize=9,color='#526173')
    fig.tight_layout(rect=(0,.10,1,.95),h_pad=3)
    for suffix in ['png','svg']:fig.savefig(out/f'comparison.{suffix}',dpi=180)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--power-w',type=float,default=1e-6)
    p.add_argument('--output-dir',type=Path,default=ROOT/'model/matched_fets/three_layer_comparison')
    args=p.parse_args()
    if not math.isfinite(args.power_w) or args.power_w<=0:p.error('Power must be positive and finite')
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    source=ROOT/'model/matched_fets'
    old=json.loads((source/'geometry.json').read_text())
    fin=old['devices']['FinFET']['boxes']
    # Explicitly pin this comparison to the user-selected existing FinFET.
    assert fin['substrate']==[[0,52.,0,24.,0,20.]]
    assert fin['channel']==[[10.,42.,9.,15.,20.,40.]]
    assert fin['source']==[[0.,10.,0.,24.,20.,50.]]
    assert fin['drain']==[[42.,52.,0.,24.,20.,50.]]
    count,width=3,10.
    area=sum((b[3]-b[2])*(b[5]-b[4]) for b in fin['channel'])
    thickness=area/(count*width);gap=(30-count*thickness)/(count+1)
    gaa={k:v for k,v in fin.items() if k!='channel'}
    gaa['channel']=[[10,42,7,17,20+gap+i*(thickness+gap),20+gap+i*(thickness+gap)+thickness] for i in range(count)]
    devices={'FinFET':fin,'GAAFET':gaa};size=[52,24,50]
    # Thin selector slabs, with final overrides for thermal contacts. Keep
    # triangles unmerged: pad/substrate side walls share a coplanar plane.
    e=1e-6
    selectors=[[-e,-e,-e,1+e,1+e,1+e],[-e,-e,-e,e,1+e,.4-e],
               [1-e,-e,-e,1+e,1+e,.4-e],[-e,-e,-e,1+e,e,.4-e],
               [-e,1-e,-e,1+e,1+e,.4-e],[-e,-e,-e,1+e,1+e,e],
               [-e,-e,1-e,10/52+e,1+e,1+e],[42/52-e,-e,1-e,1+e,1+e,1+e]]
    conditions=['R']+['P']*4+['T']*3;values=[1.,0.,0.,0.,0.,300.,300.,300.]
    periodic_pairs=selectors[1:5]
    material=ROOT/'material/Si_ph3py_m252525_T50-800'
    import h5py
    import numpy as np
    with h5py.File(material/'kappa-fbz.hdf5') as f:
        temps=f['temperature'][:]
        g,n,u=[f[k][:] for k in ['gamma','gamma_N','gamma_U']]
        assert g.shape==n.shape==u.shape and all(np.isfinite(v).all() and (v>=0).all() for v in [g,n,u])
        assert np.allclose(g,n+u,rtol=1e-7,atol=1e-14)
        material_info=dict(directory=os.path.relpath(material,out),temperature_K=temps.tolist(),rate_shape=list(g.shape),NU_sum_verified=True,sha256={x:sha(material/x) for x in ['POSCAR','kappa-fbz.hdf5']})
    report=dict(units='nm',box_order=['xmin','xmax','ymin','ymax','zmin','zmax'],
                source_geometry_sha256=sha(source/'geometry.json'),sheet_count=count,sheet_gap_nm=gap,
                bounds_size_nm=size,grid_xyz=[52,24,100],grid_spacing_nm=[1,1,.5],
                material=material_info,total_power_per_device_W=args.power_w,
                boundary=dict(selectors=selectors,conditions=conditions,values=values,
                              periodic_pair=periodic_pairs,period_nm=[52,24],periodic_height_nm=[0,20],
                              merge_coplanar_facets=False,conservative_boundaries=True,
                              roughness_is_not_fully_diffuse=True),
                devices={},transport_simulated=False,
                solver_audit_path='validation/solver_audit.json',
                solver_audit_must_be_rerun_after_generation=True)
    for name,parts in devices.items():
        slug=name.lower();boxes=[b for group in parts.values() for b in group]
        mesh=out/f'{slug}.stl'
        if name=='FinFET':
            shutil.copyfile(source/'finfet_with_substrate.stl',mesh)
            check=validate_stl(mesh,sum(map(box_volume,boxes)))
            assert sha(mesh)==sha(source/'finfet_with_substrate.stl')
        else:check=write_mesh(mesh,boxes)
        channels=parts['channel'];lo=[36,min(b[2] for b in channels),min(b[4] for b in channels)]
        hi=[42,max(b[3] for b in channels),max(b[5] for b in channels)]
        center=[39,12,(lo[2]+hi[2])/2]
        source_text=f'''[heat_source]
enabled = true
profile = "gaussian"
min = {json.dumps(lo)}
max = {json.dumps(hi)}
center = {json.dumps(center)}
sigma = [1.5, 0.0, 0.0]
total_power = {args.power_w:.15g}
spectrum = "thermal"
time_profile = "constant"
start_time = 0
'''
        for model in ['rta','callaway']:
            collision=f'[scattering]\nmodel = "{model}"\n'
            if model=='callaway':collision+='callaway_formulation = "nonlinear"\n'
            text=f'''# Prepared baseline, not a converged production result.
# Same physical boundary and total device power for both devices/models.
# Keep facets unmerged so substrate sides and upper pad sides stay separate.
[geometry]
model = "{slug}.stl"
merge_coplanar_facets = false

[simulation]
particle_count = 10000000
time_step = 0.05
iterations = 40000
random_seed = 12345
grid_xyz = [52, 24, 100]
initial_temperature = 300
background_temperature = 300
lifetime_temperature = "local"
conservative_boundaries = true
compute_kappa = false
write_cell_heat_flux = true
convergence_write_interval = 100
progress_temperature_summary_only = true

{collision}
[boundary]
boundary_position = {json.dumps(selectors)}
boundary_conditions = {json.dumps(conditions)}
boundary_values = {json.dumps(values)}
periodic_pair = {json.dumps(periodic_pairs)}

{source_text}
[io]
material_folder = "{os.path.relpath(material,out)}"
output_folder = "results/{slug}_{model}"
'''
            target=out/f'{slug}_{model}.toml';target.write_text(text)
            cfg=tomllib.loads(text)
            assert cfg['heat_source']['total_power']==args.power_w
        report['devices'][name]=dict(boxes=parts,mesh=check,mesh_sha256=sha(mesh),
                channel_area_nm2=sum((b[3]-b[2])*(b[5]-b[4]) for b in channels),
                channel_volume_nm3=sum(map(box_volume,channels)),
                layer_power_W=[args.power_w/len(channels)]*len(channels),source=tomllib.loads(source_text)['heat_source'])
    (out/'model.json').write_text(json.dumps(report,indent=2)+'\n')
    plot_model(out,devices,size,args.power_w)
    print(json.dumps({'output':str(out),'sheets':count,'sheet_height_nm':thickness,'gap_nm':gap,'total_power_W':args.power_w},indent=2))


if __name__=='__main__':main()
