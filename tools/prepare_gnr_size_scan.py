#!/usr/bin/env python3
"""Prepare (never submit) matched finite-width GNR scans with no isotopes."""
import json, math, shutil, hashlib, tomllib
from pathlib import Path
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'calculation/gnr_size_scan_20261017'
assert not OUT.exists(), 'Preserve existing campaign; choose a new version'
OUT.mkdir()
for s in ['scripts','validation','logs','analysis','results']:(OUT/s).mkdir()
shutil.copytree(ROOT/'material/Graphene_ph3py_m31311_T50-800',OUT/'material/Graphene')
shutil.copytree(ROOT/'calculation/fet_comparison_20261016/source',OUT/'source')
shutil.copy2(ROOT/'calculation/fet_comparison_20261016/scripts/build.sh',OUT/'scripts/build.sh')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(path,d):path.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
with h5py.File(OUT/'material/Graphene/kappa-fbz.hdf5') as f:
    assert not any(k in f for k in ['gamma_isotope','gamma_impurity','gamma_defect'])
    np.testing.assert_allclose(f['gamma_N'][:]+f['gamma_U'][:],f['gamma'][:],rtol=1e-7,atol=1e-14)
    vmax=float(np.max(np.linalg.norm(f['group_velocity'][:],axis=-1)))/10
    assert np.max(abs(f['group_velocity'][:,:,2]))==0
sizes=[10,50,100,500,1000,5000]
axes=sorted(set([(5000,w) for w in sizes]+[(l,1000) for l in sizes]))
cross=sorted(set((l,w) for l in [100,1000,5000] for w in [100,1000,5000])-set(axes))
cases=[]
def add(stage,L,W,model,seed=12345,delta=2,kind=None):
    nx=max(20,math.ceil(L/50));ny=max(20,math.ceil(W/50));grid=[nx,ny,1]
    dx,dy=L/nx,W/ny
    dt=max(t for t in [.001,.0025,.005,.01,.025,.05,.1] if vmax*t<=.5*min(dx,dy))
    particles=200*nx*ny
    end=500 if L<=100 else 2000 if L<=1000 else 10000
    if stage=='smoke':end=100*dt
    if stage=='pilot':end=200
    tag=f'L{L}nm_W{W}nm_{model}_dT{delta}_s{seed}'
    interval=10 if stage=='smoke' else round(10/dt)
    c=dict(id=tag,stage=stage,length_nm=L,width_nm=W,volume_nm3=L*W,
        model=model,seed=seed,delta_T_K=delta,temperature_K=300,roughness_nm=1,
        simulation_thickness_nm=1,effective_thickness_nm=.335,isotope=False,
        grid_xyz=grid,particles=particles,carriers_per_cell_initial=200,
        dt_ps=dt,iterations=round(end/dt),end_ps=end,average_start_ps=end/2,
        write_interval=interval,cell_traversal_ratio=vmax*dt/min(dx,dy),
        duration_is_provisional=True,classification=kind or stage,status='prepared_not_run')
    folder=OUT/'inputs'/stage;folder.mkdir(parents=True,exist_ok=True)
    text=f'''[geometry]
model = "box"
sizes = [{L}, {W}, 1]
merge_coplanar_facets = false
[simulation]
particle_count = {particles}
time_step = {dt}
iterations = {c['iterations']}
random_seed = {seed}
compute_kappa = false
write_cell_heat_flux = true
transport_axis = "x"
initial_temperature = "linear"
background_temperature = 300
lifetime_temperature = "local"
conservative_boundaries = true
temperature_lookup_dt = 0.1
convergence_write_interval = {interval}
profile_timers = true
progress_temperature_summary_only = true
grid_xyz = {grid}
[scattering]
model = "{model}"
'''
    if model=='callaway':text+='callaway_formulation = "nonlinear"\n'
    text+=f'''[boundary]
boundary_position = [[-0.01,0,0,0.01,1,1],[0.99,0,0,1.01,1,1],[0,-0.01,0,1,0.01,1],[0,0.99,0,1,1.01,1],[0,0,-0.01,1,1,0.01],[0,0,0.99,1,1,1.01]]
boundary_conditions = ["T","T","R","R","P","P"]
boundary_values = [{300+delta/2},{300-delta/2},1,1,0,0]
periodic_pair = [[0,0,-0.01,1,1,0.01],[0,0,0.99,1,1,1.01]]
[io]
material_folder = "../../material/Graphene"
output_folder = "../../results/{stage}/{tag}"
'''
    path=folder/(tag+'.toml');path.write_text(text)
    c['input']=str(path.relative_to(OUT));c['input_sha256']=sha(path);cases.append(c)
for L,W in [(10,10),(5000,10),(5000,5000)]:
    for m in ['rta','callaway']:
        for delta in [0,2]:add('smoke',L,W,m,delta=delta)
for L,W in [(100,100),(5000,1000)]:
    for m in ['rta','callaway']:add('pilot',L,W,m)
for stage,dims in [('axes',axes),('cross_checks',cross)]:
    for L,W in dims:
        for seed in [12345,12346,12347]:
            for m in ['rta','callaway']:add(stage,L,W,m,seed)
pairs=[]
for c in cases:
    if c['model']!='rta':continue
    other=next(o for o in cases if o['stage']==c['stage'] and o['id']==c['id'].replace('_rta_','_callaway_'))
    configs=[]
    for item in [c,other]:
        cfg=tomllib.loads((OUT/item['input']).read_text());cfg.pop('scattering');cfg['io'].pop('output_folder');configs.append(cfg)
    assert configs[0]==configs[1]
    pairs.append(dict(stage=c['stage'],rta=c['id'],callaway=other['id'],match_except_collision=True))
save(OUT/'cases.json',cases);save(OUT/'validation/pairing.json',pairs)
save(OUT/'stage_readiness.json',{s:dict(ready=False,reason='Prepared only; real-material smoke, resource and convergence review required') for s in ['smoke','pilot','axes','cross_checks']})
save(OUT/'plan.json',dict(width_scan=dict(fixed_length_nm=5000,widths_nm=sizes),
    length_scan=dict(fixed_width_nm=1000,lengths_nm=sizes),main_unique_geometries=len(axes),
    independent_seeds=[12345,12346,12347],axes_runs=len(axes)*6,
    cross_check_geometries=cross,cross_check_runs=len(cross)*6,
    smoke_runs=12,pilot_runs=4,no_jobs_submitted=True,
    primary_observables=['Q_W','G_W_K','sheet_coefficient_W_K','kappa_effective_at_0.335nm_W_mK'],
    caveats=['Two-terminal effective conductivity includes contact jumps',
             '10 nm is a semiclassical edge-model extrapolation, not atomistic confinement validation',
             '31x31x1 phonon-grid and supplied Cmmm structure not converged/idealized',
             'Nonlinear distributions/local-temperature rates differ from old linearized fixed-T scan',
             'Three seeds do not replace temporal or discretization convergence'],
    next_checks=['zero-bias and biased smoke','200 ps cost pilots (not steady-state acceptance)',
                 'extend representative duration until stationary','half dt; doubled particles; refined grid with matched particle-per-cell control',
                 'DeltaT 1/2/4 K linear-response check','forward/reverse straight-ribbon null check',
                 'q-mesh convergence independent of real-space grid']))
save(OUT/'material_audit.json',dict(no_isotopes=True,N_U_only=True,vmax_nm_per_ps=vmax,
    material_sha256=sha(OUT/'material/Graphene/kappa-fbz.hdf5'),source_material=str(ROOT/'material/Graphene_ph3py_m31311_T50-800')))
print('Prepared',len(cases),'inputs;',len(axes)*6,'axes;',len(cross)*6,'cross checks. No jobs submitted.')
