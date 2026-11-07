#!/usr/bin/env python3
"""Prepare a portable, frozen FinFET/three-sheet GAAFET comparison. No submission."""
from __future__ import annotations
import argparse, hashlib, json, math, re, shutil, subprocess, tomllib
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'model/matched_fets/three_layer_comparison'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def replace(text,key,value):
    text,n=re.subn(r'^'+re.escape(key)+r'\s*=.*$',key+' = '+value,text,flags=re.M)
    assert n==1,(key,n)
    return text

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=ROOT/'calculation/fet_comparison_20261016')
    ap.add_argument('--validator',type=Path,default=Path('/tmp/phonomc_validate_three_layer_fets'))
    args=ap.parse_args();out=args.out.resolve()
    if out.exists():raise SystemExit('Refusing to overwrite existing campaign: '+str(out))
    for directory in ['geometry','material/Si','source','scripts','inputs/smoke','inputs/pilot','inputs/production','validation','logs']:
        (out/directory).mkdir(parents=True,exist_ok=True)
    for device in ['finfet','gaafet']:shutil.copy2(BASE/(device+'.stl'),out/'geometry'/(device+'.stl'))
    shutil.copy2(BASE/'model.json',out/'geometry/model.json')
    for name in ['POSCAR','kappa-fbz.hdf5']:shutil.copy2(ROOT/'material/Si_ph3py_m252525_T50-800'/name,out/'material/Si'/name)
    for name in ['src','include']:shutil.copytree(ROOT/name,out/'source'/name)
    shutil.copy2(ROOT/'CMakeLists.txt',out/'source/CMakeLists.txt')
    shutil.copy2(ROOT/'tools/validate_three_layer_fets.cpp',out/'scripts/validate_geometry.cpp')
    shutil.copy2(ROOT/'tools/run_fet_comparison_case.py',out/'scripts/run_case.py')
    cases=[]
    plans=[('smoke',[0.,1.],[12345],500000,5,1),
           ('pilot',[1.],[12345],10000000,1000,100),
           ('production',[.2,1.,2.],[12345,12346,12347],10000000,40000,200)]
    for stage,powers,seeds,particles,steps,interval in plans:
        for power in powers:
            for seed in seeds:
                for device in ['finfet','gaafet']:
                    for model in ['rta','callaway']:
                        tag=f'{device}_{model}_P{power:g}uW_s{seed}'
                        text=(BASE/f'{device}_{model}.toml').read_text()
                        text=text.replace(f'model = "{device}.stl"',f'model = "../../geometry/{device}.stl"')
                        for k,v in {'particle_count':str(particles),'iterations':str(steps),'random_seed':str(seed),
                                    'convergence_write_interval':str(interval),
                                    'material_folder':'"../../material/Si"','output_folder':f'"../../results/{stage}/{tag}"'}.items():
                            text=replace(text,k,v)
                        # Keep the source support identical at every power; zero-source smoke disables it.
                        text=replace(text,'total_power',f'{(power if power else 1.)*1e-6:.15g}')
                        if power==0:text=replace(text,'enabled','false')
                        flux=stage=='smoke' or (stage=='production' and power==1 and seed==12345)
                        text=replace(text,'write_cell_heat_flux',str(flux).lower())
                        path=out/'inputs'/stage/(tag+'.toml');path.write_text(text)
                        cfg=tomllib.loads(text)
                        assert cfg['simulation']['conservative_boundaries'] is True
                        assert cfg['simulation']['lifetime_temperature']=='local'
                        assert cfg['boundary']['boundary_conditions']==['R','P','P','P','P','T','T','T']
                        assert math.isclose(cfg['heat_source']['total_power'],(power if power else 1.)*1e-6,rel_tol=1e-14)
                        cases.append(dict(id=tag,stage=stage,device=device,model=model,power_uW=power,seed=seed,
                            input=str(path.relative_to(out)),input_sha256=sha(path),particles=particles,iterations=steps,dt_ps=.05,
                            end_ps=.05*steps,write_interval=interval,cell_heat_flux=flux,
                            status='prepared_not_run',power_selection='provisional_pending_1uW_pilot' if stage=='production' else 'fixed'))
    # Compare models with everything except the collision operator and result path held fixed.
    paired=[]
    for c in cases:
        if c['model']!='rta':continue
        b=next(x for x in cases if all(x[k]==c[k] for k in ['stage','device','power_uW','seed']) and x['model']=='callaway')
        ca=tomllib.loads((out/c['input']).read_text());cb=tomllib.loads((out/b['input']).read_text())
        for d in [ca,cb]:d.pop('scattering');d['io'].pop('output_folder')
        assert ca==cb,(c['id'],b['id'])
        paired.append(dict(stage=c['stage'],rta=c['id'],callaway=b['id'],physical_and_numerical_inputs_equal=True))
    # Audit every production group with the actual C++ geometry/source implementation.
    audit_files=[]
    for power in [.2,1.,2.]:
        for seed in [12345,12346,12347]:
            selected=[c for c in cases if c['stage']=='production' and c['power_uW']==power and c['seed']==seed]
            audit=out/'validation'/f'geometry_P{power:g}uW_s{seed}.json'
            subprocess.run([str(args.validator.resolve()),str(audit)]+[str(out/c['input']) for c in selected],check=True,stdout=subprocess.DEVNULL)
            assert json.loads(audit.read_text())['all_passed']
            audit_files.append(str(audit.relative_to(out)))
    save(out/'cases.json',cases);save(out/'validation/model_pairing.json',paired)
    save(out/'campaign.json',dict(title='Matched FinFET and three-sheet GAAFET',
        material='Identical Si 25^3 three-phonon N/U table; 50-800 K; local-temperature RTA / nonlinear Callaway',
        device_envelope_nm=[52,24,50],substrate_nm=[52,24,20],channel_area_nm2=120,channel_volume_nm3=3840,
        gaafet_layers=3,gaafet_layer_power_fraction=[1/3]*3,roughness_nm=1,baths_K=300,
        periodic_substrate_only=True,total_power_is_per_device=True,
        source='x=36..42 nm Gaussian centered at 39 nm, sigma=1.5 nm; uniform over actual channel cross sections',
        counts={s:sum(c['stage']==s for c in cases) for s in ['smoke','pilot','production']},
        geometry_audits=audit_files,transport_convergence_verified=False,submitted=False,
        output_policy='Temperature histories every 10 ps for all production cases; vector heat-flux fields only for four 1 uW seed-12345 cases.',
        estimated_production_output_GiB_range=[12,22],
        estimate_note='Planning estimate only, excludes filesystem replication and longer runs; measure with pilot before scheduling.',
        per_bath_power_partition_available=False,
        per_bath_note='Current output has total reservoir exchange per sampled step, not complete per-bath cumulative energy. Do not infer bath fractions from sparse samples.',
        remaining=['real-material short smokes','server pilot timing/RSS/output calibration','steady-state/grid/time-step/particle/seed study','per-bath energy diagnostics if Figure 6 requires bath-resolved powers']))
    save(out/'production_readiness.json',dict(ready=False,reasons=['Server pilots not run','Power range provisional','2 ns and 10M are initial settings, not converged'],requires={'smoke_gate':True,'pilot_gate':True}))
    build='''#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
cmake -S source -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build -j "${BUILD_JOBS:-4}"
'''
    (out/'scripts/build.sh').write_text(build)
    # Slurm wrappers are examples for this account; they do not submit any jobs.
    for stage,n in [('smoke',8),('pilot',4),('production',36)]:
        limit=1
        script=f'''#!/usr/bin/env bash
#SBATCH --job-name=FET-{stage}
#SBATCH --partition=normal
#SBATCH --account=wanggang
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --array=0-{n-1}%{limit}
#SBATCH --output=logs/{stage}-%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK:-28}}" OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 OMP_WAIT_POLICY=PASSIVE
python3 scripts/run_case.py {stage} --index "$SLURM_ARRAY_TASK_ID" --binary build/PhonoMC
'''
        (out/'scripts'/f'{stage}.slurm').write_text(script)
    filemap={str(p.relative_to(out)):sha(p) for folder in ['source','geometry','material','inputs','scripts'] for p in sorted((out/folder).rglob('*')) if p.is_file()}
    save(out/'files_sha256.json',filemap)
    shutil.copy2(BASE/'campaign_README.md',out/'README.md')
    print(json.dumps(dict(output=str(out),smoke=8,pilot=4,production=36,geometry_groups_verified=len(audit_files),submitted=False),indent=2))

if __name__=='__main__':main()
