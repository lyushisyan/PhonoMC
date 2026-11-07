#!/usr/bin/env python3
"""Prepare staged convergence inputs from verified pilots. Does not submit jobs."""
import copy, json, math, shutil, subprocess, tomllib
from pathlib import Path
from prepare_fet_comparison_campaign import ROOT, replace, save, sha

base=ROOT/'calculation/fet_comparison_20261016'
results=base/'server_results'
out=ROOT/'calculation/fet_convergence_20261017'
assert not out.exists(),'Never overwrite an existing convergence campaign'
pilots=[c for c in json.loads((base/'cases.json').read_text()) if c['stage']=='pilot']
for c in pilots:
    assert json.loads((results/'validation/pilot'/(c['id']+'.json')).read_text())['state']=='passed'
    assert json.loads((results/'validation'/(c['id']+'_local_verification.json')).read_text())['matches_remote_analysis']
out.mkdir()
for name in ['geometry','material','source']:shutil.copytree(base/name,out/name)
for name in ['scripts','validation','logs']:(out/name).mkdir()
shutil.copy2(results/'validation/smoke_gate.json',out/'validation/smoke_gate.json')
shutil.copy2(base/'scripts/build.sh',out/'scripts/build.sh')
# end ps, dt ps, carriers, grid multiplier, seeds, memory GiB, wall h
plans={'time200':(200,.05,10_000_000,1,[12345],8,8),
       'time400':(400,.05,10_000_000,1,[12345],8,16),
       'dt200':(200,.025,10_000_000,1,[12345],8,16),
       'particles200':(200,.05,20_000_000,1,[12345],16,16),
       'grid200':(200,.05,80_000_000,2,[12345],48,64),
       'grid_control200':(200,.05,80_000_000,1,[12345],48,64),
       'seeds200':(200,.05,10_000_000,1,[12346,12347],8,8)}
stages=['smoke']+list(plans)
cases=[]
for stage,(end,dt,particles,refine,seeds,mem,wall) in plans.items():
    (out/'inputs'/stage).mkdir(parents=True)
    for seed in seeds:
        for old in pilots:
            c=copy.deepcopy(old);tag=f"{old['device']}_{old['model']}_P1uW_s{seed}"
            text=(base/old['input']).read_text()
            values={'iterations':str(round(end/dt)),'time_step':str(dt),'particle_count':str(particles),
                    'random_seed':str(seed),'grid_xyz':str([52*refine,24*refine,100*refine]),
                    'convergence_write_interval':str(round(5/dt)),
                    'output_folder':f'"../../results/{stage}/{tag}"'}
            for k,v in values.items():text=replace(text,k,v)
            path=out/'inputs'/stage/(tag+'.toml');path.write_text(text)
            c.update(id=tag,stage=stage,seed=seed,input=str(path.relative_to(out)),input_sha256=sha(path),
                     iterations=round(end/dt),particles=particles,dt_ps=dt,end_ps=end,write_interval=round(5/dt),
                     grid_xyz=[52*refine,24*refine,100*refine],refinement=refine,request_memory_GiB=mem,request_wall_hours=wall,
                     duration_provisional=stage not in ['time200','time400'],status='prepared_not_run')
            pilot=json.loads((results/'analysis/pilot'/(old['id']+'.json')).read_text())
            c['estimated_wall_hours']=pilot['wall_seconds_solver_and_validation']/3600*(end/50)*(.05/dt)*(particles/10_000_000)
            c['estimate_caveat']='Linear extrapolation only; temperature, geometry work, node and communication costs can change.'
            cases.append(c)
# New binaries can run the same real-material smoke suite; the existing server
# binary reuses its already-passed gate, identified by exact SHA256.
(out/'inputs/smoke').mkdir()
for c in json.loads((base/'cases.json').read_text()):
    if c['stage']=='smoke':
        shutil.copy2(base/c['input'],out/c['input']);cases.append(c)
runner=(base/'scripts/run_case.py').read_text()
runner=runner.replace("choices=['smoke','pilot','production']",'choices='+repr(stages))
runner=runner.replace("if args.stage=='production':","if args.stage!='smoke':")
runner=runner.replace("gate=json.loads((root/'production_readiness.json').read_text())","gate=json.loads((root/'stage_readiness.json').read_text())[args.stage]")
(out/'scripts/run_case.py').write_text(runner)
analyzer=(ROOT/'tools/analyze_fet_pilot.py').read_text()
analyzer=analyzer.replace("choices=['smoke','pilot']",'choices='+repr(stages))
analyzer=analyzer.replace('assert source.sum()==1440','assert np.isclose(grid[\'volume_nm3\'][source].sum(),720)')
analyzer=analyzer.replace("root/'validation/resources'/(case['id']+'.time.txt')", "root/'validation/resources'/case['stage']/(case['id']+'.time.txt')")
(out/'scripts/analyze.py').write_text(analyzer)
binary=json.loads((out/'validation/smoke_gate.json').read_text())['binary_sha256']
save(out/'stage_readiness.json',{s:{'ready':s=='time200','binary_sha256':binary,
    'reason':'First time-extension batch prepared; no job submitted' if s=='time200' else 'Review time200 first; adapt all comparison durations to verified steady-state horizon before release'} for s in plans})
save(out/'cases.json',cases)
# Audit true C++ domain/source construction also for the eightfold grid.
code=(ROOT/'tools/validate_three_layer_fets.cpp').read_text()
code=code.replace('const SimulationDomain domain(config);',
    'const SimulationDomain domain(config);\n            require(config.grid.nx%52==0 && config.grid.ny%24==0 && config.grid.nz%100==0,"Nonaligned refinement");\n            const int mult=(config.grid.nx/52)*(config.grid.ny/24)*(config.grid.nz/100);')
code=code.replace('domain.grid_count()==86400','domain.grid_count()==86400*mult')
code=code.replace('domain.grid_volumes()[i],500','domain.grid_volumes()[i],500.0/mult')
code=code.replace('heated==1440','heated==1440*mult')
cpp=out/'scripts/validate_geometry.cpp';cpp.write_text(code)
validator=Path('/tmp/phonomc_validate_convergence')
subprocess.run(['c++','-std=c++17','-O2','-I'+str(base/'source/include'),str(cpp)]+[str(base/'source'/x) for x in ['src/SimulationConfig.cpp','src/SurfaceMesh.cpp','src/SimulationDomain.cpp','src/solver/PrescribedHeatSource.cpp']]+['-o',str(validator)],check=True)
pairing=[]
for stage,(_,_,_,_,seeds,mem,wall) in plans.items():
    for seed in seeds:
        group=[c for c in cases if c['stage']==stage and c['seed']==seed]
        report=out/'validation'/f'{stage}_s{seed}_geometry.json'
        subprocess.run([str(validator),str(report)]+[str(out/c['input']) for c in group],check=True,stdout=subprocess.DEVNULL)
        assert json.loads(report.read_text())['all_passed']
        for i in [0,2]:
            configs=[tomllib.loads((out/c['input']).read_text()) for c in group[i:i+2]]
            for config in configs:config.pop('scattering');config['io'].pop('output_folder')
            assert configs[0]==configs[1]
        pairing.append({'stage':stage,'seed':seed,'geometry_and_source_passed':True,'models_match_except_collision':True})
    n=sum(c['stage']==stage for c in cases)
    slurm=f'''#!/usr/bin/env bash
#SBATCH --job-name=FET-{stage}
#SBATCH --partition=normal,normal2
#SBATCH --account=wanggang
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem={mem}G
#SBATCH --time={wall}:00:00
#SBATCH --array=0-{n-1}%1
#SBATCH --output=logs/{stage}-%A_%a.out
set -eo pipefail
export LC_ALL=C LANG=C
source /public1/home/wanggang/miniforge3/etc/profile.d/conda.sh
conda activate alamode
set -u
cd "$SLURM_SUBMIT_DIR"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 OMP_WAIT_POLICY=PASSIVE
binary="${{FET_BINARY:-/public1/home/wanggang/software/PhonoMC_runs/fet_comparison_20261016/build/PhonoMC}}"
mkdir -p validation/resources/{stage} analysis/{stage}
case_id=$(python3 -c 'import json,os; a=[c for c in json.load(open("cases.json")) if c["stage"]=="{stage}"]; print(a[int(os.environ["SLURM_ARRAY_TASK_ID"])]["id"])')
/usr/bin/time -v -o "validation/resources/{stage}/${{case_id}}.time.txt" python3 -u scripts/run_case.py {stage} --index "$SLURM_ARRAY_TASK_ID" --binary "$binary" --timeout-seconds {wall*3600-600}
python3 scripts/analyze.py --root . --stage {stage} --index "$SLURM_ARRAY_TASK_ID" > "analysis/{stage}/case_${{SLURM_ARRAY_TASK_ID}}.log"
'''
    (out/'scripts'/(stage+'.slurm')).write_text(slurm)
save(out/'validation/pairing.json',pairing)
save(out/'plan.json',{'prepared_date':'2026-10-17','submitted':False,'first_batch':'time200','convergence_input_count':32,
    'restart_policy':'Fresh 300 K runs; no checkpoint restart is assumed. Compare 50/100/200 ps windows within the 200 ps trajectory.',
    'first_batch_estimated_sequential_hours':sum(c['estimated_wall_hours'] for c in cases if c['stage']=='time200'),
    'grid_control':'Compare refined 80M grid against coarse 80M, not only coarse 10M; refine all three axes by two and preserve carrier density relative to baseline.',
    'criteria':'Provisional <1% late-window mean drift plus uncertainty well below the RTA/Callaway effect; temperature-window agreement alone does not establish independent-seed or grid convergence.',
    'sequence':['time200','time400 only if needed; continue doubling until stable','dt200 and particles200 after updating duration to stable horizon','grid_control200 versus grid200','seeds200 plus accepted seed12345 reference'],
    'output_policy':'Temperature every 5 ps; no full vector heat-flux fields in convergence runs. Thermal-boundary power partition still requires separate diagnostics.'})
save(out/'files_sha256.json',{str(p.relative_to(out)):sha(p) for folder in ['source','material','geometry','inputs','scripts'] for p in sorted((out/folder).rglob('*')) if p.is_file()})
shutil.copy2(ROOT/'model/matched_fets/three_layer_comparison/convergence_README.md',out/'README.md')
print(json.dumps({'out':str(out),'convergence_inputs':32,'smoke_inputs':8,'geometry_groups_passed':len(pairing),'submitted':False},indent=2))
