#!/usr/bin/env python3
"""Run a single frozen FET input, or the short smoke suite. Never calls sbatch."""
from __future__ import annotations
import argparse, hashlib, json, math, os, subprocess, time, tomllib
from pathlib import Path
import numpy as np

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def analyze(folder,c,log):
    with (folder/'convergence.txt').open() as f:
        cols=f.readline().lstrip('#').split();a=np.loadtxt(f,ndmin=2)
    assert a.shape[1]==len(cols) and np.isfinite(a).all(),'Malformed/nonfinite output'
    d={k:a[:,i] for i,k in enumerate(cols)}
    assert int(d['timestep'][-1])==c['iterations'],'Incomplete trajectory'
    summary={}
    for line in (folder/'summary.txt').read_text().splitlines():
        if ' = ' in line:
            k,v=line.split(' = ',1);summary[k]=v
    for k in ['collision_cache_failures_total','excessive_collision_recoveries_total']:
        assert int(summary[k])==0,(k,summary[k])
    assert 'reservoir_gen=fixed_incoming_flux' in log,'Incorrect bath implementation'
    temps=a[:,[i for i,k in enumerate(cols) if k.startswith('T_')]]
    assert temps.min()>=50 and temps.max()<=800,'Outside material rate table'
    er=float(np.max(abs(d['energy_balance_residual_ev'])))
    cr=float(np.max(abs(d['lifetime_energy_residual_ev'])))
    sr=float(np.max(abs(d['hs_ledger_residual_ev'])))
    total=float(d['hs_prescribed_total_ev'][-1]);pending=float(d['hs_pending_total_ev'][-1])
    target=c['power_uW']*1e-6*c['end_ps']*6.241509074e6
    assert math.isclose(total,target,rel_tol=1e-9,abs_tol=1e-10),'Integrated source power differs'
    assert er<1e-6 and cr<1e-6 and sr<max(1e-8,total*1e-10),'Energy residual exceeds short-run screen'
    # Empty-cell energy is retained, not silently discarded. Flag it for pilot assessment.
    pending_fraction=pending/max(total,1e-30)
    if c['power_uW']==0:assert np.max(abs(temps-300))<1e-7,'Zero-source equilibrium drift'
    else:assert float(summary['injected_energy_total_ev'])>0,'No heat deposited'
    if c['power_uW']>0:assert pending_fraction<.01,'More than 1% source energy pending'
    fallback_keys=['fallback_missing_rough_data','fallback_missing_spec_match','fallback_outgoing_pool','fallback_global_random','residual_flux_frequency_window_fallback']
    fallbacks={k:int(summary.get(k,0)) for k in fallback_keys}
    assert not any(fallbacks.values()),('Rough-boundary fallback',fallbacks)
    return dict(final_step=int(d['timestep'][-1]),temperature_min_K=float(temps.min()),temperature_max_K=float(temps.max()),
                max_energy_residual_ev=er,max_collision_residual_ev=cr,max_source_ledger_residual_ev=sr,
                prescribed_energy_ev=total,expected_prescribed_energy_ev=target,pending_energy_fraction=pending_fraction,
                particle_count_min=float(d['particle_count'].min()),particle_count_max=float(d['particle_count'].max()),
                rough_fallbacks=fallbacks,transport_converged=False,
                checks='Short-run operability and accounting only; no steady-state or discretization claim')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['smoke','pilot','production'])
    p.add_argument('--index',type=int)
    p.add_argument('--all',action='store_true')
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--threads',type=int,default=int(os.environ.get('OMP_NUM_THREADS','4')))
    p.add_argument('--timeout-seconds',type=float,default=None)
    args=p.parse_args();root=args.root.resolve();binary=args.binary.resolve()
    assert binary.is_file();bh=sha(binary)
    for name,expected in json.loads((root/'files_sha256.json').read_text()).items():
        assert sha(root/name)==expected,'Frozen file changed: '+name
    cases=json.loads((root/'cases.json').read_text());selected=[c for c in cases if c['stage']==args.stage]
    assert args.all != (args.index is not None),'Choose --all or --index'
    if args.all:assert args.stage=='smoke','Only short smokes can run as a local batch'
    else:selected=[selected[args.index]]
    if args.stage!='smoke':
        gate=json.loads((root/'validation/smoke_gate.json').read_text())
        assert gate['passed'] and gate['binary_sha256']==bh,'Run smokes with this exact binary first'
    if args.stage=='production':
        gate=json.loads((root/'production_readiness.json').read_text())
        assert gate['ready'] and gate.get('binary_sha256')==bh,'Pilot/resource/scientific readiness not established'
    for c in selected:
        receipt=root/'validation'/args.stage/(c['id']+'.json')
        if receipt.exists():raise RuntimeError('Existing receipt; preserve it and use a fresh campaign for reruns: '+str(receipt))
        cfg=tomllib.loads((root/c['input']).read_text())
        prefix=(root/c['input']).parent/cfg['io']['output_folder'];folder=Path(str(prefix.resolve())+'_0')
        assert not folder.exists(),'Never overwrite old trajectory'
        logpath=root/'logs'/(args.stage+'_'+c['id']+'.log')
        state=dict(case=c['id'],stage=args.stage,input_sha256=sha(root/c['input']),binary_sha256=bh,
                   threads=args.threads,started_unix=time.time(),state='running')
        save(receipt,state);start=time.monotonic()
        env=dict(os.environ,OMP_NUM_THREADS=str(args.threads),OMP_DYNAMIC='FALSE',OMP_WAIT_POLICY='PASSIVE',OPENBLAS_NUM_THREADS='1')
        env.pop('PHONOMC_ALLOW_SYNTHETIC_MATERIAL',None)
        try:
            print('START',args.stage,c['id'],flush=True)
            with logpath.open('w') as f:
                # MaterialFactory searches cwd before the input directory. Use
                # the input directory to prevent another material/Si winning.
                proc=subprocess.run([str(binary),str(root/c['input'])],cwd=(root/c['input']).parent,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=args.timeout_seconds)
            state['returncode']=proc.returncode
            if proc.returncode:raise RuntimeError('Solver returned '+str(proc.returncode)+'; see '+str(logpath))
            actual_log=logpath.read_text()
            expected_material=((root/c['input']).parent/cfg['io']['material_folder']).resolve()
            assert ('PhononMaterial material folder: '+str(expected_material)+'\n') in actual_log,'Solver loaded an unexpected material directory'
            state['metrics']=analyze(folder,c,actual_log)
            state['state']='passed';state['output_bytes']=sum(f.stat().st_size for f in folder.rglob('*') if f.is_file())
        except Exception as e:
            state['state']='failed';state['error']=str(e)
            raise
        finally:
            state['wall_seconds']=time.monotonic()-start;save(receipt,state)
        print('PASS',c['id'],round(state['wall_seconds'],2),'s',flush=True)
    if args.stage=='smoke':
        expected=[c for c in cases if c['stage']=='smoke'];receipts=[]
        for c in expected:
            f=root/'validation/smoke'/(c['id']+'.json')
            if f.exists():receipts.append(json.loads(f.read_text()))
        passed=len(receipts)==len(expected) and all(r['state']=='passed' and r['binary_sha256']==bh for r in receipts)
        save(root/'validation/smoke_gate.json',dict(passed=passed,completed=len(receipts),expected=len(expected),binary_sha256=bh,
             note='Short real-material runs only. Server binary must pass its own smoke suite; production remains gated.'))

if __name__=='__main__':main()
