#!/usr/bin/env python3
"""Summarize FET transients and observed costs; never certify convergence or submit."""
import argparse, json, re, tomllib
from pathlib import Path
import numpy as np

def analyze(root, case):
    cfg=tomllib.loads((root/case['input']).read_text())
    output=((root/case['input']).parent/cfg['io']['output_folder']).resolve()
    folder=Path(str(output)+'_0')
    with (folder/'convergence.txt').open() as f:
        columns=f.readline().lstrip('#').split();a=np.loadtxt(f,ndmin=2)
    assert a.shape[1]==len(columns) and np.isfinite(a).all()
    d={name:a[:,i] for i,name in enumerate(columns)}
    temp=a[:,[i for i,name in enumerate(columns) if name.startswith('T_')]]
    grid=np.genfromtxt(folder/'grid_centers.csv',delimiter=',',names=True)
    assert np.array_equal(grid['index'],np.arange(temp.shape[1]))
    x,y,z=grid['x_nm'],grid['y_nm'],grid['z_nm']
    fin=case['device']=='finfet'
    bounds=[(20,40)] if fin else [(24.5,28.5),(33,37),(41.5,45.5)]
    width=(9,15) if fin else (7,17)
    layers=[(x>10)&(x<42)&(y>width[0])&(y<width[1])&(z>lo)&(z<hi) for lo,hi in bounds]
    channel=np.logical_or.reduce(layers)
    source=channel&(x>36)&(x<42)
    weights=np.exp(-.5*((x[source]-39)/1.5)**2)*grid['volume_nm3'][source]
    assert source.sum()==1440
    series={'peak_rise_K':temp.max(axis=1)-300,
            'source_weighted_rise_K':np.average(temp[:,source],axis=1,weights=weights)-300,
            'channel_mean_rise_K':np.average(temp[:,channel],axis=1,weights=grid['volume_nm3'][channel])-300}
    for i,mask in enumerate(layers):series[f'layer_{i+1}_mean_rise_K']=temp[:,mask].mean(axis=1)-300
    end=float(d['time_ps'][-1]);late=d['time_ps']>end/2
    previous=late&(d['time_ps']<=.75*end);last=d['time_ps']>.75*end
    metrics={}
    for name,v in series.items():
        first=float(v[previous].mean()) if previous.any() else None
        final=float(v[last].mean()) if last.any() else None
        drift=final-first if first is not None and final is not None else None
        metrics[name]={'last_K':float(v[-1]),'last_quarter_mean_K':final,'previous_quarter_mean_K':first,
                       'late_half_drift_K':drift,'late_half_drift_percent':100*drift/abs(float(v[late].mean())) if drift is not None and abs(v[late].mean())>.01 else None}
    receipt=json.loads((root/'validation'/case['stage']/(case['id']+'.json')).read_text())
    resource=root/'validation/resources'/(case['id']+'.time.txt')
    rss=None
    if resource.exists():
        match=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',resource.read_text())
        if match:rss=int(match[1])
    out=root/'analysis'/case['stage'];out.mkdir(parents=True,exist_ok=True)
    np.savetxt(out/(case['id']+'_transient.csv'),np.column_stack([d['time_ps']]+list(series.values())),delimiter=',',
               header=','.join(['time_ps']+list(series)),comments='')
    summary={'case':case['id'],'device':case['device'],'model':case['model'],'end_ps':end,
             'complete':int(d['timestep'][-1])==case['iterations'],'stored_time_points':len(a),
             'late_time_points':int(late.sum()),'metrics':metrics,'steady_state_proven':False,
             'steady_state_note':'Short transient screen only: sparse and correlated samples, one seed; extend duration before certifying steady state.',
             'wall_seconds_solver_and_validation':receipt['wall_seconds'],'max_rss_kib_runner_and_solver':rss,
             'output_bytes':receipt['output_bytes'],'linear_2ns_wall_estimate_seconds':receipt['wall_seconds']*2000/end,
             'cost_note':'Simple extrapolation includes initialization; later nonlinear costs and full-field output can differ.',
             'next_action':'extend_duration_before_discretization_scan',
             'fixed_seed':case['seed'],'numerical_checks_passed':receipt['state']=='passed'}
    (out/(case['id']+'.json')).write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--index',type=int,required=True);p.add_argument('--stage',default='pilot',choices=['smoke','pilot'])
    args=p.parse_args();root=args.root.resolve()
    cases=[c for c in json.loads((root/'cases.json').read_text()) if c['stage']==args.stage]
    analyze(root,cases[args.index])
