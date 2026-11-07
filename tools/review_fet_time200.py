#!/usr/bin/env python3
"""Verify a time200 export and assess drift without treating blocks as seeds."""
import argparse, hashlib, json, tarfile
from pathlib import Path
import numpy as np

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def save(p,d):
    p.parent.mkdir(parents=True,exist_ok=True)
    text=json.dumps(d,indent=2,allow_nan=False)+'\n'
    if p.exists():assert p.read_text()==text,'Preserve previous review: '+str(p)
    else:p.write_text(text)

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('archive',type=Path);p.add_argument('--sha256',required=True)
p.add_argument('--root',type=Path,required=True);p.add_argument('--pilot-root',type=Path,required=True)
args=p.parse_args();root=args.root.resolve();case=args.archive.name.removesuffix('.tar.gz')
assert sha(args.archive)==args.sha256
with tarfile.open(args.archive) as t:
    for m in t:
        assert m.isfile() and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts
        path=root/m.name;assert path.resolve().is_relative_to(root)
        data=t.extractfile(m).read()
        if path.exists():assert path.read_bytes()==data,'Existing file differs: '+m.name
        else:path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
manifest=json.loads((root/'export_manifests/time200'/(case+'.json')).read_text())
for name,h in manifest['files'].items():assert sha(root/name)==h,name
c=next(c for c in json.loads((root/'cases.json').read_text()) if c['id']==case and c['stage']=='time200')
folder=root/'results/time200'/(case+'_0')
with (folder/'convergence.txt').open() as f:
    cols=f.readline().lstrip('#').split();a=np.loadtxt(f,ndmin=2)
assert a.shape[1]==len(cols) and np.isfinite(a).all()
assert int(a[-1,0])==c['iterations'] and abs(a[-1,1]-200)<1e-7
d={k:a[:,i] for i,k in enumerate(cols)}
temp=a[:,[i for i,k in enumerate(cols) if k.startswith('T_')]]
grid=np.genfromtxt(folder/'grid_centers.csv',delimiter=',',names=True)
assert np.array_equal(grid['index'],np.arange(temp.shape[1]))
x,y,z,v=[grid[k] for k in ['x_nm','y_nm','z_nm','volume_nm3']]
fin=c['device']=='finfet';width=(9,15) if fin else (7,17)
bounds=[(20,40)] if fin else [(24.5,28.5),(33,37),(41.5,45.5)]
layers=[(x>10)&(x<42)&(y>width[0])&(y<width[1])&(z>lo)&(z<hi) for lo,hi in bounds]
channel=np.logical_or.reduce(layers);source=channel&(x>36)&(x<42)
assert np.isclose(v[source].sum(),720)
series={'peak_rise_K':temp.max(axis=1)-300,
        'source_weighted_rise_K':np.average(temp[:,source],axis=1,weights=v[source]*np.exp(-.5*((x[source]-39)/1.5)**2))-300,
        'channel_mean_rise_K':np.average(temp[:,channel],axis=1,weights=v[channel])-300}
for i,mask in enumerate(layers):series[f'layer_{i+1}_mean_rise_K']=np.average(temp[:,mask],axis=1,weights=v[mask])-300
remote=json.loads((root/'analysis/time200'/(case+'.json')).read_text())
old=json.loads((args.pilot_root/'analysis/pilot'/(case+'.json')).read_text())
receipt=json.loads((root/'validation/time200'/(case+'.json')).read_text())
assert receipt['state']=='passed'
for name,s in series.items():assert abs(s[-1]-remote['metrics'][name]['last_K'])<1e-9
for key in ['energy_balance_residual_ev','lifetime_energy_residual_ev']:
    assert np.max(abs(d[key]))<1e-6
assert np.max(abs(d['hs_ledger_residual_ev']))<max(1e-8,float(d['hs_prescribed_total_ev'][-1])*1e-10)
assert np.isclose(d['hs_prescribed_total_ev'][-1],1e-6*200*6.241509074e6,rtol=1e-9,atol=1e-10)
metrics={};t=d['time_ps']
for name,s in series.items():
    late=t>100;tl=t[late];sl=s[late]
    slope,intercept=np.polyfit(tl,sl,1);res=sl-(slope*tl+intercept)
    rms=float(np.sqrt(np.mean(res**2)))
    blocks=[float(s[(t>lo)&(t<=lo+25)].mean()) for lo in [100,125,150,175]]
    mean=float(sl.mean());drift=float(s[t>150].mean()-s[(t>100)&(t<=150)].mean())
    metrics[name]=dict(mean_100_200_K=mean,mean_150_200_K=float(s[t>150].mean()),
                       mean_last_quarter_old50ps_K=old['metrics'][name]['last_quarter_mean_K'],
                       old_comparison_window_note='Old value is the last quarter, 37.5–50 ps, sampled every 5 ps.',
                       rise_change_from_old_last_quarter_K=float(s[t>150].mean()-old['metrics'][name]['last_quarter_mean_K']),
                       four_25ps_block_means_K=blocks,half_window_drift_K=drift,
                       half_window_drift_percent=100*drift/abs(mean),slope_100_200_K_per_ps=float(slope),
                       fitted_change_100ps_K=float(100*slope),detrended_temporal_rms_K=rms,
                       fitted_change_to_temporal_rms=float(abs(100*slope)/rms) if rms>0 else None,
                       slope_150_200_K_per_ps=float(np.polyfit(t[t>150],s[t>150],1)[0]),
                       temporal_blocks_are_independent_seeds=False)
review=dict(case=case,archive_sha256=args.sha256,verified_files=len(manifest['files']),raw_steps=int(a[-1,0]),
            numerical_checks=receipt['metrics'],metrics=metrics,steady_state_certified=False,
            interpretation='Assess trends against correlated temporal scatter; this single-seed screen does not certify time or discretization convergence.',
            wall_seconds_measured=receipt['wall_seconds'],peak_rss_kib=remote['max_rss_kib_runner_and_solver'],
            fresh_400ps_wall_hours_linear_estimate=receipt['wall_seconds']/3600*2,
            estimate_caveat='Linear extrapolation, not measured 400 ps cost; no 400 ps submission authorized by this review.')
save(root/'validation/local_verification'/(case+'.json'),review)
print(json.dumps(review,indent=2))
