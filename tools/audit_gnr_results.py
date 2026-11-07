#!/usr/bin/env python3
"""Read-only GNR evidence audit; never relabel a snapshot as final results."""
import json, hashlib, collections, html, tomllib
from pathlib import Path
import numpy as np
import h5py

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/gnr_audit_20261017'
OUT.mkdir(parents=True,exist_ok=True)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data): (OUT/name).write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

snapshot=ROOT/'output/gnr_300K_width_comparison/snapshot/results.json'
s=json.loads(snapshot.read_text()); done=s['completed']; paired=collections.defaultdict(dict)
for c in done: paired[c['length_nm'],c['width_nm']][c['model']]=c
finite=dict(expected_cases=s['expected_cases'],snapshot_completed=len(done),snapshot_failed=len(s['failed']),
    snapshot_flagged=sum(bool(c['review_required']) for c in done),
    complete_model_pairs=sum(set(x)=={'rta','callaway'} for x in paired.values()),
    paired_screen_pass=sum(set(x)=={'rta','callaway'} and not any(v['review_required'] for v in x.values()) for x in paired.values()),
    lengths_present_nm=sorted(set(c['length_nm'] for c in done)),
    final_carriers_range=[min(c['final_carriers'] for c in done),max(c['final_carriers'] for c in done)],
    wall_seconds_range=[min(c['wall_seconds'] for c in done),max(c['wall_seconds'] for c in done)],
    snapshot_path=str(snapshot),snapshot_sha256=sha(snapshot),raw_recomputed=False,
    reason='Local production trajectories absent; this is an incomplete historical snapshot, not current campaign completion.')
save('finite_width_snapshot.json',dict(summary=finite,cases=done))

base=ROOT/'output/graphene_length_scan_100_300K_20260910'
records=[]
for old in json.loads((base/'mc_results.json').read_text()):
    name=old['case']; cfg=tomllib.loads((base/(name+'.toml')).read_text())
    path=base/(name+'_0')/'convergence.txt'
    with path.open() as f:
        columns=f.readline().lstrip('#').split(); a=np.loadtxt(f,ndmin=2)
    assert a.shape[1]==len(columns) and np.isfinite(a).all()
    d={k:a[:,i] for i,k in enumerate(columns)}
    mask=(d['time_ps']>old['analysis_start_ps'])&(d['time_ps']<=old['end_ps']+1e-8)
    k=d['kappa_eff'][mask]; mean=float(k.mean()); mid=len(k)//2
    delta=abs(mean-old['kappa_raw_1nm_W_mK'])/max(abs(mean),1e-30)
    records.append(dict(case=name,model=old['model'],temperature_K=old['temperature_K'],length_nm=old['length_nm'],
        raw_complete=int(d['timestep'][-1])==cfg['simulation']['iterations'],
        actual_end_ps=float(d['time_ps'][-1]),kappa_0335nm_W_mK=mean/.335,
        summary_relative_difference=delta,summary_matches=delta<1e-7,
        late_half_change_percent=100*float(k[mid:].mean()-k[:mid].mean())/abs(mean),
        max_energy_residual_ev=float(np.max(abs(d['energy_balance_residual_ev']))),
        max_collision_residual_ev=float(np.max(abs(d['lifetime_energy_residual_ev']))),
        y_periodic=True,finite_width_nanorribbon=False,
        input_sha256=sha(base/(name+'.toml')),trajectory_sha256=sha(path)))
save('periodic_sheet_raw_audit.json',records)

material=ROOT/'material/Graphene_ph3py_m31311_T50-800'
with h5py.File(material/'kappa-fbz.hdf5') as f:
    keys=list(f);nu=float(np.max(abs(f['gamma_N'][:]+f['gamma_U'][:]-f['gamma'][:])))
    extras={k:float(np.max(abs(f[k][:]))) for k in keys if k in ['gamma_isotope','gamma_impurity','gamma_defect']}
    ma=dict(path=str(material),sha256=sha(material/'kappa-fbz.hdf5'),mesh=f['mesh'][:].tolist(),
        frequency_shape=list(f['frequency'].shape),temperature_K=f['temperature'][:].tolist(),
        N_plus_U_max_error_THz=nu,extra_scattering=extras,isotope_present='gamma_isotope' in keys,
        max_abs_vz_A_per_ps=float(np.max(abs(f['group_velocity'][:,:,2]))),
        vmax_nm_per_ps=float(np.max(np.linalg.norm(f['group_velocity'][:],axis=-1)))/10)
    assert not extras and nu<1e-12
save('material_audit.json',ma)
missing=[]
for name in ['gnr_rectification_20260913','gnr_width_2K_20260913','gnr_20K_W1000_20260914']:
    root=ROOT/'calculation'/name; cases=json.loads((root/'cases.json').read_text())
    available=[c['case'] for c in cases if (root/'results'/(c['case']+'_0')/'convergence.txt').exists()]
    missing.append(dict(campaign=name,configured_cases=len(cases),production_raw_present=len(available),available_cases=available))
summary=dict(finite_width=finite,periodic_sheet=dict(raw_traces=len(records),complete=sum(r['raw_complete'] for r in records),
    summary_matches=sum(r['summary_matches'] for r in records),over_1percent_drift=sum(abs(r['late_half_change_percent'])>1 for r in records),
    interpretation='Transverse-periodic sheet length dependence; not finite-width ribbon validation. Independent reference uses the same approximate collision model, not full three-phonon LBTE.'),
    nonlinear_campaigns=missing,material=ma,T7_mounted=Path('/Volumes/T7 Shield').exists(),
    finite_width_steady_state_certified=0,
    required=['Recover finite-width raw trajectories and final campaign receipts from T7',
              'Keep sheet and finite-width ribbon results separate',
              'Use identical boundary kernel, source material and thickness for model pairs',
              'Check time, grid, timestep, carrier and independent-seed convergence',
              'Validate q-grid and supplied Cmmm-identified slab before ideal atomistic-ribbon claims'])
save('summary.json',summary)
rows=[('有限宽扫描快照',f"{len(done)}/72 条；{finite['snapshot_flagged']} 条带复核标记；原始生产轨迹待归档核对"),
      ('横向周期片层',f"{len(records)} 条原始轨迹已重算，{sum(r['summary_matches'] for r in records)} 条与汇总一致；不能作有限宽纳米带曲线"),
      ('非线性整流/宽度系列','本地只有配置及开发短测，生产轨迹待 T7'),
      ('材料','31×31×1，三声子 N/U，无同位素；模拟高度 1 nm，展示厚度 0.335 nm')]
page='<!doctype html><meta charset="utf-8"><title>石墨烯纳米带结果核查</title><style>body{font:16px system-ui;max-width:1100px;margin:40px auto;color:#20344a}td,th{border:1px solid #ccd;padding:12px}table{border-collapse:collapse;width:100%}pre{white-space:pre-wrap}</style><h1>石墨烯纳米带结果核查</h1><p>2026-10-17。只读核查；快照完成数量不代表当前最终完成状态。</p><table>'+''.join('<tr><th>'+html.escape(k)+'</th><td>'+html.escape(v)+'</td></tr>' for k,v in rows)+'</table><h2>后续核验</h2><ul>'+''.join('<li>'+html.escape(x)+'</li>' for x in summary['required'])+'</ul>'
(OUT/'report.html').write_text(page)
print(json.dumps({k:v for k,v in summary.items() if k not in ['material','required']},ensure_ascii=False,indent=2))
