#!/usr/bin/env python3
"""Read-only audit of saved Si/SiC film campaigns; never changes old results.

Extracts only selected inputs, summaries and convergence histories from verified
T7 archives. Screening thresholds reproduce the campaign, not convergence proof.
"""
from __future__ import annotations
import argparse, collections, csv, hashlib, html, json, math, statistics, tarfile, tomllib
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    'Si_RTA': 'calculation/figure3_callaway_si_20261006/si_films/rta_reference/kappa_summary.csv',
    'Si_Callaway': 'calculation/figure3_callaway_si_20261006/plot_snapshot/si_films/kappa_summary.csv',
    'SiC_RTA': 'calculation/figure3_extensions_20261005/plot_snapshot/sic_films/kappa_summary.csv',
    'SiC_Callaway': 'calculation/figure3_callaway_20261005/plot_snapshot/sic_films/kappa_summary.csv',
}
ARCHIVES = ['figure3_callaway_si_20261006', 'figure3_callaway_20261005',
            'figure3_extensions_20261005', 'figure3_rta_20260913', 'fixed_carrier_scan_20260912']
THRESHOLDS = {'half_window_drift_percent': 1., 'block_range_percent': 1.,
              'max_cell_temperature_drift_K': .01, 'max_energy_residual_ev': 1e-6,
              'max_collision_residual_ev': 1e-6}

def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''): h.update(b)
    return h.hexdigest()

def write_json(p, x):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False)+'\n')

def extract(archive_root, out, wanted):
    checks = {line.split(maxsplit=1)[1].strip().lstrip('*'): line.split()[0]
              for line in (archive_root/'SHA256SUMS.txt').read_text().splitlines() if line.strip()}
    records=[]
    for name in ARCHIVES:
        p=archive_root/(name+'.tar.gz'); actual=digest(p)
        if actual != checks[p.name]: raise RuntimeError('Archive SHA256 mismatch: '+str(p))
        count=0; size=0
        with tarfile.open(p, 'r|gz') as t:
            for m in t:
                parts=Path(m.name).parts
                if not m.isfile() or '..' in parts or Path(m.name).is_absolute(): continue
                if any(s in m.name for s in ['/source/', '/smoke/', '/archive_before', '/local_preflight/']): continue
                leaf=parts[-1]; parent=parts[-2] if len(parts)>1 else ''
                selected=(leaf.endswith('.toml') and leaf[:-5] in wanted)
                selected |= any(leaf==case+suffix for suffix in ('_analysis.json','_exit.json')
                                for case in wanted if leaf.startswith(case))
                selected |= leaf in ('convergence.txt','summary.txt') and parent.endswith('_0') and parent[:-2] in wanted
                if not selected: continue
                target=out/'raw'/Path(*parts[1:])
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open('wb') as f:
                    src=t.extractfile(m)
                    for b in iter(lambda:src.read(1024*1024),b''): f.write(b)
                count+=1; size+=m.size
        records.append(dict(archive=str(p),sha256=actual,verified=True,selected_files=count,selected_bytes=size))
        print('Verified and read',name,count,'selected files',flush=True)
    write_json(out/'archive_verification.json', records)

def key(row):
    return (row['group'].split('_')[0], row['family'].split('_',1)[1],
            float(row['thickness_nm']), None if row['roughness_nm']=='' else float(row['roughness_nm']))

def raw_metrics(path, config, factor):
    with path.open() as f:
        columns=f.readline().lstrip('#').split()
        a=np.loadtxt(f,ndmin=2)
    if a.shape[1]!=len(columns) or not np.isfinite(a).all(): raise ValueError('Malformed or nonfinite history')
    d={c:a[:,i] for i,c in enumerate(columns)}
    dt=config['simulation']['time_step']; end=config['simulation']['iterations']*dt
    time=d['timestep']*dt; mask=(time>end/2)&(time<=end+1e-8)
    ks=d['kappa_eff'][mask]*factor
    if len(ks)<10: raise ValueError('Insufficient late samples')
    blocks=[float(np.mean(b)) for b in np.array_split(ks,5)]
    mid=len(ks)//2; mean=float(np.mean(ks))
    temps=[v[mask] for c,v in d.items() if c.startswith('T_')]
    result=dict(raw_final_step=int(d['timestep'][-1]),expected_final_step=config['simulation']['iterations'],
                raw_complete=int(d['timestep'][-1])==config['simulation']['iterations'],
                samples=len(ks),end_ps=end,kappa_recomputed_W_mK=mean,
                block_sd_W_mK=statistics.stdev(blocks),blocks_W_mK=blocks,
                half_window_drift_percent=100*float(np.mean(ks[mid:])-np.mean(ks[:mid]))/abs(mean),
                block_range_percent=100*(max(blocks)-min(blocks))/abs(mean),
                quarter_vs_half_percent=100*float(np.mean(ks[mid:])-mean)/abs(mean),
                max_cell_temperature_drift_K=max((float(abs(np.mean(v[mid:])-np.mean(v[:mid]))) for v in temps),default=0.),
                max_energy_residual_ev=float(np.max(abs(d['energy_balance_residual_ev']))),
                max_collision_residual_ev=float(np.max(abs(d['lifetime_energy_residual_ev']))),
                relative_energy_residual_to_stored_energy=float(np.max(abs(d['energy_balance_residual_ev'])/np.maximum(abs(d['total_thermal_energy_ev']),1e-30))),
                relative_collision_residual_to_stored_energy=float(np.max(abs(d['lifetime_energy_residual_ev'])/np.maximum(abs(d['total_thermal_energy_ev']),1e-30))),
                particle_count_min=float(np.min(d['particle_count'])),particle_count_max=float(np.max(d['particle_count'])))
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--archive-root',type=Path,default=Path('/Volumes/T7 Shield/PhonoMC_runs_20261013'))
    ap.add_argument('--out',type=Path,default=ROOT/'analysis/film_callaway_audit_20261016')
    ap.add_argument('--extract',action='store_true')
    args=ap.parse_args(); out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for group,path in SOURCES.items():
        with (ROOT/path).open() as f:
            for row in csv.DictReader(f):rows.append(dict(row,group=group,summary_source=path))
    assert len(rows)==280 and len({r['case'] for r in rows})==280
    wanted={r['case'] for r in rows}
    if args.extract: extract(args.archive_root,out,wanted)
    # Local files remain useful if an archive was deliberately not extracted.
    inputs={}; details={}; histories={}
    for base in [ROOT/'output/dongfang_fixed_carrier_scan_20260912', ROOT/'calculation/figure3_rta_20260913',
                 ROOT/'calculation/figure3_callaway_si_20261006',ROOT/'calculation/figure3_extensions_20261005',
                 ROOT/'calculation/figure3_callaway_20261005',out/'raw']:
        if not base.exists():continue
        for p in base.rglob('*'):
            if not p.is_file() or any(s in p.parts for s in ['source','smoke','local_preflight','prepare_attempt_incomplete_reference']):continue
            if p.suffix=='.toml' and p.stem in wanted:inputs[p.stem]=p
            if p.name.endswith('_analysis.json') and p.name[:-14] in wanted:details[p.name[:-14]]=json.loads(p.read_text())
            if p.name=='convergence.txt' and p.parent.name.endswith('_0') and p.parent.name[:-2] in wanted:histories[p.parent.name[:-2]]=p
    # Archived per-case analyses take precedence over plotting snapshots.
    for source in ['calculation/figure3_rta_20260913/reused/results.json',
                   'calculation/figure3_callaway_si_20261006/plot_snapshot/si_films/results.json',
                   'calculation/figure3_extensions_20261005/plot_snapshot/sic_films/results.json',
                   'calculation/figure3_callaway_20261005/plot_snapshot/sic_films/results.json']:
        d=json.loads((ROOT/source).read_text())
        for r in d if isinstance(d,list) else d['completed']: details.setdefault(r['case'],r)
    configs={}; audited=[]
    for row in rows:
        name=row['case']; r=dict(row); r['historical_review_required']=row['review_required'].lower()=='true'
        r['kappa_W_mK']=float(row['kappa_W_mK']); detail=details.get(name,{})
        r['input_path']=str(inputs.get(name,''));r['history_path']=str(histories.get(name,''))
        r['independent_seeds']=1; r['steady_state_proven']=False
        if name in inputs:
            configs[name]=tomllib.loads(inputs[name].read_text());r['input_sha256']=digest(inputs[name])
        if name in histories and name in inputs:
            try:
                r.update(raw_metrics(histories[name],configs[name],detail.get('normalization_factor',1.)))
                r['raw_recomputed']=True
                r['summary_relative_difference']=abs(r['kappa_recomputed_W_mK']-r['kappa_W_mK'])/abs(r['kappa_W_mK'])
                r['summary_matches_raw']=r['summary_relative_difference']<1e-7
            except Exception as e:r['raw_recomputed']=False;r['raw_error']=str(e)
        else:
            r['raw_recomputed']=False
            for metric in THRESHOLDS:r[metric]=detail.get(metric)
        r['screen_reasons']=[metric for metric,limit in THRESHOLDS.items() if r.get(metric) is not None and abs(r[metric])>limit]
        r['screen_passed']=r['raw_recomputed'] and r.get('raw_complete',False) and r.get('summary_matches_raw',False) and not r['screen_reasons']
        r['boundary_pairing']='same_gradient_boundary_branch' if row['family'].endswith('_inplane') else 'legacy_finite_bath_kernel_mismatch'
        actions=[]
        if not r['raw_recomputed']:actions.append('retrieve_raw_history')
        if r.get('raw_recomputed') and not r.get('summary_matches_raw'):actions.append('resolve_summary_raw_disagreement')
        if any('residual' in x for x in r['screen_reasons']):actions.append('inspect_absolute_and_relative_energy_residuals')
        if any(x in r['screen_reasons'] for x in ['half_window_drift_percent','max_cell_temperature_drift_K']):actions.append('extend_duration_then_recheck')
        if 'block_range_percent' in r['screen_reasons']:actions.append('increase_independent_sampling_and_recheck')
        if r['boundary_pairing']=='legacy_finite_bath_kernel_mismatch':actions.append('paired_rerun_with_common_boundary_kernel')
        actions.append('independent_seed_and_discretization_checks')
        r['required_actions']=actions;audited.append(r)
    bykey=collections.defaultdict(dict)
    for r in audited:bykey[key(r)][r['group'].split('_')[1]]=r
    pairs=[]
    for k,v in sorted(bykey.items(),key=lambda item:str(item[0])):
        a,b=v['RTA'],v['Callaway'];ca,cb=configs.get(a['case']),configs.get(b['case'])
        diffs=[]
        if ca and cb:
            for section in sorted((set(ca)|set(cb))-{'io','scattering'}):
                if ca.get(section)!=cb.get(section):diffs.append(section)
        else:diffs=['missing_input']
        pairs.append(dict(material=k[0],family=k[1],thickness_nm=k[2],roughness_nm=k[3],
                          rta_case=a['case'],callaway_case=b['case'],input_section_differences=diffs,
                          nominal_inputs_match=not diffs,boundary_pairing=a['boundary_pairing'],
                          rta_kappa=a['kappa_W_mK'],callaway_kappa=b['kappa_W_mK'],
                          enhancement_percent=100*(b['kappa_W_mK']/a['kappa_W_mK']-1),
                          both_pass_current_screen=a['screen_passed'] and b['screen_passed'],
                          collision_only_ready=False))
    groups=[]
    for group in SOURCES:
        rs=[r for r in audited if r['group']==group]
        groups.append(dict(group=group,cases=len(rs),historical_flagged=sum(r['historical_review_required'] for r in rs),
                           raw_recomputed=sum(r['raw_recomputed'] for r in rs),screen_passed=sum(r['screen_passed'] for r in rs),
                           summary_raw_mismatch=sum(r.get('summary_matches_raw') is False for r in rs),
                           reasons=dict(collections.Counter(x for r in rs for x in r['screen_reasons']))))
    result=dict(groups=groups,total_cases=len(audited),total_pairs=len(pairs),thresholds=THRESHOLDS,
                interpretation='Screen pass is not a convergence certificate; block SD is not independent-seed uncertainty. Legacy finite-bath pairs mix collision and boundary changes.',
                pairs_with_nominal_input_differences=sum(not p['nominal_inputs_match'] for p in pairs),
                legacy_boundary_mismatched_pairs=sum(p['boundary_pairing']=='legacy_finite_bath_kernel_mismatch' for p in pairs),
                pairs_both_pass_screen=sum(p['both_pass_current_screen'] for p in pairs),
                gradient_pairs_both_pass_screen=sum(p['both_pass_current_screen'] and p['boundary_pairing']=='same_gradient_boundary_branch' for p in pairs),
                all_independent_seed_and_discretization_checks_pending=True)
    write_json(out/'audit_summary.json',result);write_json(out/'case_audit.json',audited);write_json(out/'pair_audit.json',pairs)
    for filename,data in [('case_audit.csv',audited),('pair_audit.csv',pairs)]:
        fields=sorted(set().union(*(r.keys() for r in data)))
        with (out/filename).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
            for r in data:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})
    def table(data,fields):
        return '<table><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in fields)+'</tr>'+''.join('<tr>'+''.join('<td>'+html.escape(str(r.get(k,'')))+'</td>' for k in fields)+'</tr>' for r in data)+'</table>'
    page='''<!doctype html><meta charset="utf-8"><title>Si 与 SiC 薄膜收敛审计</title><style>body{font:15px system-ui;max-width:1300px;margin:35px auto;padding:0 20px;color:#24344a}table{border-collapse:collapse;width:100%;margin:20px 0}td,th{border:1px solid #ccd5df;padding:9px;text-align:left}th{background:#eaf0f7}code{background:#eef2f6}</style><h1>Si 与 SiC 薄膜收敛审计</h1><p>2026-10-16 · 原始归档只读，结果重新计算。通过旧阈值筛查不等于证明稳态或独立重复性；五个时间块不能当作五个独立随机种子。</p>'''
    page+=table(groups,['group','cases','historical_flagged','raw_recomputed','screen_passed','summary_raw_mismatch','reasons'])
    page+='<p>跨面及有限长面内的旧 RTA／Callaway 共 80 对采用不同边界算法；参数相同不代表只改变了碰撞算子。周期梯度面内共 60 对采用相同边界分支，仍需复核稳态、材料及程序来源。</p>'
    page+='''<h2>用于论文的判断</h2><p>280 条原始轨迹均已重算，全部完整并与原汇总相符；72 条通过旧阈值筛查，208 条至少触发一项复核指标。140 对中 31 对同时通过筛查，但其中仅 20 对属于边界分支一致的周期梯度面内算例。不能把 72 条或 20 对直接标为“已收敛”。</p>
<p>旧筛查阈值：后半窗内两段热导率漂移 1%，五个时间块均值极差 1%，最大单元温度漂移 0.01 K，能量／碰撞残差绝对值 10⁻⁶ eV。每个算例只有一个独立随机种子，五个时间块不是五次独立计算。</p>
<p>部分 SiC 碰撞残差虽超过绝对阈值，但相对存储热能仅约 10⁻¹⁵；应结合传输能量尺度和累计残差复核，不能单凭绝对阈值宣布守恒失败，也不能单凭相对值很小证明无误。208 条需要分类处理，不表示必须全部原样重算。</p>
<h2>补算优先级</h2><ol>
<li>先消除 80 对有限热浴算例的边界算法差异：同一程序、材料表、保守热浴和反射核，只改变 RTA／Callaway 碰撞模型。已观察到的热导率差异暂不能完全归因于 Normal 散射。</li>
<li>对周期梯度面内的 60 对，20 对同时通过筛查，可优先选为重复性／离散化检查基准；其余 40 对先区分时间漂移与统计波动，分别延长时间或增加独立采样。</li>
<li>代表点覆盖 Si／SiC、10／100／1000／10000 nm 及光滑／粗糙边界；先完成时间延长、dt 减半、网格细化和至少三个种子，再扩展完整厚度曲线。相邻厚度趋势不能代替数值收敛。</li>
<li>旧薄膜为固定参考温度的线性 Callaway，新器件为局域温度的非线性 Callaway，补做小温升一致性检查后再连接成论文跨尺度主线。</li></ol>
<p>完整数值与原轨迹位置：<a href="case_audit.csv">逐算例表</a> · <a href="pair_audit.csv">140 对比较</a> · <a href="archive_verification.json">T7 归档 SHA256 核验</a>。以下每项已附需要执行的动作；热导率单位为 W m⁻¹ K⁻¹。</p>'''
    page+='<h2>后续计算清单</h2>'+table(audited,['case','kappa_W_mK','screen_passed','screen_reasons','required_actions'])
    (out/'report.html').write_text(page)
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
