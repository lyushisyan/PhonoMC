#!/usr/bin/env python3
"""Produce a scientific transient/cost report from four verified server pilots."""
import html,json,os,re
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/phonomc-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root=Path(__file__).resolve().parents[1]
data=root/'calculation/fet_comparison_20261016/server_results'
out=root/'analysis/fet_pilot_review_20261017';out.mkdir(parents=True,exist_ok=True)
order=['finfet_rta','finfet_callaway','gaafet_rta','gaafet_callaway']
labels=['FinFET RTA','FinFET Callaway','GAAFET RTA','GAAFET Callaway']
colors=['#2468a2','#2468a2','#bc5627','#bc5627']
rows=[];traces=[]
for short,label in zip(order,labels):
    case=short+'_P1uW_s12345'
    verification=json.loads((data/'validation'/(case+'_local_verification.json')).read_text())
    assert verification['matches_remote_analysis']
    a=json.loads((data/'analysis/pilot'/(case+'.json')).read_text())
    r=json.loads((data/'validation/pilot'/(case+'.json')).read_text())
    assert r['state']=='passed' and a['complete']
    resource=(data/'validation/resources'/(case+'.time.txt')).read_text()
    cpu=float(re.search(r'Percent of CPU this job got:\s*([\d.]+)%',resource)[1])/100
    m=a['metrics'];summary=dict(case=case,label=label,end_ps=a['end_ps'],wall_minutes=a['wall_seconds_solver_and_validation']/60,
        peak_rss_GiB=a['max_rss_kib_runner_and_solver']/1024**2,average_active_CPU_cores=cpu,reserved_cores=28,
        source_rise_K=m['source_weighted_rise_K']['last_K'],source_drift_percent=m['source_weighted_rise_K']['late_half_drift_percent'],
        channel_rise_K=m['channel_mean_rise_K']['last_K'],channel_drift_percent=m['channel_mean_rise_K']['late_half_drift_percent'],
        peak_rise_K=m['peak_rise_K']['last_K'],layer_rises_K=[v['last_K'] for k,v in m.items() if k.startswith('layer_')],
        max_energy_residual_ev=r['metrics']['max_energy_residual_ev'],max_collision_residual_ev=r['metrics']['max_collision_residual_ev'],
        max_source_ledger_residual_ev=r['metrics']['max_source_ledger_residual_ev'],pending_energy_fraction=r['metrics']['pending_energy_fraction'],
        estimated_200ps_hours=a['wall_seconds_solver_and_validation']*4/3600,
        estimated_2ns_hours=a['linear_2ns_wall_estimate_seconds']/3600,output_MiB=a['output_bytes']/1024**2,
        archive_sha256=verification['archive_sha256'],steady_state_proven=False)
    rows.append(summary)
    traces.append(np.genfromtxt(data/'analysis/pilot'/(case+'_transient.csv'),delimiter=',',names=True))
(out/'summary.json').write_text(json.dumps({'cases':rows,'independent_seeds_per_case':1,'all_raw_files_verified':True,'all_transport_complete':True,'steady_state_proven':False},indent=2)+'\n')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(11,7.8),constrained_layout=True)
for i,(a,label,color) in enumerate(zip(traces,labels,colors)):
    for ax,key in [(axes[0,0],'source_weighted_rise_K'),(axes[0,1],'channel_mean_rise_K')]:
        ax.plot(a['time_ps'],a[key],color=color,ls='--' if i%2 else '-',marker='o',ms=3,label=label)
for ax,title in [(axes[0,0],'(a) Source-weighted temperature rise'),(axes[0,1],'(b) Channel-mean temperature rise')]:
    ax.set(title=title,xlabel='Time (ps)',ylabel='Temperature rise (K)');ax.grid(alpha=.18);ax.legend(fontsize=8)
for i in [2,3]:
    for layer,color in zip([1,2,3],['#2468a2','#bc5627','#4e9273']):
        a=traces[i];axes[1,0].plot(a['time_ps'],a[f'layer_{layer}_mean_rise_K'],color=color,ls='--' if i==3 else '-',label=f'Layer {layer} '+('Callaway' if i==3 else 'RTA'))
axes[1,0].set(title='(c) GAAFET layer means (1 = bottom)',xlabel='Time (ps)',ylabel='Temperature rise (K)');axes[1,0].legend(fontsize=8,ncol=2);axes[1,0].grid(alpha=.18)
bars=axes[1,1].bar(range(4),[r['wall_minutes'] for r in rows],color=['#2468a2','#74a1c5','#bc5627','#d69d7b'])
axes[1,1].bar_label(bars,fmt='%.1f min',padding=3,fontsize=9)
axes[1,1].set(title='(d) Measured runtime: 50 ps, 10M carriers',ylabel='Wall time (min)',xticks=range(4),xticklabels=['Fin\nRTA','Fin\nCallaway','GAA\nRTA','GAA\nCallaway'],ylim=(0,68))
fig.suptitle('Matched 1 µW devices: transient pilots, not steady-state results',fontsize=14)
fig.savefig(out/'pilot_transients.png',dpi=180);fig.savefig(out/'pilot_transients.pdf');plt.close(fig)
def table(fields):
    head='<tr>'+''.join('<th>'+html.escape(label)+'</th>' for _,label in fields)+'</tr>'
    body=''
    for r in rows:
        body+='<tr>'+''.join('<td>'+html.escape((f'{r[k]:.3e}' if 0<abs(r[k])<1e-4 else f'{r[k]:.3f}') if isinstance(r[k],float) else str(r[k]))+'</td>' for k,_ in fields)+'</tr>'
    return '<table>'+head+body+'</table>'
page='''<!doctype html><meta charset="utf-8"><title>FinFET / GAAFET 1 μW 试算核对</title>
<style>body{font:16px system-ui;line-height:1.7;max-width:1150px;margin:35px auto;padding:0 20px;color:#20334a}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border:1px solid #cad3df;padding:8px}th{background:#eaf1f7}img{width:100%}.note{padding:15px;background:#fff3d9}</style>
<h1>FinFET／三层 GAAFET：1 μW 试算核对</h1><p>2026-10-17（Europe/Moscow）· 数组 453039_0–3 全部 COMPLETED，退出码 0。原始轨迹、输入、日志、资源记录、提交回执和分析均已回传；四个归档及逐文件 SHA256 已核验，并在本地独立重算热源加权温升。</p>
<p class="note"><b>结论：计算可运行、能量检查通过，但四项均未证明达到稳态。</b>每项仅一个随机种子，后半程只有 5 个稀疏输出点。不能把单时刻峰值或末段两个均值接近当作稳态与统计收敛证明。</p>
<img src="pilot_transients.png" alt="四个试算的温升轨迹、GAAFET分层温升与运行时间">
<h2>50 ps 末时刻结果与后段漂移</h2><p>初始及热浴均为 300 K。热源加权温升按固定高斯体积热源权重计算；沟道及分层均温覆盖各自完整沟道。漂移比较 25–37.5 ps 与 37.5–50 ps 的存储样本均值，除以后半窗均值。离散最大温升对采样噪声敏感。</p>'''
page+=table([('label','算例'),('source_rise_K','热源温升 K'),('source_drift_percent','热源漂移 %'),('channel_rise_K','沟道温升 K'),('channel_drift_percent','沟道漂移 %'),('peak_rise_K','单格最大温升 K')])
page+='<p>FinFET 两模型末时刻热源温升相差约 0.033 K，GAAFET 相差约 0.032 K；各自时间漂移约 0.44–0.84 K，远大于模型间差值。现在不能据此宣称稳态 Callaway 降温效果。GAAFET 在该暂态时刻温升较高，也不能直接当作稳态热阻比。</p>'
page+='<h2>三层均温与守恒检查</h2>'+table([('label','算例'),('layer_rises_K','层平均温升 K（下→上）'),('max_energy_residual_ev','最大能量残差 eV'),('max_collision_residual_ev','最大碰撞残差 eV')])
page+='<p>上述残差的科学计数法精确值见 summary.json：最大能量残差小于 8.2×10⁻¹¹ eV，碰撞残差小于 4.7×10⁻¹⁰ eV，记录的热源账本残差和末时刻待沉积能量均为零；无粗糙边界回退采样。每项累计规定输入能量 312.0754537 eV，与 1 μW×50 ps 一致。这些为记录点与最终回执检查，不替代长时间累计守恒检查。</p>'
page+='<h2>实测运行成本与外推</h2>'+table([('label','算例'),('wall_minutes','实测分钟'),('peak_rss_GiB','峰值 GiB'),('average_active_CPU_cores','平均占用核数'),('estimated_200ps_hours','200 ps 外推小时'),('estimated_2ns_hours','2 ns 外推小时')])
page+='''<p>时间来自运行回执，内存及 CPU 来自 GNU time，Slurm 记录单独保留。每项申请 28 核；FinFET 在 c431、GAAFET 在 c311 运行，因此跨器件成本包含节点差异。28 核并不表示所有时间均使用满核。外推按时长线性缩放，包含初始化，未验证长时间非线性碰撞成本；不能视作承诺的完成时间。</p>
<h2>下一批收敛安排（输入已准备，尚未提交）</h2><ol>
<li><b>先跑四项 200 ps：</b>保持 1 μW、10M 载体、dt=0.05 ps 和原网格，28 核／8 GiB／每项 8 小时上限；预计顺序计算合计约 11.3 小时。每 5 ps 保存温度，不输出全场矢量热流。从 300 K 重新运行，不假定现有程序支持断点续算。</li>
<li>若 200 ps 后段仍漂移，再使用已准备的 400 ps 输入，必要时继续延长。先对固定区域均温、分层均温及热源温升做时间窗检查，再讨论峰值。</li>
<li>时间稳定后，更新所有对照的总时长，再进行 dt=0.025 ps、20M 载体检查。阈值暂以均值变化 &lt;1% 为起点；若模型差异小于 1%，需要更严格的误差预算。</li>
<li>网格三向各细化两倍至 0.5×0.5×0.25 nm³，使用 80M 载体；另准备同为 80M 的粗网格对照，避免同时改变网格与载体数而混淆判断。此组成本较高，后置安排，48 GiB 申请尚需实测确认。</li>
<li>补种子 12346、12347，与已接受时长的 12345 构成三次独立重复；误差条不能使用五个相邻时间块冒充独立样本。</li></ol>
<p>32 个收敛输入与 8 个备用短测输入位于 calculation/fet_convergence_20261017。只有 time200 阶段标为可运行，其余阶段需要根据时间延长结果审查。32 个输入均通过实际 C++ 几何、边界与分层功率检查；未提交原来的 36 项功率扫描。</p>
<p><a href="summary.json">完整数值 JSON</a> · <a href="pilot_transients.pdf">可导出 PDF 图</a></p>'''
(out/'report.html').write_text(page)
print(json.dumps({'report':str(out/'report.html'),'cases':len(rows),'first_200ps_total_estimate_hours':sum(r['estimated_200ps_hours'] for r in rows)},indent=2))
