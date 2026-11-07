#!/usr/bin/env python3
"""Read archived GNR traces without extracting or changing the T7 archives."""
import hashlib,json,tarfile,time
from pathlib import Path
import numpy as np

archive_root=Path('/Volumes/T7 Shield/PhonoMC_runs_20261013')
out=Path(__file__).resolve().parents[1]/'analysis/gnr_audit_20261017/t7'
out.mkdir(parents=True,exist_ok=True)
expected={line.split()[1]:line.split()[0] for line in (archive_root/'SHA256SUMS.txt').read_text().splitlines() if line.strip()}
class HashReader:
    def __init__(self,f):self.f=f;self.hash=hashlib.sha256()
    def read(self,n=-1):
        data=self.f.read(n);self.hash.update(data);return data

for campaign in ['graphene_LW_300K_20260911','gnr_width_2K_20260913','gnr_rectification_20260913','gnr_20K_W1000_20260914']:
    dest=out/(campaign+'.json')
    if dest.exists():print('Already audited',campaign,flush=True);continue
    path=archive_root/(campaign+'.tar.gz');start=time.monotonic()
    traces=[];metadata={};summaries={};members=0
    print('START',campaign,flush=True)
    with path.open('rb') as f:
        reader=HashReader(f)
        with tarfile.open(fileobj=reader,mode='r|gz') as t:
            for m in t:
                if not m.isfile():continue
                members+=1
                if m.name.endswith('convergence.txt'):
                    stream=t.extractfile(m);h=hashlib.sha256();header=stream.readline();h.update(header)
                    names=header.decode().lstrip('#').split()
                    wanted=[k for k in names if k in ['timestep','time_ps','heatflux','kappa_int','kappa_eff','particle_count','energy_balance_residual_ev','lifetime_energy_residual_ev']]
                    indices=[names.index(k) for k in wanted];rows=[];bad=0
                    for line in stream:
                        h.update(line)
                        if not line.strip():continue
                        parts=line.split()
                        if len(parts)!=len(names):bad+=1;continue
                        try:rows.append([float(parts[i]) for i in indices])
                        except ValueError:bad+=1
                    a=np.asarray(rows);record=dict(path=m.name,sha256=h.hexdigest(),bytes=m.size,columns=wanted,malformed_rows=bad,rows=len(rows))
                    if len(rows):
                        d={k:a[:,i] for i,k in enumerate(wanted)}
                        record.update(last_step=float(d['timestep'][-1]),end_ps=float(d['time_ps'][-1]),finite_selected_columns=bool(np.isfinite(a).all()))
                        end=d['time_ps'][-1];late=d['time_ps']>end/2;first=late&(d['time_ps']<=end*.75);last=d['time_ps']>end*.75
                        stats={}
                        for k in ['heatflux','kappa_eff','kappa_int','particle_count']:
                            if k not in d:continue
                            v=d[k];mean=float(v[late].mean())
                            drift=float(v[last].mean()-v[first].mean()) if first.any() and last.any() else None
                            stats[k]=dict(mean_late_half=mean,last=float(v[-1]),min=float(v.min()),max=float(v.max()),half_window_drift_percent=100*drift/abs(mean) if drift is not None and mean!=0 else None)
                        for k in ['energy_balance_residual_ev','lifetime_energy_residual_ev']:
                            if k in d:stats[k]={'max_abs':float(np.max(abs(d[k])))}
                        record['statistics']=stats
                    traces.append(record)
                elif (m.name.endswith(('.json','.toml')) and m.size<2_000_000) or m.name.endswith('summary.txt'):
                    raw=t.extractfile(m).read()
                    if m.name.endswith('.json'):
                        try:metadata[m.name]=json.loads(raw)
                        except (ValueError,UnicodeDecodeError):metadata[m.name]={'parse_error':True}
                    elif m.name.endswith('summary.txt'):summaries[m.name]=raw.decode(errors='replace')
                    else:metadata[m.name]=raw.decode(errors='replace')
        while reader.read(1024*1024):pass
        digest=reader.hash.hexdigest()
    assert digest==expected[path.name],('Archive checksum mismatch',path)
    result=dict(campaign=campaign,archive=str(path),archive_sha256=digest,archive_sha256_matches_migration_manifest=True,members=members,
                raw_traces=traces,metadata=metadata,summaries=summaries,elapsed_seconds=time.monotonic()-start,
                scope='Recomputed raw time/heat-flux/conductivity and available accounting columns; this is not a mesh/seed convergence certificate.')
    dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('DONE',campaign,'traces',len(traces),'seconds',round(result['elapsed_seconds'],1),flush=True)
