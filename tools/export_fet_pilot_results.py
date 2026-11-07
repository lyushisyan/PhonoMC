#!/usr/bin/env python3
"""Archive only completed, analyzed FET pilots with per-file checksums."""
import hashlib, json, sys, tarfile
from pathlib import Path

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

root=Path(sys.argv[1]).resolve()
exports=root/'exports';exports.mkdir(exist_ok=True)
manifests=root/'export_manifests';manifests.mkdir(exist_ok=True)
records=[]
for c in json.loads((root/'cases.json').read_text()):
    if c['stage']!='pilot':continue
    case=c['id'];receipt=root/'validation/pilot'/(case+'.json')
    analysis=root/'analysis/pilot'/(case+'.json')
    if not receipt.exists() or not analysis.exists():continue
    if json.loads(receipt.read_text())['state']!='passed':continue
    files=[root/'cases.json',root/'files_sha256.json',root/'validation/server_submission.json',root/'validation/smoke_gate.json',root/c['input']]
    files+=list((root/'results/pilot'/(case+'_0')).rglob('*'))
    files+=list((root/'analysis/pilot').glob(case+'*'))
    files+=list((root/'validation/resources').glob(case+'.*'))
    files+=[receipt,root/'logs'/('pilot_'+case+'.log')]
    files=[p for p in files if p.is_file()]
    manifest=manifests/(case+'.json')
    contents={'case':case,'files':{str(p.relative_to(root)):sha(p) for p in sorted(set(files))}}
    archive=exports/(case+'.tar.gz')
    if archive.exists():
        assert json.loads(manifest.read_text())==contents,'Completed results changed after export'
    else:
        manifest.write_text(json.dumps(contents,indent=2)+'\n')
        temporary=archive.with_suffix(archive.suffix+'.partial')
        assert not temporary.exists(),'Previous export is incomplete; preserve and inspect it first'
        with tarfile.open(temporary,'w:gz',compresslevel=1) as t:
            for p in files+[manifest]:t.add(p,arcname=str(p.relative_to(root)),recursive=False)
        temporary.rename(archive)
    records.append({'case':case,'archive':str(archive),'sha256':sha(archive),'bytes':archive.stat().st_size})
print(json.dumps(records,indent=2))
