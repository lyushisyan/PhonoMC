#!/usr/bin/env python3
"""Preserve completed FET stage outputs in immutable per-case checksum archives."""
import hashlib, json, sys, tarfile
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

root=Path(sys.argv[1]).resolve();stage=sys.argv[2]
assert stage in ['time200','time400','dt200','particles200','grid200','grid_control200','seeds200']
exports=root/'exports'/stage;exports.mkdir(parents=True,exist_ok=True)
manifests=root/'export_manifests'/stage;manifests.mkdir(parents=True,exist_ok=True)
records=[]
for c in json.loads((root/'cases.json').read_text()):
    if c['stage']!=stage:continue
    case=c['id'];receipt=root/'validation'/stage/(case+'.json')
    analysis=root/'analysis'/stage/(case+'.json')
    if not receipt.exists() or not analysis.exists():continue
    if json.loads(receipt.read_text())['state']!='passed':continue
    files=[root/'cases.json',root/'files_sha256.json',root/'stage_readiness.json',
           root/'validation'/('submission_'+stage+'.json'),root/'validation/smoke_gate.json',root/c['input']]
    files+=list((root/'results'/stage/(case+'_0')).rglob('*'))
    files+=list((root/'analysis'/stage).glob(case+'*'))
    files+=list((root/'validation/resources'/stage).glob(case+'.*'))
    files+=[receipt,root/'logs'/(stage+'_'+case+'.log')]
    files=sorted(set(p for p in files if p.is_file()))
    manifest=manifests/(case+'.json')
    contents={'case':case,'stage':stage,'files':{str(p.relative_to(root)):sha(p) for p in files}}
    archive=exports/(case+'.tar.gz')
    if archive.exists():
        assert json.loads(manifest.read_text())==contents,'Completed results changed after export'
    else:
        assert not manifest.exists(),'Prior incomplete export; preserve and inspect'
        manifest.write_text(json.dumps(contents,indent=2)+'\n')
        temporary=archive.with_suffix(archive.suffix+'.partial')
        assert not temporary.exists()
        with tarfile.open(temporary,'w:gz',compresslevel=1) as t:
            for p in files+[manifest]:t.add(p,arcname=str(p.relative_to(root)),recursive=False)
        temporary.rename(archive)
    records.append({'case':case,'stage':stage,'archive':str(archive),'sha256':sha(archive),'bytes':archive.stat().st_size})
print(json.dumps(records,indent=2))
