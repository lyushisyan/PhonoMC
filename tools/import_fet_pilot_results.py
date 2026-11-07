#!/usr/bin/env python3
"""Verify a completed-pilot export and independently recompute source temperature."""
import argparse,hashlib,json,tarfile
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('archive',type=Path);p.add_argument('--sha256',required=True)
p.add_argument('--root',type=Path,required=True)
a=p.parse_args();root=a.root.resolve();archive=a.archive.resolve()
assert hashlib.sha256(archive.read_bytes()).hexdigest()==a.sha256,'Archive checksum differs'
case=archive.name.removesuffix('.tar.gz')
with tarfile.open(archive) as t:
    for m in t:
        assert m.isfile() and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts
        data=t.extractfile(m).read();path=root/m.name
        assert path.resolve().is_relative_to(root),'Export path escapes result directory'
        if path.exists():assert path.read_bytes()==data,'Existing result differs: '+m.name
        else:path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
manifest=json.loads((root/'export_manifests'/(case+'.json')).read_text())
assert manifest['case']==case
for name,h in manifest['files'].items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,name
folder=root/'results/pilot'/(case+'_0')
with (folder/'convergence.txt').open() as f:
    header=f.readline().lstrip('#').split();data=np.loadtxt(f,ndmin=2)
grid=np.genfromtxt(folder/'grid_centers.csv',names=True,delimiter=',')
temp=data[:,[i for i,n in enumerate(header) if n.startswith('T_')]]
assert np.array_equal(grid['index'],np.arange(temp.shape[1]))
x,y,z=grid['x_nm'],grid['y_nm'],grid['z_nm']
if case.startswith('finfet_'):mask=(y>9)&(y<15)&(z>20)&(z<40)
else:mask=(y>7)&(y<17)&(((z>24.5)&(z<28.5))|((z>33)&(z<37))|((z>41.5)&(z<45.5)))
mask&=(x>36)&(x<42)
assert mask.sum()==1440
weights=np.exp(-.5*((x[mask]-39)/1.5)**2)*grid['volume_nm3'][mask]
local_last=float(np.average(temp[-1,mask],weights=weights)-300)
remote=json.loads((root/'analysis/pilot'/(case+'.json')).read_text())
assert abs(local_last-remote['metrics']['source_weighted_rise_K']['last_K'])<1e-9
record={'case':case,'archive_sha256':a.sha256,'verified_files':len(manifest['files']),
        'raw_last_step':int(data[-1,0]),'source_weighted_rise_K_recomputed':local_last,'matches_remote_analysis':True}
(root/'validation'/(case+'_local_verification.json')).write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
