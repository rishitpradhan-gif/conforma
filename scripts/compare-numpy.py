"""Independent optional check: NumPy SVD vs the JavaScript quaternion fit."""
import json,subprocess
from pathlib import Path
import numpy as np
import scipy
root=Path(__file__).resolve().parents[1]
# Export the JavaScript inputs and outputs rather than reusing its rotation in the oracle.
code="""import {readFile} from 'node:fs/promises';import {compareStructures,rigidFit} from './src/science.mjs';const d=JSON.parse(await readFile('public/data/structures.json'));const reports=['core','all'].map(fit=>compareStructures(d.open,d.closed,{fit}));let seed=713;const rnd=()=>{seed=(1664525*seed+1013904223)>>>0;return seed/4294967296;};const cases=[];for(let k=0;k<30;k++){const a=Array.from({length:20},()=>[rnd()*20,rnd()*20,rnd()*20]);const angle=rnd()*6.28,offset=[rnd()*100,rnd()*100,rnd()*100],b=a.map(([x,y,z])=>[x*Math.cos(angle)-y*Math.sin(angle)+offset[0],x*Math.sin(angle)+y*Math.cos(angle)+offset[1],z+offset[2]]);const f=rigidFit(a,b);cases.push({a,b,rmsd:f.rmsd,aligned:a.map(f.apply)});}console.log(JSON.stringify({reports,cases}));"""
run=subprocess.run(['node','--input-type=module','-e',code],cwd=root,capture_output=True,text=True,check=True)
data=json.loads(run.stdout)
def fit(moving,fixed):
    moving=np.array(moving,float);fixed=np.array(fixed,float)
    cm=moving.mean(axis=0);cf=fixed.mean(axis=0)
    u,s,vt=np.linalg.svd((moving-cm).T@(fixed-cf))
    d=np.eye(3);d[2,2]=np.linalg.det(u@vt)
    rotation=u@d@vt
    aligned=(moving-cm)@rotation+cf
    return rotation,cm,cf,aligned,float(np.sqrt(np.mean(np.sum((aligned-fixed)**2,axis=1))))
results=[]
source=json.loads((root/'public/data/structures.json').read_text())
open_ca={a['key']:a for a in source['open'] if a['name']=='CA'}
closed_ca={a['key']:a for a in source['closed'] if a['name']=='CA'}
for report in data['reports']:
    keys=[k for k,v in open_ca.items() if k in closed_ca and (report['fit']=='all' or v['domain']=='CORE')]
    rot,cm,cf,aligned,rmsd=fit([closed_ca[k]['position'] for k in keys],[open_ca[k]['position'] for k in keys])
    js={a['key']:a for a in report['residues']}
    error=max(float(np.max(np.abs((np.array(closed_ca[k]['position'])-cm)@rot+cf-np.array(js[k]['closed'])))) for k in open_ca if k in js)
    assert error<1e-8 and abs(rmsd-report['fit_rmsd'])<1e-9
    results.append(dict(fit=report['fit'],fit_residues=len(keys),fit_rmsd=rmsd,max_coordinate_difference_angstrom=error))
max_error=0
for case in data['cases']:
    _,_,_,aligned,rmsd=fit(case['a'],case['b'])
    error=float(np.max(np.abs(aligned-np.array(case['aligned']))));max_error=max(max_error,error)
    assert error<1e-8 and abs(rmsd-case['rmsd'])<1e-9
out=dict(schema='conforma-independent-check/1',numpy=np.__version__,scipy=scipy.__version__,method='NumPy Kabsch SVD with proper-rotation correction; independent of JS Horn quaternion eigensolver.',real_structure_checks=results,seeded_rigid_transform_checks=30,seed=713,max_synthetic_coordinate_error_angstrom=max_error,passed=True,limitation='Numerical agreement only. Does not validate physical intermediate states, biological mechanism, crystal structures, or visual accuracy.')
(root/'docs/numpy-comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
