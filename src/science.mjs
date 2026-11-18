// Distances are angstroms. Rigid fitting uses Horn's quaternion solution.
export const VERSION='1.0.0';
export const DOMAIN_RANGES={NMP:[[30,67]],LID:[[118,160]],CORE:[[1,29],[68,117],[161,214]]};
export const domainFor=id=>id>=30&&id<=67?'NMP':id>=118&&id<=160?'LID':'CORE';
const finitePoint=p=>Array.isArray(p)&&p.length===3&&[0,1,2].every(i=>Number.isFinite(p[i])&&Math.abs(p[i])<=1e6);
export function distance(a,b){if(!finitePoint(a)||!finitePoint(b))throw Error('Expected finite bounded 3D coordinates.');return Math.hypot(...a.map((v,i)=>v-b[i]));}
export function centroid(points){if(!Array.isArray(points)||!points.length||Array.from(points).some(p=>!finitePoint(p)))throw Error('Expected finite 3D coordinates.');return [0,1,2].map(i=>points.reduce((s,p)=>s+p[i],0)/points.length);}
export function parsePDB(text,{chain='A'}={}){
 if(typeof text!=='string'||text.length>10_000_000)throw Error('PDB input must be text under 10 MB.');
 const atoms=new Map();let firstModel=true,modelSeen=false;
 for(const line of text.split(/\r?\n/)){
  if(line.startsWith('MODEL ')){if(modelSeen)break;modelSeen=true;firstModel=true;}if(line.startsWith('ENDMDL')&&firstModel)break;
  if(!line.startsWith('ATOM  ')||line[21]!==chain)continue;
  const alt=line[16];if(alt!==' '&&alt!=='A')continue;
  const id=Number(line.slice(22,26)),name=line.slice(12,16).trim(),residue=line.slice(17,20).trim(),insertion=line[26]||' ';
  const values=[line.slice(30,38),line.slice(38,46),line.slice(46,54)];
  if(!line.slice(22,26).trim()||!Number.isSafeInteger(id)||!name||!residue||values.some(v=>!v.trim())||!finitePoint(values.map(Number)))throw Error('Malformed PDB atom coordinates.');
  const element=line.slice(76,78).trim()||name.replace(/[0-9]/g,'')[0];if(['H','D'].includes(element))continue;
  const key=[id,insertion,name].join(':');const atom={key,id,insertion,name,residue,element,position:values.map(Number),domain:domainFor(id)};
  if(!atoms.has(key)||alt===' ')atoms.set(key,atom);
 }
 if(!atoms.size)throw Error('No supported heavy atoms found in the requested chain.');return [...atoms.values()];
}
function largestEigenvector(matrix){
 const a=matrix.map(r=>r.slice()),v=Array.from({length:4},(_,i)=>Array.from({length:4},(_,j)=>+(i===j)));
 for(let iteration=0;iteration<100;iteration++){
  let p=0,q=1,max=0;for(let i=0;i<4;i++)for(let j=i+1;j<4;j++)if(Math.abs(a[i][j])>max){max=Math.abs(a[i][j]);p=i;q=j;}
  if(max<1e-12)break;const angle=.5*Math.atan2(2*a[p][q],a[q][q]-a[p][p]),c=Math.cos(angle),s=Math.sin(angle);
  const app=a[p][p],aqq=a[q][q],apq=a[p][q];for(let k=0;k<4;k++)if(k!==p&&k!==q){const kp=a[k][p],kq=a[k][q];a[k][p]=a[p][k]=c*kp-s*kq;a[k][q]=a[q][k]=s*kp+c*kq;}
  a[p][p]=c*c*app-2*s*c*apq+s*s*aqq;a[q][q]=s*s*app+2*s*c*apq+c*c*aqq;a[p][q]=a[q][p]=0;
  for(let k=0;k<4;k++){const kp=v[k][p],kq=v[k][q];v[k][p]=c*kp-s*kq;v[k][q]=s*kp+c*kq;}
 }
 let best=0;for(let i=1;i<4;i++)if(a[i][i]>a[best][best])best=i;return v.map(row=>row[best]);
}
export function rigidFit(moving,fixed){
 if(moving.length!==fixed.length||moving.length<3)throw Error('At least three paired points are required.');const cm=centroid(moving),cf=centroid(fixed),S=Array.from({length:3},()=>[0,0,0]);
 moving.forEach((p,k)=>{for(let i=0;i<3;i++)for(let j=0;j<3;j++)S[i][j]+=(p[i]-cm[i])*(fixed[k][j]-cf[j]);});
 const [[xx,xy,xz],[yx,yy,yz],[zx,zy,zz]]=S;
 const [w,x,y,z]=largestEigenvector([[xx+yy+zz,yz-zy,zx-xz,xy-yx],[yz-zy,xx-yy-zz,xy+yx,zx+xz],[zx-xz,xy+yx,-xx+yy-zz,yz+zy],[xy-yx,zx+xz,yz+zy,-xx-yy+zz]]);
 const rotation=[[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],[2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],[2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]];
 const rotate=p=>rotation.map(r=>r.reduce((sum,v,i)=>sum+v*p[i],0));const rc=rotate(cm),translation=cf.map((v,i)=>v-rc[i]);const apply=p=>{if(!finitePoint(p))throw Error('Invalid coordinate.');return rotate(p).map((v,i)=>v+translation[i]);};
 const rmsd=Math.sqrt(moving.reduce((sum,p,i)=>sum+distance(apply(p),fixed[i])**2,0)/moving.length);return {rotation,translation,rmsd,apply};
}
export function compareStructures(open,closed,{fit='core'}={}){
 if(!['core','all'].includes(fit))throw Error('Unknown alignment selection.');
 const b=new Map(closed.map(a=>[a.key,a]));const pairs=open.flatMap(a=>{const match=b.get(a.key);return match&&match.residue===a.residue?[{...a,open:a.position,closed:match.position}]:[];});
 const ca=pairs.filter(a=>a.name==='CA');if(ca.length<3)throw Error('Insufficient matched alpha carbons.');const chosen=fit==='core'?ca.filter(a=>a.domain==='CORE'):ca;const transform=rigidFit(chosen.map(a=>a.closed),chosen.map(a=>a.open));
 const atoms=pairs.map(a=>({...a,closed:transform.apply(a.closed)}));const residues=atoms.filter(a=>a.name==='CA').map(a=>({id:a.id,insertion:a.insertion,key:a.key,name:a.residue,domain:a.domain,open:a.open,closed:a.closed,displacement:distance(a.open,a.closed)}));
 const domains={};for(const name of ['CORE','NMP','LID']){const rows=residues.filter(r=>r.domain===name);domains[name]={residues:rows.length,rmsd:rows.length?Math.sqrt(rows.reduce((s,r)=>s+r.displacement**2,0)/rows.length):null,centroid_open:rows.length?centroid(rows.map(r=>r.open)):null,centroid_closed:rows.length?centroid(rows.map(r=>r.closed)):null};}
 return {schema:'conforma-comparison/1',engine:VERSION,fit,units:'angstrom',matched_atoms:atoms.length,matched_residues:residues.length,unmatched_open:open.length-atoms.length,unmatched_closed:closed.length-atoms.length,fit_residues:chosen.length,fit_rmsd:transform.rmsd,all_ca_rmsd:Math.sqrt(residues.reduce((s,r)=>s+r.displacement**2,0)/residues.length),rotation:transform.rotation,translation:transform.translation,domains,residues,atoms,limitation:'Matched heavy atoms, chain A, first model, blank/A alternate locations. C-alpha correspondence by residue number, insertion code, atom name and residue identity; not a general sequence alignment. Displacements depend on the chosen rigid fit. Interpolated coordinates are illustrative geometry, not a measured trajectory, molecular dynamics, kinetics or mechanistic proof.'};
}
export function interpolate(a,b,t){if(!finitePoint(a)||!finitePoint(b)||!Number.isFinite(t)||t<0||t>1)throw Error('Invalid interpolation.');return a.map((v,i)=>v+(b[i]-v)*t);}
export function measurePair(report,first,second){const aa=report.residues.filter(r=>r.id===first),bb=report.residues.filter(r=>r.id===second);if(aa.length!==1||bb.length!==1)throw Error('Choose two uniquely numbered modeled residues.');const [a]=aa,[b]=bb;return {first,second,open:distance(a.open,b.open),closed:distance(a.closed,b.closed),units:'angstrom',definition:'C-alpha to C-alpha distance in each observed endpoint.'};}
