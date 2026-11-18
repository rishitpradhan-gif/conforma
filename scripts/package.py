from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import hashlib
root=Path(__file__).resolve().parents[1]
target=root/'public/conforma-source.zip'
files=[]
for name in ['.gitignore','.gitattributes','.vercelignore','package.json','package-lock.json','README.md','LICENSE','server.mjs','cli.mjs','vercel.json']:
    files.append(root/name)
for folder in ['src','scripts','tests','docs','public']:
    for p in sorted((root/folder).rglob('*')):
        rel=p.relative_to(root)
        if p.is_file() and not p.is_symlink() and not any(x.startswith('.') or x in {'node_modules','__pycache__'} for x in rel.parts) and p.suffix!='.pyc' and p.name not in {'conforma-source.zip','conforma-source.sha256'}:
            files.append(p)
with ZipFile(target,'w',ZIP_DEFLATED) as z:
    for p in files:z.write(p,'conforma/'+str(p.relative_to(root)))
(root/'public/conforma-source.sha256').write_text(hashlib.sha256(target.read_bytes()).hexdigest()+'  conforma-source.zip\n')
print(f'Packaged {len(files)} files, {target.stat().st_size:,} bytes; no dependencies, credentials or environment files.')
