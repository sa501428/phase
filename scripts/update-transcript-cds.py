# Add genomic CDS bounds to the embedded RefSeq transcript tuples from the same
# public source downloads recorded in references/manifest.json.
import gzip,json,base64,hashlib
from pathlib import Path
for build in ['hg19','hg38']:
 source=Path('/tmp/phase-'+build+'-genes.gz')
 manifest=json.loads(Path('references/manifest.json').read_text())
 assert hashlib.sha256(source.read_bytes()).hexdigest()==manifest[build]['sha256']
 text=Path('references/'+build+'.js').read_text();packed=text.split("'",2)[1]
 genes=json.loads(gzip.decompress(base64.b64decode(packed)));lookup={}
 for line in gzip.open(source,'rt'):
  a=line.rstrip().split('\t');lookup[(a[1],a[2],int(a[4]),int(a[5]))]=[int(a[6]),int(a[7]),a[13],a[14]]
 for g in genes:
  for t in g[6]:t[4:]=lookup[(t[0],g[1],t[1],t[2])]
 packed=base64.b64encode(gzip.compress(json.dumps(genes,separators=(',',':')).encode(),mtime=0)).decode()
 Path('references/'+build+'.js').write_text("export default '"+packed+"';\n")
