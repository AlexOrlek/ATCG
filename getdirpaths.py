import sys,os,re
from pythonmods import runsubprocess

dirpath=sys.argv[1]  #args.sequences directory path
filepathinfo=sys.argv[2]
blastdbdir=sys.argv[3]  #actually where blastdbs are stored
blasttype=sys.argv[4]

runsubprocess(['mkdir -p %s'%blastdbdir],shell=True)

directory=str(dirpath).rstrip('/')
dircontents=os.listdir(directory)

samples=set()
f2=open(filepathinfo,'w')
for dircontent in dircontents:
    filepath='%s/%s'%(directory,dircontent)
    if os.path.isfile(filepath): #check for fasta files...
        if filepath.endswith('.gz'):
            gunzipfilepath=re.sub(r'\.gz$','',filepath)
            extension=os.path.splitext(gunzipfilepath)[1]
            sample=os.path.splitext(os.path.basename(gunzipfilepath))[0]
        else:
            extension=os.path.splitext(filepath)[1]
            sample=os.path.splitext(os.path.basename(filepath))[0]
        if extension in {'.fa','.fasta','.fna'}:
            if sample not in samples: #skip duplicates e.g. sample.fa and sample.fa.gz
                samples.add(sample)
                blastdbpath='%s/%s_db'%(blastdbdir,sample)
                f2.write('%s\t%s\t%s\n'%(sample,filepath,blastdbpath))
f2.close()

if blasttype=='allvallpairwise':
    assert len(samples)>1,'Error: at least 2 fasta files must be provided in the directory: %s'%directory
else:
    assert len(samples)>0,'Error: at least 1 fasta file must be provided in the directory'

