import sys,os
sys.path.append('./')
from pythonmods import makeBLASTdb, runsubprocess

#sys.argv[1] is outputpath; sys.argv[2] is fastafilepaths file; default: outpath/fastafilepaths.tsv (created by renamefastas.py) 
args=['mkdir -p %s/splitfastas'%sys.argv[1]] #in case previous scripts not run
runsubprocess(args,shell=True,verbose=False)

f2=open('%s/blastdbfilepaths.tsv'%sys.argv[1],'w')
with open(sys.argv[2]) as f:
    for line in f:
        sample,fastafilepath=line.strip().split('\t')
        blastdbpath=os.path.abspath('%s/splitfastas/%s_db'%(sys.argv[1],sample))
        makeBLASTdb(fastafilepath, blastdbpath,'nucl')
        f2.write('%s\t%s\n'%(sample, blastdbpath))
f2.close()
