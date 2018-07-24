import sys,os
sys.path.append('./')
from pythonmods import makeBLASTdb, runsubprocess


#sys.argv[1] is outputpath; sys.argv[2] is fastadir; sys.argv[3] is fastafilepaths file; default: outpath/fastafilepaths.tsv (created by renamefastas.py); arg[4] is blastdb filepaths file 


args=['mkdir -p %s'%sys.argv[2]] #in case previous scripts not run
runsubprocess(args,shell=True,verbose=False)

f2=open(sys.argv[4],'w')
with open(sys.argv[3]) as f:
    for line in f:
        sample,fastafilepath=line.strip().split('\t')
        blastdbpath=os.path.abspath('%s/%s_db'%(sys.argv[2],sample))
        makeBLASTdb(fastafilepath, blastdbpath,'nucl')
        f2.write('%s\t%s\n'%(sample, blastdbpath))
f2.close()
