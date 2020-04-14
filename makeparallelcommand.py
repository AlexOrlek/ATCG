import sys,os,re
sourcedir=sys.argv[1]
sourcedir=str(sourcedir).rstrip('/')
filepathinfo=sys.stdin

for filepaths in filepathinfo:
    filepaths=filepaths.strip().split('\t')
    sample=filepaths[0]
    fastafilepath=filepaths[1]
    blastdbfilepath=filepaths[2]

    command='python %s/editfastaheaders.py %s %s | makeblastdb -dbtype nucl -out %s -title %s -logfile /dev/null'%(sourcedir,'arg',fastafilepath,blastdbfilepath,sample)
    print(command)
    



