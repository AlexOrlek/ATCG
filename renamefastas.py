import sys,os,re

#arg[1] is output path; arg[2] is fasta file

#renmame fastas (remove filename prefix) and create file with fasta filepaths

prefix=os.path.basename(sys.argv[2])
prefix=re.sub('.fasta$','',prefix) #remove .fasta suffix

f2=open('%s/fastafilepaths.tsv'%sys.argv[1],'w')
for filename in os.listdir('%s/splitfastas/'%sys.argv[1]):
    if filename.endswith('.fasta'):
        newfilename=re.sub('^%s.id_'%prefix,'',filename) #remove prefix.id prefix
        if filename!=newfilename:
            os.rename('%s/splitfastas/%s'%(sys.argv[1],filename),'%s/splitfastas/%s'%(sys.argv[1],newfilename))
        sample=re.sub('.fasta$','',newfilename)
        fastafilepath=os.path.abspath('%s/splitfastas/%s'%(sys.argv[1],newfilename))
        f2.write('%s\t%s\n'%(sample,fastafilepath))
f2.close()
