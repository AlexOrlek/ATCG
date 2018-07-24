import sys,os,re


#arg[1] is output path; arg[2] is fasta directory; arg[3] is fasta file name; arg[4] is fastafile paths file

#renmame fastas (remove filename prefix) and create file with fasta filepaths

prefix=os.path.basename(sys.argv[3])
prefix=re.sub('.fasta$','',prefix) #remove .fasta suffix

f2=open(sys.argv[4],'w')
for filename in os.listdir(sys.argv[2]):
    if filename.endswith('.fasta'):
        newfilename=re.sub('^%s.id_'%prefix,'',filename) #remove prefix.id prefix
        if filename!=newfilename:
            os.rename('%s/%s'%(sys.argv[2],filename),'%s/%s'%(sys.argv[2],newfilename))
        sample=re.sub('.fasta$','',newfilename)
        fastafilepath=os.path.abspath('%s/%s'%(sys.argv[2],newfilename))
        f2.write('%s\t%s\n'%(sample,fastafilepath))
f2.close()
