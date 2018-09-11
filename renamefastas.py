import sys,os,re


#arg[1] is output path; arg[2] is fasta directory; arg[3] is fasta file name; arg[4] is fastafile paths file; args[5] is blastdb paths file

#renmame fastas (remove filename prefix) and create file with fasta filepaths

prefix=os.path.basename(sys.argv[3])
prefix=re.sub('.fasta$|.fa$|.fna$','',prefix) #remove .fasta suffix

f2=open(sys.argv[4],'w')
f3=open(sys.argv[5],'w')
for filename in os.listdir(sys.argv[2]):
    if filename.endswith('.fasta') or filename.endswith('.fa') or filename.endswith('.fna'):
        newfilename=re.sub('^%s.id_'%prefix,'',filename) #remove prefix.id prefix
        newfilename=re.sub(r'(\S+)( .*?)(.fasta$|.fa$|.fna$)', r'\1.fasta', newfilename) #remove description and ensure ending is .fasta
        newfilename=re.sub(r'(\S+)(.fasta$|.fa$|.fna$)', r'\1.fasta', newfilename) #if no description, ensure ending is .fasta   
        if filename!=newfilename:
            os.rename('%s/%s'%(sys.argv[2],filename),'%s/%s'%(sys.argv[2],newfilename))
        sample=re.sub('.fasta$','',newfilename)
        fastafilepath=os.path.abspath('%s/%s'%(sys.argv[2],newfilename))
        blastdbfilepath=os.path.splitext(fastafilepath)[0]
        blastdbfilepath='%s_db'%blastdbfilepath
        f2.write('%s\t%s\n'%(sample,fastafilepath))
        f3.write('%s\t%s\n'%(sample,blastdbfilepath))
f2.close()
f3.close()
