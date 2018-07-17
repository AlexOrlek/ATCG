from mymod import runblastn,runsubprocess
import sys

#argv[1] is filepath to pipeline output folder for a given dataset
#argv[2] is query plasmids
#argv[3] is filepaths to plasmid databases
#argv[4] is evalue #1e-8 default
#argv[5] is wordsize #38 default
#argv[6] is threads


#create output directories
args=['mkdir -p %s/blast'%sys.argv[1]]
runsubprocess(args,shell=True)
args=['mkdir -p %s/output'%sys.argv[1]]
runsubprocess(args,shell=True)


#blast query plasmids against subject plasmid database
query='%s'%sys.argv[2] #plasmid query fasta file (all samples)
databasefiles='%s'%sys.argv[3]  #plasmid blast database; format: sample \t filepath to database
evalue=sys.argv[4]
wordsize=sys.argv[5]
threads=sys.argv[6]


with open(databasefiles) as f:
    for line in f:
        data=line.strip().split('\t')
        sample=data[0]
        database=data[1]
        blastoutput='%s/blast/%s_plasmidalignments.tsv'%(sys.argv[1],sample)
        runblastn(query, database, blastoutput, num_threads=threads, outfmt='custom qcov', evalue=evalue, word_size=wordsize, task='blastn',culling_limit='5')
print('runblastn finshed')


with open('%s/blast/blastsettings.txt'%sys.argv[1],'w') as f:
    f.write('e-value cutoff: %s\n'%evalue)
    f.write('word size: %s\n' %wordsize)
