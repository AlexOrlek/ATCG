#!/bin/bash
set -e
set -u
set -o pipefail

#argv[1] is filepath to pipeline output folder for a given dataset
#argv[2] is splitfasta directory
#argv[3] is filepaths to subject sequence databases
#argv[4] is evalue #1e-8 default
#argv[5] is wordsize #38 default
#argv[6] is threads

mkdir -p ${1}/blast
mkdir -p ${1}/output

splitfastadir=${2}
databasefiles=${3} #subject sequence blast database; format: sample \t filepath to database
evalue=${4}
wordsize=${5}
task=${6}
threads=${7}

sort -k1,1V -o ${databasefiles} ${databasefiles}  #-k command must precede -o command


while IFS=$'\t' read -r -a line
do
    sample="${line[0]}"
    database="${line[1]}"
    blastoutput="${1}/blast/${sample}_alignments.tsv"
    #cat ${query} | seqkit grep -r -p ^${sample} -v  | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt ${outfmt} -task 'blastn' -num_threads ${threads} -max_target_seqs '500' -word_size ${wordsize} -culling_limit '5'
    cat `find ${splitfastadir}/ -maxdepth 1 -mindepth 1 -name "*.fasta"` | seqkit grep -r -p ^${sample} -v  | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit '5'
done < ${databasefiles}

> ${1}/blast/blastsettings.txt
echo "e-value cutoff: ${evalue}" >> ${1}/blast/blastsettings.txt
echo "word size: ${wordsize}" >> ${1}/blast/blastsettings.txt




#OLD CODE

#outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen'
#sample='H112800346_9309cef4-089e-4b7f-938f-d6555e74e174'
#database='/home/alex/Documents/testing/ATCGpipeline/mankpc/samplelevel/subsetoutput/splitfastas/H112800346_9309cef4-089e-4b7f-938f-d6555e74e174_db'
#blastoutput="${1}/blast/${sample}_plasmidalignments.tsv"
#cat ${query} | seqkit grep -r -p ^${sample} -v  | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task 'blastn' -num_threads ${threads} -max_target_seqs '500' -word_size ${wordsize} -culling_limit '5'


# while IFS=$'\t' read -r -a line
# do
#     sample="${line[0]}"
#     database="${line[1]}"
#     blastoutput="${1}/blast/${sample}_plasmidalignments.tsv"
#     cat ${query} | seqkit grep -r -p ^${sample} -v | python runblast.py ${database} ${blastoutput} ${threads} ${evalue} ${wordsize}
# done < ${databasefiles}

