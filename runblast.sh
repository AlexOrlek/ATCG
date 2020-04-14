#!/bin/bash
set -e
set -u
set -o pipefail

#argv[1] is filepath to pipeline output folder for a given dataset
#argv[2] is blastdb directory
#argv[3] is filepathinfo.tsv file with filepaths to subject sequence databases
#argv[4] is evalue #1e-8 default
#argv[5] is wordsize #38 default
#argv[6] is blast task
#argv[7] is threads
#argv[8] is blast bi-directional boolean flag (default False)
#argv[9] is blast type (-s flag or -1/-2 flags)

mkdir -p ${1}/blast
mkdir -p ${1}/output

blastdbdir=${2}
filepathinfo=${3} #format: sample \t filepath to fasta \t filepath to blast database
evalue=${4}
wordsize=${5}
task=${6}
threads=${7}
bidirectionalblast=${8} #default is False
blasttype=${9}  #allvallpairwise (-s flag) or pairwise (-1/-2 flags)

sort -k1,1V -o ${filepathinfo} ${filepathinfo}  #-k command must precede -o command

if [ ${blasttype} != 'pairwiserun2' ]; then #prevents writing to file twice in the case of pairwise blasttype
    > ${1}/allsubjects.txt  #used in reformatblastoutput.sh
fi

if [ ${blasttype} == 'allvallpairwise' ]; then
    if [ ${bidirectionalblast} == 'False' ]; then
        numsamples=$( cat ${filepathinfo} | wc -l )
        counter=0
        while IFS=$'\t' read -r -a line
        do
            let counter=$counter+1
            if [[ $counter -eq $numsamples ]]; then
                continue
            fi
            sample="${line[0]}"
            database="${line[2]}"
            echo "${sample}" >> ${1}/allsubjects.txt
            blastoutput="${1}/blast/${sample}_alignments.tsv"
            if [[ $counter -gt 1 ]]; then
                cat `find ${blastdbdir}/ -maxdepth 1 -mindepth 1 -name "*.fasta" ! -name "${sample}.fasta" | grep -v -F -f <(printf "%s\n" "${excludesamples[@]}")` | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit '5'
                excludesamples=(${excludesamples[@]} "${blastdbdir}/${sample}.fasta")
            else
                cat `find ${blastdbdir}/ -maxdepth 1 -mindepth 1 -name "*.fasta" ! -name "${sample}.fasta"` | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit '5'
                excludesamples=("${blastdbdir}/${sample}.fasta")
            fi
        done < ${filepathinfo}
    else
        while IFS=$'\t' read -r -a line
        do
            sample="${line[0]}"
            database="${line[2]}"
            echo "${sample}" >> ${1}/allsubjects.txt
            blastoutput="${1}/blast/${sample}_alignments.tsv"
            cat `find ${blastdbdir}/ -maxdepth 1 -mindepth 1 -name "*.fasta" ! -name "${sample}.fasta"` | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit '5'
        done < ${filepathinfo}
    fi
else   #blasttype==pairwise (or pairwiserun2 in the case of bidirectionalblast second run)
    while IFS=$'\t' read -r -a line
    do
        sample="${line[0]}"
        database="${line[2]}"
        echo "${sample}" >> ${1}/allsubjects.txt
        blastoutput="${1}/blast/${sample}_alignments.tsv"
        cat `find ${blastdbdir}/ -maxdepth 1 -mindepth 1 -name "*.fasta"` | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit '5'
    done < ${filepathinfo}
fi


if [ ${blasttype} != 'pairwiserun2' ]; then #prevents writing to file twice in the case of pairwise blasttype
    > ${1}/blastsettings.txt
    echo "e-value cutoff: ${evalue}" >> ${1}/blastsettings.txt
    echo "word size: ${wordsize}" >> ${1}/blastsettings.txt
    echo "bi-directional BLAST: ${bidirectionalblast}" >> ${1}/blastsettings.txt
fi


#OLD CODE

    #cat ${query} | seqkit grep -r -p ^${sample} -v  | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt ${outfmt} -task 'blastn' -num_threads ${threads} -max_target_seqs '500' -word_size ${wordsize} -culling_limit '5'
    #cat `find ${blastdbdir}/ -maxdepth 1 -mindepth 1 -name "*.fasta"` | seqkit grep -r -p ^${sample} -v  | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit '5'


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
# done < ${filepathinfo}

