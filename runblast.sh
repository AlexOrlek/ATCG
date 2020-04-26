#!/bin/bash
set -e
set -u
set -o pipefail

#arg[1] is filepath to pipeline output folder for a given dataset
#arg[2] is blastdb directory
#arg[3] is filepathinfo.tsv file with filepaths to subject sequence databases
#arg[4] is evalue #1e-8 default
#arg[5] is wordsize #38 default
#arg[6] is blast task
#arg[7] is culling limit
#arg[8] is threads
#arg[9] is blast bi-directional boolean flag (default False)
#arg[10] is blast type (-s flag or -1/-2 flags)

outputpath=${1}
blastdbdir=${2}
filepathinfo=${3} #format: sample \t filepath to fasta \t filepath to blast database
evalue=${4}
wordsize=${5}
task=${6}
cullinglimit=${7}
threads=${8}
bidirectionalblast=${9} #default is False
blasttype=${10}  #allvallpairwise (-s flag) or pairwise (-1/-2 flags)

mkdir -p ${outputpath}/blast
mkdir -p ${outputpath}/output

sort -k1,1V -o ${filepathinfo} ${filepathinfo}  #-k command must precede -o command

if [ ${blasttype} != 'pairwiserun2' ]; then #prevents writing to file twice in the case of pairwise blasttype
    > ${outputpath}/allsubjects.txt  #used in reformatblastoutput.sh
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
            echo "${sample}" >> ${outputpath}/allsubjects.txt
            blastoutput="${outputpath}/blast/${sample}_alignments.tsv"
            if [[ $counter -gt 1 ]]; then
                cat `find ${blastdbdir}/ -maxdepth 1 -mindepth 1 -name "*.fasta" ! -name "${sample}.fasta" | grep -v -F -f <(printf "%s\n" "${excludesamples[@]}")` | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit ${cullinglimit}
                excludesamples=(${excludesamples[@]} "${blastdbdir}/${sample}.fasta")
            else
                cat `find ${blastdbdir}/ -maxdepth 1 -mindepth 1 -name "*.fasta" ! -name "${sample}.fasta"` | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit ${cullinglimit}
                excludesamples=("${blastdbdir}/${sample}.fasta")
            fi
        done < ${filepathinfo}
    else
        while IFS=$'\t' read -r -a line
        do
            sample="${line[0]}"
            database="${line[2]}"
            echo "${sample}" >> ${outputpath}/allsubjects.txt
            blastoutput="${outputpath}/blast/${sample}_alignments.tsv"
            cat `find ${blastdbdir}/ -maxdepth 1 -mindepth 1 -name "*.fasta" ! -name "${sample}.fasta"` | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit ${cullinglimit}
        done < ${filepathinfo}
    fi
else   #blasttype==pairwise (or pairwiserun2 in the case of bidirectionalblast second run)
    while IFS=$'\t' read -r -a line
    do
        sample="${line[0]}"
        database="${line[2]}"
        echo "${sample}" >> ${outputpath}/allsubjects.txt
        blastoutput="${outputpath}/blast/${sample}_alignments.tsv"
        cat `find ${blastdbdir}/ -maxdepth 1 -mindepth 1 -name "*.fasta"` | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit ${cullinglimit}
    done < ${filepathinfo}
fi


if [ ${blasttype} != 'pairwiserun2' ]; then #prevents writing to file twice in the case of pairwise blasttype
    > ${outputpath}/blastsettings.txt
    echo "BLAST task: ${task}" >> ${outputpath}/blastsettings.txt
    echo "e-value cutoff: ${evalue}" >> ${outputpath}/blastsettings.txt
    echo "word size: ${wordsize}" >> ${outputpath}/blastsettings.txt
    echo "bi-directional BLAST: ${bidirectionalblast}" >> ${outputpath}/blastsettings.txt
fi


#OLD CODE

    #cat ${query} | seqkit grep -r -p ^${sample} -v  | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt ${outfmt} -task 'blastn' -num_threads ${threads} -max_target_seqs '500' -word_size ${wordsize} -culling_limit ${cullinglimit}
    #cat `find ${blastdbdir}/ -maxdepth 1 -mindepth 1 -name "*.fasta"` | seqkit grep -r -p ^${sample} -v  | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit ${cullinglimit}


#outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen'
#sample='H112800346_9309cef4-089e-4b7f-938f-d6555e74e174'
#database='/home/alex/Documents/testing/ATCGpipeline/mankpc/samplelevel/subsetoutput/splitfastas/H112800346_9309cef4-089e-4b7f-938f-d6555e74e174_db'
#blastoutput="${1}/blast/${sample}_plasmidalignments.tsv"
#cat ${query} | seqkit grep -r -p ^${sample} -v  | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task 'blastn' -num_threads ${threads} -max_target_seqs '500' -word_size ${wordsize} -culling_limit ${cullinglimit}


# while IFS=$'\t' read -r -a line
# do
#     sample="${line[0]}"
#     database="${line[1]}"
#     blastoutput="${1}/blast/${sample}_plasmidalignments.tsv"
#     cat ${query} | seqkit grep -r -p ^${sample} -v | python runblast.py ${database} ${blastoutput} ${threads} ${evalue} ${wordsize}
# done < ${filepathinfo}

