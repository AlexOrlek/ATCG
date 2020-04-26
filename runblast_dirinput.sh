#!/bin/bash
set -e
set -u
set -o pipefail

#argv[1] is filepath to pipeline output folder for a given dataset
#argv[2] is query sequence filepaths
#argv[3] is subject sequence filepaths (the same as query in the case of allvallpairwise blasttype)
#argv[4] is evalue #1e-8 default
#argv[5] is wordsize #38 default
#argv[6] is blast task
#argv[7] is threads
#argv[8] is blast bi-directional boolean flag (default False)
#argv[9] is blast type (-s flag or -1/-2 flags)

outputpath=${1}
sourcedir=${2}
filepathinfoquery=${3}
filepathinfosubject=${4} #subject sequence blast database; format: sample \t filepath to fasta \t filepath to database
evalue=${5}
wordsize=${6}
task=${7}
cullinglimit=${8}
threads=${9}
bidirectionalblast=${10} #default is False
blasttype=${11}  #allvallpairwise (-s flag) or pairwise (-1/-2 flags)

mkdir -p ${outputpath}/blast
mkdir -p ${outputpath}/output

if [ ${blasttype} != 'pairwiserun2' ]; then #prevents writing to file twice in the case of pairwise blasttype
    > ${outputpath}/allsubjects.txt  #used in reformatblastoutput.sh
    sort -k1,1V -o ${filepathinfoquery} ${filepathinfoquery}  #-k command must precede -o command
    sort -k1,1V -o ${filepathinfosubject} ${filepathinfosubject} 
fi

if [ ${blasttype} == 'allvallpairwise' ]; then
    if [ ${bidirectionalblast} == 'False' ]; then
        numsamples=$( cat ${filepathinfosubject} | wc -l )
        counter=0
        while IFS=$'\t' read -r -a line
        do
            let counter=$counter+1
            if [[ $counter -eq $numsamples ]]; then
                continue
            fi
            sample="${line[0]}"
            fastafile="${line[1]}"
            database="${line[2]}"
            echo "${sample}" >> ${outputpath}/allsubjects.txt
            blastoutput="${outputpath}/blast/${sample}_alignments.tsv"
            if [[ $counter -gt 1 ]]; then
                cat ${filepathinfoquery} | cut -f2 | grep -v -F ${fastafile} | grep -v -F -f <(printf "%s\n" "${excludesamples[@]}") | python ${sourcedir}/editfastaheaders.py 'stdin' | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit ${cullinglimit}
                excludesamples=(${excludesamples[@]} ${fastafile})
            else
                cat ${filepathinfoquery} | cut -f2 | grep -v -F ${fastafile} | python ${sourcedir}/editfastaheaders.py 'stdin' | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit ${cullinglimit}
                excludesamples=(${fastafile})
            fi
        done < ${filepathinfosubject}
    else
        while IFS=$'\t' read -r -a line
        do
            sample="${line[0]}"
            fastafile="${line[1]}"
            database="${line[2]}"
            echo "${sample}" >> ${outputpath}/allsubjects.txt
            blastoutput="${outputpath}/blast/${sample}_alignments.tsv"
            cat ${filepathinfoquery} | cut -f2 | grep -v -F ${fastafile} | python ${sourcedir}/editfastaheaders.py 'stdin' | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit ${cullinglimit}
        done < ${filepathinfosubject}
    fi
else   #blasttype==pairwise (or pairwiserun2 in the case of bidirectionalblast second run)
    while IFS=$'\t' read -r -a line
    do
        sample="${line[0]}"
        fastafile="${line[1]}"
        database="${line[2]}"
        echo "${sample}" >> ${outputpath}/allsubjects.txt
        blastoutput="${outputpath}/blast/${sample}_alignments.tsv"
        cat ${filepathinfoquery} | cut -f2 | python ${sourcedir}/editfastaheaders.py 'stdin' | blastn -db ${database} -out ${blastoutput} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit ${cullinglimit}
    done < ${filepathinfosubject}
fi


if [ ${blasttype} != 'pairwiserun2' ]; then #prevents writing to file twice in the case of pairwise blasttype
    > ${outputpath}/blastsettings.txt
    echo "BLAST task: ${task}" >> ${outputpath}/blastsettings.txt
    echo "e-value cutoff: ${evalue}" >> ${outputpath}/blastsettings.txt
    echo "word size: ${wordsize}" >> ${outputpath}/blastsettings.txt
    echo "bi-directional BLAST: ${bidirectionalblast}" >> ${outputpath}/blastsettings.txt
fi

