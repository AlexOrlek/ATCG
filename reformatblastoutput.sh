#!/bin/bash
set -e
set -u
set -o pipefail


> ${1}/included.txt
#$1 is filepath to pipeline output folder
bidirectionalblast=${4}

#reformat blast output file
samples=($(cut -f1 ${2} | sort -V | uniq)) #${2} is subject sequence blast database file; format: sample \t filepath to database

counter=0
for sample in ${samples[@]}; do
    let counter=$counter+1
    if [ ${bidirectionalblast} == 'False' ] && [ $counter -eq ${#samples[@]} ]; then
        continue
    fi
    mkdir -p ${1}/blast/${sample}
    cat ${1}/blast/${sample}_alignments.tsv | awk '{OFS="\t"; if($9 < $10) {qcovprop=($13/100); qcovhspprop=($14/100); $13=qcovprop; $14=qcovhspprop; print $0"\t""+"} else {qcovprop=($13/100); qcovhspprop=($14/100); $13=qcovprop; $14=qcovhspprop; sstart=$10; send=$9; $9=sstart; $10=send; print $0"\t""-"}}' | sed '1iqname\tsname\tpid\talnlen\tmismatches\tgapopens\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcov\tqcovhsp\tqlength\tslength\tstrand' | tee ${1}/blast/${sample}/alignments.tsv | if [ $(wc -l) -gt 0 ]; then echo "${sample}" >> ${1}/included.txt; fi
    rm ${1}/blast/${sample}_alignments.tsv
done


python ${3}/checkfornoblasthits.py ${1}

#NOTES
#if($9 < $10)... means if sstart < send reassign qcov/qcovhsp as proportions; print line and append "+"; else do same but also swap sstart/send columns and append '-' instead of '+'

#sed '1i ...' adds header line, but only if there is data to add line to so wc -l -gt 0 is still a test for whether there is data (as opposed to -gt 1 to account for header line)
