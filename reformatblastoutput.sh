#!/bin/bash
set -e
set -u
set -o pipefail


> ${1}/included.txt
#$1 is filepath to pipeline output folder

#reformat blast output file
samples=($(cut -f1 ${2} | sort | uniq)) #${2} is subject sequence blast database file; format: sample \t filepath to database

for sample in ${samples[@]}; do
    mkdir -p ${1}/blast/${sample}
    cat ${1}/blast/${sample}_alignments.tsv | awk '{OFS="\t"; if($9 < $10) {qcovprop=($13/100); qcovhspprop=($14/100); $13=qcovprop; $14=qcovhspprop; print $0"\t""+"} else {qcovprop=($13/100); qcovhspprop=($14/100); $13=qcovprop; $14=qcovhspprop; sstart=$10; send=$9; $9=sstart; $10=send; print $0"\t""-"}}' | tee ${1}/blast/${sample}/alignments.tsv | if [ $(wc -l) -gt 0 ]; then echo "${sample}" >> ${1}/included.txt; fi 
done


python getexcludednames.py ${1} ${2}


#NOTES
#if($9 < $10)... means if sstart < send reassign qcov/qcovhsp as proportions; print line and append "+"; else do same but also swap sstart/send columns and append '-' instead of '+'
