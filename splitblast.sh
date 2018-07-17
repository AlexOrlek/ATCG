#!/bin/bash
set -e
set -u
set -o pipefail

> ${1}/samplenames.txt
#$1 is filepath to pipeline output folder

#reformat blast output file; output reformatted file into sample folder as queryplasmidalignments.tsv; initialise subjectplasmidalignments.tsv
while IFS=$'\t' read -ra line
do
    sample=${line[0]}
    mkdir -p ${1}/blast/${sample}
    cat ${1}/blast/${sample}_plasmidalignments.tsv | python filterselfselfhits.py | awk '{OFS="\t"; if($9 < $10) {qcovprop=($13/100); qcovhspprop=($14/100); $13=qcovprop; $14=qcovhspprop; print $0"\t""+"} else {qcovprop=($13/100); qcovhspprop=($14/100); $13=qcovprop; $14=qcovhspprop; sstart=$10; send=$9; $9=sstart; $10=send; print $0"\t""-"}}' | tee ${1}/blast/${sample}/plasmidalignments.tsv | if [ $(wc -l) -gt 0 ]; then echo "${sample}" >> ${1}/samplenames.txt; fi 
done < ${2} #plasmid blast database file; format: sample \t filepath to database


python getexcludedsamples.py ${1} ${2}


#NOTES
#if($9 < $10)... means if sstart < send reassign qcov/qcovhsp as proportions; print line and append "+"; else do same but also swap sstart/send columns and append '-' instead of '+'
