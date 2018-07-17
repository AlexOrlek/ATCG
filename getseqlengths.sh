#!/bin/bash
set -e
set -u
set -o pipefail

#args[1] is outputpath args[2] is query plasmid fasta

cat ${2} | bioawk -c fastx '{print $name,length($seq)}' > ${1}/seqlengths.tsv
