#!/bin/bash
set -e
set -u
set -o pipefail

#args[1] is outputpath args[2] is blast type (all v all or non all v all); arg[3] is fasta; arg[4] is fasta 2 if running non all v all

if [ "${2}" == 'allvallpairwise' ]; then
    cat ${3} | bioawk -c fastx '{print $name,length($seq)}' > ${1}/seqlengths.tsv
else
    cat ${3} ${4} | bioawk -c fastx '{print $name,length($seq)}' > ${1}/seqlengths.tsv
fi

> ${1}/output/distancestats.tsv
