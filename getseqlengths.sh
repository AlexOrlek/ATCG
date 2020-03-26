#!/bin/bash
set -e
set -u
set -o pipefail

#args[1] is outputpath args[2] is blast type (all v all or non all v all); arg[3] is fasta filepath; arg[4] is fasta filepath 2 if running non all v all

if [ "${2}" == 'allvallpairwise' ]; then
    cat `find ${3}/ -maxdepth 1 -mindepth 1 -name "*.fasta"` | bioawk -c fastx '{print $name,length($seq)}' | sort -k1,1V | python ${4}/formatseqlen.py > ${1}/seqlengths.tsv
else
    cat `find ${3}/ ${4}/ -maxdepth 1 -mindepth 1 -name "*.fasta"` | bioawk -c fastx '{print $name,length($seq)}' | sort -k1,1V | python ${5}/formatseqlen.py > ${1}/seqlengths.tsv
fi


