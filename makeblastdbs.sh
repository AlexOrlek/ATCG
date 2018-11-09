#!/bin/bash
set -e
set -u
set -o pipefail

#arg[1] is fasta directory; arg[2] is fastafilepaths.tsv file; arg[3] is threads

mkdir -p ${1}

cat ${2} | cut -f2 | python ${4}/removeextension.py | parallel -k -j ${3} "makeblastdb -dbtype nucl -in {}.fasta -out {}_db"

