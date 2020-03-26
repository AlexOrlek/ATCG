#!/bin/bash
set -e
set -u
set -o pipefail

#arg[1] is fastafilepaths.tsv file; arg[2] is threads; arg[3] is sourcedir

cat ${1} | cut -f2 | python ${3}/removeextension.py | parallel -k -j ${2} "makeblastdb -dbtype nucl -in {}.fasta -out {}_db -logfile /dev/null"

