#!/bin/bash
set -e
set -u
set -o pipefail

#arg[1] is filepathinfo.tsv file; arg[2] is threads; arg[3] is sourcedir
filepathinfo=${1}
threads=${2}
sourcedir=${3}

cat ${filepathinfo} | python ${sourcedir}/makeparallelcommand.py ${sourcedir} | parallel -k -j 1 "{}"