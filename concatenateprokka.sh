#!/bin/bash
set -e
set -u
set -o pipefail

inputdir=${1}
#outdir=${2} #stdout to python script
gfffiles=( $(find ${inputdir} -maxdepth 1 -mindepth 1 -type f -name "*.gff") )
#echo ${gfffiles[@]}

endstring='##FASTA'

for file in ${gfffiles[@]}
do
    #echo ${file}
    sedcommand="/${endstring}/"',$d'
    #echo "${sedcommand}"
    sed -e "${sedcommand}" ${file}
done
