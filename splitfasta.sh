#!/bin/bash
set -e
set -u
set -o pipefail


#arg[1] is output path; arg[2] is fastadir; arg[3] is fasta; arg[4] is threads
   
#make output directories                                                        
mkdir -p ${1}/blast
mkdir -p ${1}/output
mkdir -p ${2}

#split multifasta into individual genomes/samples                               
seqkit split ${3} -i --id-regexp "^(.+?)\|" -O ${2} --quiet --threads ${4} --force
