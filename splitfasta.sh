#!/bin/bash
set -e
set -u
set -o pipefail

#arg[1] is output path; arg[2] is multifasta file; arg[3] is threads            

#make output directories                                                        
mkdir -p ${1}/blast
mkdir -p ${1}/splitfastas
mkdir -p ${1}/output

#split multifasta into individual genomes/samples                               
seqkit split ${2} -i --id-regexp "^(.+?)\|" -O ${1}/splitfastas --quiet --threads ${3} --force
