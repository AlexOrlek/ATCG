#!/bin/bash

#$ -S /bin/bash
#$ -N driver-b
#$ -e b.err-driver
#$ -o b.out-driver
#$ -P bag.prjb
#$ -q short.qb -cwd -V
#$ -t 1-18

i=$(expr $SGE_TASK_ID)
sample=$(cat ./../../../../../data/apha/fastqs_longread/finalgenomes/output_plasmids/individual/samplenames.txt | sed -n ${i}p)
python runpipeline.py -s1 ./../../../../../data/apha/fastqs_longread/finalgenomes/output_plasmids/individual/${sample}.fasta -s2 ./../../../../../data/apha/fastqs_longread/finalgenomes/output_plasmids/allplasmids.fa -o ./output/${sample} -t 1 -b 100 --breakpoint --alnlenstats
