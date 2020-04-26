outputpath=${1}
sourcedir=${2}
filepathinfo=${3}
comparisonsfile=${4}
evalue=${5}
wordsize=${6}
task=${7}
cullinglimit=${8}
threads=${9}
bidirectionalblast=${10} #default is False
blasttype=${11} #sequencepairs or sequencepairsrun2 (second run of bidirectional blast)

mkdir -p "${outputpath}/blast"
mkdir -p "${outputpath}/output"

declare -A comparisonsdict
counter=0
while IFS=$'\t' read -r -a line
do
    let counter=$counter+1
    if [ ${blasttype} != 'sequencepairsrun2' ]; then
        filepath1="${line[0]}"
        filepath2="${line[1]}"
    else
        filepath1="${line[1]}"
        filepath2="${line[0]}"
    fi

    if [[ -v "${comparisonsdict[$filepath2]}" ]]; then   #key not yet added - initialise array key
        comparisonsdict["${filepath2}"]="${filepath1}"
    else #append to array key
        myarray=${comparisonsdict["${filepath2}"]}
        myarray=$(echo ${myarray} ${filepath1})
        comparisonsdict[${filepath2}]=${myarray}
    fi
done < ${comparisonsfile}



#now need to run blast, looping through keys as subjects and using values as queries; for bidirectional, on run2 a reciprocal dictionary will be created above

while IFS=$'\t' read -r -a line
do
    sample="${line[0]}"
    filepath="${line[1]}"
    database="${line[2]}"
    blastoutput="${outputpath}/blast/${sample}_alignments.tsv"
    #echo "${comparisonsdict[${filepath}]}" | tr ' ' '\n'
    #echo 'database:' ${database}
    echo "${comparisonsdict[${filepath}]}" | tr ' ' '\n' | python ${sourcedir}/editfastaheaders.py 'stdin' | blastn -db ${database} -evalue ${evalue} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task ${task} -num_threads ${threads} -word_size ${wordsize} -culling_limit ${cullinglimit} >> ${blastoutput}
done < ${filepathinfo}


if [ ${blasttype} != 'sequencepairsrun2' ]; then #prevents writing to file twice in the case of pairwise blasttype
    > ${outputpath}/blastsettings.txt
    echo "BLAST task: ${task}" >> ${outputpath}/blastsettings.txt
    echo "e-value cutoff: ${evalue}" >> ${outputpath}/blastsettings.txt
    echo "word size: ${wordsize}" >> ${outputpath}/blastsettings.txt
    echo "bi-directional BLAST: ${bidirectionalblast}" >> ${outputpath}/blastsettings.txt
fi


