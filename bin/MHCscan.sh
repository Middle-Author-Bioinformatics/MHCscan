#! /bin/bash

function usage() {
    cat <<USAGE

    Usage: $0 [-i input] [-o out] [-t thr] [-x server] [--prot]

    Options:
        -i, --input:  input genome or protein sequences in FASTA format
        -o, --output: output directory namee
        -g, --gff:    GFF file corresponding to protein sequences (ignore if unannotated contigs are provided)
        -d, --id:     attribute name of the locus tag in the GFF file
        -t, --thr:    parallel threads to use

        --prot        input file contains protein sequences
USAGE
    exit 1
}

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

THR=2
OUTPUT="MHCs.csv"
PROT=false
ID="ID"
while [ "$1" != "" ]; do
    case $1 in
    --prot)
        PROT=true
        ;;
    -i | --input)
        shift
        INPUT=$1
        ;;
    -g | --gff)
        shift
        GFF=$1
        ;;
    -d | --id)
        shift
        ID=$1
        ;;
    -o | --output)
        shift
        OUTPUT=$1
        ;;
    -t | --thr)
        shift
        THR=$1
        ;;
    -h | --help)
        usage
        ;;
    *)
        usage
        exit 1
        ;;
    esac
    shift
done

if [[ ${INPUT} == "" ]]; then
    echo 'Please provide input file';
    exit 1;
fi

if [[ ${PROT} != true ]]; then
    echo "--Running Prodigal"
    prodigal -i ${INPUT} -a ${INPUT%.*}.faa
    INPUT=${INPUT%.*}.faa
else
    echo '--Skipping Prodigal'
fi

remove_asterisk.py -f ${INPUT} -o ${INPUT%.*}-fixed.faa

phobius.pl -short ${INPUT%.*}-fixed.faa > ${INPUT%.*}.phobius

if [[ ${PROT} == true ]]; then
    if [[ ${GFF} == "" ]]; then
        echo 'Please provide GFF file';
        exit 1;
    fi


    cytoscan.py -p ${INPUT%.*}.phobius -f ${INPUT%.*}-fixed.faa -o ${OUTPUT} -g ${GFF} -i ${ID} --add_cds

    eval "$(conda shell.bash hook)"
    conda activate hhsuite_env
    mhc_annotate.py -i ${OUTPUT} -o ${OUTPUT}-annotated.csv -g ${GFF} --add_cds -id ${ID} -t ${THR}
#    mv ${OUTPUT}-annotated.csv ${OUTPUT}
    conda deactivate


else
    cytoscan.py -p ${INPUT%.*}.phobius -f ${INPUT%.*}-fixed.faa -o ${OUTPUT} --prodigal

    eval "$(conda shell.bash hook)"
    conda activate hhsuite_env
    mhc_annotate.py -i ${OUTPUT} -o ${OUTPUT}-annotated.csv -id ${ID} -t ${THR}
#    mv ${OUTPUT}-annotated.csv ${OUTPUT}
    conda deactivate
fi

#rm ${INPUT%.*}-fixed.faa
#rm ${INPUT%.*}.phobius




