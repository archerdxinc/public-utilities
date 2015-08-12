#!/bin/bash
#===================================================================================
#
# FILE: merge_nextseq_lanes.sh
#
# DESCRIPTION: Merges mutiple NextSeq lanes into a single Fastq
#
# OPTIONS: see function ’usage’ below
# REQUIREMENTS: None
# BUGS: ---
# NOTES: ---
# AUTHOR: Aaron Berlin, aberlin@archerdx.com
# COMPANY: ArcherDx
# VERSION: 1.0
# CREATED: 07JAN2015
# REVISION:
#===================================================================================


#usage
usage()
{
cat << EOF >> /dev/stderr
 usage: $0 options

 OPTIONS:
   -h Show this message
   -x Extension to be added (default: combined)
   -d Directory to merge

EOF
}

#workflow
workflow()
{

    cd $directory &&
    read_types=("R1 R2")
    for read in $read_types
    do
        ls *L00?*${read}*gz |
        grep -v combined | # Ignore previously combined fastqs
        grep -v Undet | # Ignore the Undetermined
        awk '{if(NR%4==0){print }else{printf "%s\t",$0}}' |

        while read group
        do
            curr=`echo $group |cut -d ' ' -f1`
            `cat $group > ${curr%.fastq.gz}.${extension}.fastq.gz`
        done
    done

}

#parse_options
parse_options()
{

    #Extension to be added
    extension="combined"

    # Directory to work in
    directory=./

    while getopts "hx:d:" OPTION
    do
    case $OPTION in
        h)
        usage
        exit 1
        ;;
        x)
        extension=$OPTARG
        ;;
        d)
        directory=$OPTARG
        ;;
    esac
    done

    return 1
}

#Main Function
main()
{
    #Parse the command line options
    # "$*" are the command line arguments
    parse_options $*

    #Ensure halting if non-zero exit status detected
    set -e
    set -o pipefail

    #Call workflow function
    workflow
}

#Call the main function, pass in command line arguments
main $* 2>> /dev/stderr
