#!/bin/bash
# Title: trim_polyphred.sh
# Version: 1.2
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2015-06-24
# Modified in: 2024-07-09
# License: GPL v3



#======#
# Aims #
#======#

aim="Generate a summary table of genotypes, reads and scores for each sample from a polyphred output file."



#==========#
# Versions #
#==========#

# v1.2 - 2024-07-09: correct bug with the progress bar / improve speed
# v1.1 - 2021-05-05: correct bug in trap / correct bug in sort / correct bug in field differences between SNP and indel report
# v1.0 - 2021-02-07: add functions and options / rewrite input processing / remove temporary files
# v0.0 - 2015-06-24: creation

version=$(grep -i -m 1 "version" "$0" | cut -d ":" -f 2 | sed "s/^ *//g")



#===========#
# Functions #
#===========#

# Usage message
function usage {
    echo -e "
    \e[32m ${0##*/} \e[00m -i|--input file -o|--output file -p|--positions integer -d|--delimiter value -h|--help

Aim: $aim

Version: $version

Options:
    -i, --input     path of the polyphred report
    -o, --output    path of the output file
    -p, --positions positions within sample names to identify samples uniquely. Must be one or two integers.
                        - This can be used similarly to the -s option of polyphred.
                        - If -d is set, positions correspond to field delimitation.
    -d, --delimiter field delimiter to help identifying samples uniquely
    -h, --help      this message
    "
}


# Info message
function info {
    if [[ -t 1 ]]
    then
        echo -e "\e[32mInfo:\e[00m $1"
    else
        echo -e "Info: $1"
    fi
}


# Warning message
function warning {
    if [[ -t 1 ]]
    then
        echo -e "\e[33mWarning:\e[00m $1"
    else
        echo -e "Warning: $1"
    fi
}


# Error message
## usage: error "message" exit_code
## exit code optional (no exit allowing downstream steps)
function error {
    if [[ -t 1 ]]
    then
        echo -e "\e[31mError:\e[00m $1"
    else
        echo -e "Error: $1"
    fi

    if [[ -n $2 ]]
    then
        exit $2
    fi
}


# Dependency test
function test_dep {
    which $1 &> /dev/null
    if [[ $? != 0 ]]
    then
        error "Package $1 is needed. Exiting..." 1
    fi
}


# Progress bar
## Usage: ProgressBar $mystep $myend
function ProgressBar {
    if [[ -t 1 ]]
    then
        # Process data
        let _progress=(${1}*100/${2}*100)/100
        let _done=(${_progress}*4)/10
        let _left=40-$_done
        # Build progressbar string lengths
        _fill=$(printf "%${_done}s")
        _empty=$(printf "%${_left}s")

        # Build progressbar strings and print the ProgressBar line
        # Output example:
        # Progress : [########################################] 100%
        #printf "\rProgress : [${_fill// /=}${_empty// / }] ${_progress}%%"
        printf "\r\e[32mProgress:\e[00m [${_fill// /=}${_empty// / }] ${_progress}%%"

        [[ ${_progress} == 100 ]] && echo ""
    fi
}


# Clean up function for trap command
## Usage: clean_up file1 file2 ...
function clean_up {
    rm -rf $@
    echo ""
    exit 1
}



#==============#
# Dependencies #
#==============#

test_dep sed



#===========#
# Variables #
#===========#

set -e

# Options
while [[ $# -gt 0 ]]
do
    case $1 in
        -i|--input     ) input="$2" ; shift 2 ;;
        -o|--output    ) output="$2" ; shift 2 ;;
        -p|--positions ) pos="$2" ; shift 2
                           while [[ ! -z "$1" && $(echo "$1"\ | grep -qv "^-" ; echo $?) == 0 ]]
                           do
                               pos+="-$1"
                               shift
                           done ;;
        -d|--delimiter ) delim="$2"    ; shift 2 ;;
        -h|--help      ) usage ; exit 0 ;;
        *              ) error "Invalid option: $1\n$(usage)" 1 ;;
    esac
done


# Check the existence of obligatory options
[[ -z "$input" ]]  && error "The option input is required. Exiting...\n$(usage)" 1
[[ -z "$output" ]] && error "The option output is required. Exiting...\n$(usage)" 1

# Check if output exists
[[ -e "$output" ]] && error "$output exists. Exiting..." 1

# Check position
[[ -z "$pos" ]] && pos="1-"
[[ $(awk -F "-" 'END {print NF}' <<< "$pos") -gt 2 ]] && error "More than two positions entered. Exiting..." 1
[[ $(sed "s/[-0-9]//g" <<< "$pos" | grep .) ]] && error "Position must be only integer. Exiting..." 1



#============#
# Processing #
#============#

# Trap
trap "clean_up \"$output\"" SIGINT SIGTERM    # Clean_up function to remove tmp files
wait

# Isolating genotype section from polyphred output
ppo=$(sed -n "/BEGIN_GENOTYPE/,/END_GENOTYPE/{//d;p}" "$input" | sed "s/ * /\t/g")

# Sites with genotypes
sites=$(cut -f 1 <<< "$ppo" | sort -n | uniq)

# Output header
header="Sample\t"
for s in $sites
do
    s_hd=$(printf "$s-%s\t" GT reads score)
    header+="$s_hd"
done
mytable="$header"

# Adjust field numbers (differences between indel and SNP)
if [[ $(sed -n "/BEGIN_COMMAND_LINE/,/END_COMMAND_LINE/{//d;p}" "$input" | grep " -indel ") ]]
then
    spl_fd=3
    sc_fd=6
    gt_fd=4-5
else
    spl_fd=4
    sc_fd=7
    gt_fd=5-6
fi

# Listing samples
if [[ -z "$delim" ]]
then
    samples=($(cut -f $spl_fd <<< "$ppo" | cut -c $pos | sort | uniq))
else
    samples=($(cut -f $spl_fd <<< "$ppo" | cut -d "$delim" -f $pos | sort | uniq))
fi

[[ -z "$samples" ]] && error "No sample detected. Exiting..." 1

# Analyzing each sample
for ((i = 0 ; i < ${#samples[@]} ; i++))
do

    # Counter
    j=$(($i + 1))
    ProgressBar $j ${#samples[@]} || :

    sample="${samples[$i]}"

    # Isolate sample specific block
    ppo_spl=$(awk -v spl_fd=$spl_fd -v sample="$sample" ' $spl_fd ~ sample' <<< "$ppo")

    # Output line
    myline="$sample"

    for s in $sites
    do
        if $(grep -q "$s" <(cut -f 1 <<< "$ppo_spl"))
        then
            ppo_blk=$(awk -v s=$s '$1 == s' <<< "$ppo_spl" |  grep "$sample")

            # Nb of read (forward and reverse or just one of them)
            nb_read=$(wc -l <<< "$ppo_blk")

            # Score (the minimum for a given sample)
            myscore=$(cut -f $sc_fd <<< "$ppo_blk" | sort -n | head -1)

            # Genotype
            mygt=$(cut -f $gt_fd <<< "$ppo_blk" | tr "\t" "/")
            if [[ $(echo "$mygt" | wc -l) == 1 || $(echo "$mygt" | wc -l) == 2 && $(echo "$mygt" | sed -n "1p") == $(echo "$mygt" | sed -n "2p") ]]
            then
                mygt=$(echo "$mygt" | sed -n "1p")
            else
                warning "Several genotype detected for sample $sample at site $s. Skipping..."
                mygt="NA"
            fi
        else
            nb_read="NA"
            myscore="NA"
            mygt="NA"
        fi

        # Table line
        myline+="\t$mygt\t$nb_read\t$myscore"

    done

    mytable+="\n$myline"

done

echo -e "$mytable" >> "$output"

exit
