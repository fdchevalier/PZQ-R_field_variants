#!/bin/bash
# Title: trim_polyphred.sh
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2015-06-24
# Modified in:
# Licence : GPL v3



#======#
# Aims #
#======#

aim="Read polyphred output file to generate corresponding genotype table that indicates for each sample if the variant is supported by one or both reads."



#==========#
# Versions #
#==========#

# v0.0 - 2015-06-24: creation



#===========#
# Variables #
#===========#

input="$1"
output="$2"

mydate=$(date +%N)
tmp="/tmp/${input##*/}_$mydate"


# Check if output exists
if [[ -e "$output" ]]
then
    echo "$output exists. Exiting..."
    exit 1
fi



#============#
# Processing #
#============#

# Cleaning polyphred output 
sed "s/ * /\t/g" "$input" > "$tmp"


# Get the position of the polymorphic sites
cat "$tmp" | cut -f 1 | sort | uniq | \
while read site
do

    # Create a tmp file for a given site then read each sample
    ## NB: sample name is cut to take the most common part between the two reads. Grep uses "." to make the difference between x.1 and x.10 (initial naming is x.1.y and x.10.y)
    awk -v site="$site" ' $1 == site {print $0}' "$tmp" > "$tmp-a"
    cut -f 4 "$tmp-a" | cut -d "." -f -3 | sort | uniq | \
    while read sample
    do
        # Nb of read (forward and reverse or just one of them)
        nb_read=$(grep -c "$sample\." "$tmp-a")

        # Score (the minimum for a given sample)
        myscore=$(grep "$sample\." "$tmp-a" | cut -f 7 | sort -n | head -1)

        # Genotype
        mygt=$(grep "$sample\." "$tmp-a" | cut -f 5-6 | tr "\t" "/")
        if [[ $(echo "$mygt" | wc -l) == 1 || $(echo "$mygt" | wc -l) == 2 && $(echo "$mygt" | sed -n "1p") == $(echo "$mygt" | sed -n "2p") ]]
        then
            mygt=$(echo "$mygt" | sed -n "1p")
        else
            echo "Different genotype for a same sample. Skipping..."
            mygt=""
        fi

        # Table line
        myline="$site\t$sample\t$mygt\t$myscore\t$nb_read"
        echo -e "$myline" >> "$output"

    done

done

