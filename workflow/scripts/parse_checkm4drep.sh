#! /bin/bash

# extract 1st, 12th and 13th columns (genome names, completness and contamination, respectevely)
# make a csv file
# add headers
# $1 is input table
# $2 is output file

new_head="genome,completeness,contamination"
echo $new_head > $2
awk '{if (NR!=1) {print $1 ".fa," $13 "," $14}}' $1 >> $2