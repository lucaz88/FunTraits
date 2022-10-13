#! /bin/bash

# select 1st (contig names) and 3rd (mapped reads) columns
# remove last row (containing unmapped reads)
# $1 is input table
# $2 is output file

awk '{print $1 "\t" $3}' $1 | sed '/^\*/d' > $2