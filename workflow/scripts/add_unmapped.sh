#! /bin/bash

# add count of unmapped reads
# $1 is output of 'samtools idxstats'
# $2 is output of 'samtools idxstats' after filering BAM with 'filtersam' script

# compute unmapped reads in filtered file
old_unmap=$(tail -n 1 $1 | cut -f 4)
old_map=$(cat $1 | cut -f 3 | paste -sd+ - | bc)
new_map=$(cat $2 | cut -f 3 | paste -sd+ - | bc)
new_unmap=$(expr $old_unmap + $old_map - $new_map)

# remove 
head -n -1 $2 > temp.txt
mv temp.txt $2

# append new line
new_line="*       0       0       "
echo "$new_line"$new_unmap >> $2