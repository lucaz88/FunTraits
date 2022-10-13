#! /bin/bash

# filter table using values in the first column and directory containig filenames to match
# $1 table to filter
# $2 file with filenames to retain
# $3 output file (optional, otherwise it overwrites input file)

header=$(awk 'NR==1{print $1}' $1)
filename_tokeep=$(ls $2/* | while read file; do basename $file; done)
echo $header $filename_tokeep | tr ' ' '\n' > tmp_vector.tsv

grep -f tmp_vector.tsv $1 > tmp_out.txt

if [[ -v $3 ]];
then
    mv tmp_out.txt $3
else
    mv tmp_out.txt $1
fi

rm tmp_vector.tsv