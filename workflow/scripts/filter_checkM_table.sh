#! /bin/bash

# filter table using values in the first column and directory containig filenames to match
# $1 table to filter
# $2 file with filenames to retain in 1st column and new names in 2nd (e.g. output of rename_sequentially.sh)
# $3 output file (optional, otherwise it overwrites input file)

header=$(awk 'NR==1{print $1}' $1)
echo $header > tmp_vector.tsv
cat $2 | awk '{ print $1 }' >> tmp_vector.tsv

echo "MAG_ids" > tmp_vector2.tsv
cat $2 | awk '{ print $2 }' >> tmp_vector2.tsv

grep -f tmp_vector.tsv $1 > tmp_out.txt
paste tmp_vector2.tsv tmp_out.txt > tmp_out2.txt

if [[ -v $3 ]];
then
    mv tmp_out2.txt $3
else
    mv tmp_out2.txt $1
fi

rm tmp_vector.tsv tmp_vector2.tsv tmp_out.txt