#! /bin/bash

# rename sequentially files in provided directory as "basename_string_XXXX"
# $1 dir cotaining files to rename
# $2 extension of file to be renamed, e.g. ".fa"
# $3 basename string to add at the beginning of each new name
# $4 2-columns table with old filenames and new filenames
# $5 logic, retain root of filenames, false by default

def_5=false
keep_root=${5:-$def_5}

cnt=1
touch $4

if $keep_root
then
  ls $1/*$2 | while read file; do
    name_root=$(basename $file | sed -e 's/-.*//')
    new_name=$1/$(printf "%s_%04d%s" $name_root "$cnt" "$2") #04 pad to length of 4
    mv $file $new_name
    printf "%s\t%s\n" $(basename $file) $(basename $new_name) >> $4
    let cnt=cnt+1
  done
else
  ls $1/*$2 | while read file; do
    new_name=$1/$(printf "%s_%04d%s" $3 "$cnt" "$2") #04 pad to length of 4
    mv $file $new_name
    printf "%s\t%s\n" $(basename $file) $(basename $new_name) >> $4
    let cnt=cnt+1
  done
fi