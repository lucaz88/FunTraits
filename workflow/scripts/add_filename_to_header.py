#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os
from posixpath import dirname
import re
import glob


def function_arguments():
    parser = argparse.ArgumentParser(description="Add the fasta filename to all headers: \
        >contig_1 ==> >filename__contig_1.", usage="%(prog)s [options]")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i, --in_fasta", 
                        metavar="", dest="in_fasta", nargs='+',
                        help="single fasta file or space separated list of fasta files to parse", type=str)
    group.add_argument("-d, --in_dir", 
                        metavar="", dest="in_dir", nargs=1,
                        help="or path to dir containing multiple fasta files (! with .fa extesnion)", type=str)
    parser.add_argument("-x, --file_ext", 
                        metavar="", dest="file_ext", nargs=1,
                        help="file extension to look up when using '-d'", type=str)
    parser.add_argument("-o, --out_fasta", 
                        metavar="", dest="out_fasta", nargs='+', default=None,
                        help="name of parsed fasta file/s. \
                                If omitted input file will be overwritten. \
                                If multiple fasta files are provided with only one output name, fasta files will be concatenated.", type=str)
    args = parser.parse_args()
    return(args)


def add_filename_to_header():
    #! parse function arguments
    args = function_arguments()
    in_fasta = args.in_fasta
    in_dir = args.in_dir[0]
    file_ext = args.file_ext[0]
    out_fasta = args.out_fasta
    
    
    if in_dir:
        in_fasta = glob.glob(os.path.join(in_dir, "**/*"+file_ext), recursive=True)

    
    try:
        if out_fasta is None:
            for i_in in in_fasta:
                i_out = os.path.join(os.path.dirname(i_in), "tmp.fa")
                if len(os.path.dirname(i_out)) > 0: # check that path is different from PWD
                    os.makedirs(os.path.dirname(i_out), exist_ok=True)
                i_out_f = open(i_out, 'w')
                i_in_name = os.path.basename(i_in)
                with open(i_in, 'r') as i_in_lines:
                    for line in i_in_lines:
                        if line.startswith('>'):
                            line_new = re.sub(">", ">"+i_in_name+"__", line)
                            i_out_f.write(line_new)
                        else:
                            i_out_f.write(line)
                os.replace(i_out, i_in)     
                
        else:
            
            if len(in_fasta) >= 1 and len(out_fasta) >= 1  and len(in_fasta) == len(out_fasta):
                for i_in,i_out in zip(in_fasta, out_fasta):
                    if len(os.path.dirname(i_out)) > 0: # check that path is different from PWD
                        os.makedirs(os.path.dirname(i_out), exist_ok=True)
                    i_out_f = open(i_out, 'w')
                    i_in_name = os.path.basename(i_in)
                    with open(i_in, 'r') as i_in_lines:
                        for line in i_in_lines:
                            if line.startswith('>'):
                                line_new = re.sub(">", ">"+i_in_name+"__", line)
                                i_out_f.write(line_new)
                            else:
                                i_out_f.write(line)
                                
            elif len(in_fasta) > 1 and len(out_fasta) == 1: # concatenated output
                i_out = out_fasta[0]
                if len(os.path.dirname(i_out)) > 0: # check that path is different from PWD
                    os.makedirs(os.path.dirname(i_out), exist_ok=True)
                i_out_f = open(i_out, 'w')
                for i_in in in_fasta:
                    i_in_name = os.path.basename(i_in)
                    with open(i_in, 'r') as i_in_lines:
                        for line in i_in_lines:
                            if line.startswith('>'):
                                line_new = re.sub(">", ">"+i_in_name+"__", line)
                                i_out_f.write(line_new)
                            else:
                                i_out_f.write(line)
                        
            else:
                raise ValueError('Different number of input and output provided.')
    
    except (ValueError, IndexError):
        exit('Could not complete the request.')  
    


if __name__ == "__main__":
    add_filename_to_header()