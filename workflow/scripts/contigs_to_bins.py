#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__      = "Luca Zoccarato"


import argparse
import os
import glob
import re


def function_arguments():
    parser = argparse.ArgumentParser(description="Parse bin fasta files to create a (tab separated) contigs to bin table.\
        Can digest multiple binning strategy at once, and return a file containing comma-separated paths of each bin table generated", usage="%(prog)s [options]")
    parser.add_argument("-b, --bin_dirs", 
                        metavar="", dest="bin_dirs", nargs='+',
                        help="list of directories where to find the bin fasta files", type=str)
    parser.add_argument("-o, --out_dir", 
                        metavar="''", dest="out_dir",
                        help="output directory for the contigs to bin tables", type=str)
    args = parser.parse_args()
    return(args)


def contigs_to_bins():
    #! parse function arguments
    args = function_arguments()
    bin_dirs = [re.sub("/$", "", i) for i in args.bin_dirs]
    out_dir = args.out_dir
    
    out_paths = []
    os.makedirs(out_dir, exist_ok=True)

    for i_dir in bin_dirs:        
        #! find bin files
        i_bins = glob.glob(os.path.join(i_dir, "*.fa"), recursive=True)
        if len(i_bins) == 0:
            i_bins = glob.glob(os.path.join(i_dir, "*.fasta"), recursive=True)
        
        if len(i_bins) > 0: # skip binner when no bin has been retrieved  
            #! create output file
            out_tab = os.path.join(out_dir, os.path.basename(i_dir)+"_contigs2bins.tsv")
            out_f = open(out_tab, 'w')
            
            #! loop through bins
            for ij_bin in i_bins:
                ij_bin_name = os.path.basename(ij_bin)
                ij_bin_name = os.path.splitext(ij_bin_name)[0]
                
                #! use fasta headers as link to contigs
                with open(ij_bin, 'r') as ij_fasta:
                    for line in ij_fasta:
                        if line.startswith('>'):
                            ij_line = re.sub(">|\n", "", line)
                            ij_line = re.sub(" .*", "", ij_line)
                            out_f.write('%s\t%s\n' % (ij_line, ij_bin_name))

            out_f.close()
            out_paths.append(out_tab)
    
    #! save list table paths in console
    out_paths = ",".join(out_paths)
    paths_file = os.path.join(out_dir, "contig2bin_list.txt")
    paths_f = open(paths_file, 'w')
    paths_f.write('%s\n' % (out_paths)) 
    paths_f.close()


if __name__ == "__main__":
    contigs_to_bins()


    
    
    



# def contigs_to_bins(bin_dirs, out_dir):
#     import os
#     import glob
#     import re
#     import pandas as pd
    
#     out_files = []
#     for i_dir in bin_dirs:
#         #! create output file
#         out_file = os.path.join(out_dir, os.path.basename(i_dir)+"_contigs2bins.tsv")
#         os.makedirs(os.path.dirname(out_file), exist_ok=True)
#         out = open(out_file, 'w')
        
#         #! find bin files
#         i_bins = glob.glob(os.path.join(i_dir, "*.fa"), recursive=True)
#         if len(i_bins) == 0:
#             i_bins = glob.glob(os.path.join(i_dir, "*.fasta"), recursive=True)
        
#         # i_conts2bins = pd.DataFrame()
#         for ij_bin in i_bins:
#             ij_bin_name = os.path.basename(ij_bin)
#             ij_bin_name = os.path.splitext(ij_bin_name)[0]
#             # ij_contigs = []
            
#             with open(ij_bin, 'r') as ij_fasta:
#                 for line in ij_fasta:
#                     if line.startswith('>'):
#                         ij_line = re.sub(">|\n", "", line)
#                         out.write('%s\t%s\n' % (ij_line, ij_bin_name))
#             #             ij_contigs.append(re.sub(">|\n", "", line))
#         out.close()
#         out_files.append(out_file)
#             # i_bin_name = [i_bin_name] * len(ij_contigs)
#             # i_df = pd.DataFrame([i_bin_name, ij_contigs]).transpose()
#             # i_conts2bins = i_conts2bins.append(i_df)
    
#     return out_files