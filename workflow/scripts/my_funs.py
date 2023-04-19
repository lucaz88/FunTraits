#!/usr/bin/env python
# -*- coding: utf-8
 
     
## find inputs & save path and filenames as table
def find_inputs(dir, filename_pattern, double_ext=False, uniq_PE=False):
    import pandas as pd
    import glob
    import os
    import re
    
    file_list = glob.glob(os.path.join(dir, "**/*"+filename_pattern), recursive=True)
    filename_list = [os.path.basename(i) for i in file_list]
    # filename_noext_list = [filename[:filename.index('.')] for filename in filename_list] # chop at the first dot position, i.e. no "dots" allowed in filename
    if double_ext: # e.g. '.tar.gz', '.fq.gz', 'fastq.gz'
        # filename_noext_list = [re.sub('.tar.[a-z]+$', '', filename) for filename in filename_list]
        # ext_list = [re.sub('.*.tar.', 'tar.', filename) for filename in filename_list]
        filename_noext_list = [filename.split(".")[0] for filename in filename_list]
        ext_list = ['.'+'.'.join(filename.split(".")[1:]) for filename in filename_list]
    else:
        filename_noext_list = [os.path.splitext(filename)[0] for filename in filename_list]
        ext_list = [os.path.splitext(filename)[1] for filename in filename_list]
    path_list = [os.path.dirname(i) for i in file_list]
    
    df = pd.DataFrame ([file_list, path_list, filename_list, filename_noext_list, ext_list]).transpose()
    df.columns = ['full','dir', 'file', 'file_noext', 'ext']
    df = df.sort_values(by=['file'])

    if uniq_PE:
        df['uniq_PE'] = [re.sub('_R[12]$|_R[12]_.*', '', i) for i in df['file_noext']]
        df = df[~df.duplicated(subset=['uniq_PE'])]
        df = df.reset_index(drop=True)
    
    return df
## ---


## convert column of pandas DataFrame to list
def col2list(input_table, col_id=False):
    if not col_id or col_id == "last":
        col_id = input_table.shape[1]-1
    if isinstance(col_id, int) and col_id <= (input_table.shape[1]-1):
        col_list = input_table.iloc[:, col_id].values.flatten().tolist()
    elif isinstance(col_id, str) and col_id in input_table:
        col_list = input_table.loc[:, col_id].values.flatten().tolist()
    else:
        print("Something is wrong with the specified column. Use existing col names or ids\n") 
        exit() 
    
    return col_list
## ---


## set time resource of a rule
def time_per_gb(wildcards, input):
    half_gb = input.size_mb/500

    if half_gb < 0.5:
        time = ("00:%02d:00" % (59))
    else:
        time = ("%02d:00:00" % (round(half_gb, 0)))

    return time
## ---