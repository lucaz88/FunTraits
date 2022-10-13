#!/usr/bin/env python
# -*- coding: utf-8


#!!! discontinued - see myPy_* !!!


## dummy vars for testing
from xmlrpc.client import boolean


file_names = ["smpl_1", "smpl_2"]
config = {}
config["parsing_modules"] = ["demux","ccs"]
config["assembly_long_mods"] = ["canu","flye"]
config["binning_mods"] = ["maxbin2","metabat2"]
config["ann_modules"] = ["prokka","KEGG_KM","BioV_transp"]
## ---


## find inputs & save path and filenames as table
def find_inputs(dir, filename_pattern, double_ext=False, uniq_PE=False):
    import pandas as pd
    import glob
    import os
    import re
    
    file_list = glob.glob(os.path.join(dir, "**/*"+filename_pattern), recursive=True)
    filename_list = [os.path.basename(i) for i in file_list]
    # filename_noext_list = [filename[:filename.index('.')] for filename in filename_list] # chop at the first dot position, i.e. no "dots" allowed in filename
    if double_ext: # e.g. '.tar.gz', '.fq.gz'
        # filename_noext_list = [re.sub('.tar.[a-z]+$', '', filename) for filename in filename_list]
        # ext_list = [re.sub('.*.tar.', 'tar.', filename) for filename in filename_list]
        filename_noext_list = [filename.split(".")[0] for filename in filename_list]
        ext_list = ['.'+'.'.join(filename.split(".")[1:]) for filename in filename_list]
    else:
        filename_noext_list = [os.path.splitext(filename)[0] for filename in filename_list]
        ext_list = [os.path.splitext(filename)[1] for filename in filename_list]
    path_list = [os.path.dirname(i) for i in file_list]
    
    df = pd.DataFrame ([path_list, filename_list, filename_noext_list, ext_list]).transpose()
    df.columns = ['dir', 'file', 'file_noext', 'ext']
    df = df.sort_values(by=['file'])

    if uniq_PE:
        df['uniq_PE'] = [re.sub('_R[12]$|_R[12]_.*', '', i) for i in df['file_noext']]
        df = df[~df.duplicated(subset=['uniq_PE'])]
        df = df.reset_index(drop=True)
    
    return df
## ---


## extract input names from barcode file and retrun a list
def barcode_names(barcode_file, pacbio=False):
    import re
    import pandas as pd

    fna = open(barcode_file, "rt")
    inputs = []
    for line in fna:
        if re.match('^>', line):
            inputs.append(re.sub('^>|\n', '', line))
    
    if pacbio: # demux tool lima modify file names as 'Sample1' -> 'Sample1--Sample1'
        inputs = pd.DataFrame([inputs, ["--"] * len(inputs), inputs]).transpose()
        inputs = inputs.apply(lambda x: "".join(x), axis=1)
        inputs = inputs.values.flatten().tolist()
    
    return inputs
## ---


## table (v2.0) with output file paths for parsing options activated
def result_table2(*arguments, sep = "_"):
    import itertools
    import pandas as pd
    
    in_lists = []
    for i in arguments:
        if isinstance(i, list):
            in_lists.append(i)
        else:
            in_lists.append([i])
    # in_lists = [*arguments]
    
    in_lists_comb = list(itertools.product(*in_lists))
    comb_df = pd.DataFrame(in_lists_comb)
    comb_df = comb_df.apply(lambda x: sep.join(x), axis=1)
    comb_list = comb_df.values.flatten().tolist()
    
    return comb_list
## ---


## table with output file paths for parsing options activated
def result_table(input_list, module_list, out_dir, pacbio=False):
    import os
    import pandas as pd
    
    out_list = []
    input_list2 = [os.path.basename(x) for x in input_list]
    
    if pacbio and "demux" in module_list:
        input_list2 = ["demux." + x for x in input_list2]
    
    module_list2 = module_list 
    for i in module_list:      
        
        #! parse_PacBio.smk options
        if i == "demux": extra_str = ".bam"
        elif i == "ccs": extra_str = ".fq.gz"
        elif i == "decontam": extra_str = ".fq.gz"
        
        #! MAG_discovery.smk options
        elif i == "assembly_canu": extra_str = "_canu"
        elif i == "assembly_flye": extra_str = "_flye"
        elif i == "binning_maxbin2": extra_str = "_maxbin2"
        elif i == "binning_metabat2": extra_str = "_metabat2"
        elif i == "dastool": extra_str = "_dastool"
        elif i == "drep": extra_str = "_drep"
        
        #! MAG_analysis.smk options
        elif i == "coverm": extra_str = "_coverm"
        elif i == "instrain": extra_str = "_profile"
        
        #! FunLuca.smk options
        elif i == "prokka": extra_str = ""
        elif i == "KEGG_KO": extra_str = ""
        elif i == "KEGG_KM":
            module_list2 = module_list2.remove("KEGG_KM")
            continue
        elif i == "KEGG_manual": 
            module_list2 = module_list2.remove("KEGG_manual")
            continue
        elif i == "antiSMASH": extra_str = ""
        elif i == "BioV_transp": extra_str = ""
        elif i == "blast_phytohormones": extra_str = ".tsv" 
        elif i == "blast_vibrioferrin": extra_str = ".tsv"
        elif i == "blast_DMSP": extra_str = ".tsv"
        elif i == "dbCAN_CAZy": extra_str = ""
        
        else: extra_str = ""

        
        i_dir = [os.path.join(out_dir, i)] * len(input_list2)
        i_df = pd.DataFrame([i_dir, [x+extra_str for x in input_list2]]).transpose()
                                                               
        i_out_list = i_df.apply(lambda x: os.path.join(*x), axis=1)
        out_list.append(i_out_list)

    df = pd.DataFrame([*out_list]).transpose()
    df.columns = module_list2
    
    return df
## ---


##
def path_from_match(gnm_table):
    import os 

    # gnm_table_f = gnm_table[gnm_table.dir.str.match(config["query_dir"].strip("/"))]
    gnm_table_f = gnm_table[gnm_table.dir == config["query_dir"].strip("/")]
    query_list = gnm_table_f.apply(lambda x : os.path.join(*x), axis=1).values.tolist()
    
    return query_list
## ---


## convert a text table to list
def file2list(input_file, col_id=False, sep="\t", header=None):
    import pandas as pd
    
    input_table = pd.read_csv(input_file, sep=sep, header=header)
    
    if col_id == "last":
        col_id = input_table.shape[1]-1

    if not col_id:
        list_out = input_table.values.flatten().tolist()
    elif isinstance(col_id, int) and col_id <= (input_table.shape[1]-1):
        list_out = input_table.iloc[:, col_id].values.flatten().tolist()
    elif isinstance(col_id, str) and col_id in input_table:
        list_out = input_table.loc[:, col_id].values.flatten().tolist()
    else:
        print("Something is wrong with the specified column. Use existing col names or ids\n") 
        exit()
        
    return list_out
## ---


## convert a whole pandas DataFrame to list
def table2list(input_table):
    
    input_table2 = input_table
    if ('file' in input_table2):
        input_table2 = input_table2.drop(['file'], axis=1)

    all_out = input_table.values.flatten().tolist()
    
    return all_out
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


## set time resource of a rule - VSC5
def time_MEGAHIT(wildcards, input):
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size < 4:
        time = ("00:%02d:00" % (59))
    elif in_size < 20:
        time = ("%02d:00:00" % (12))
    elif in_size < 45:
        time = ("%02d:00:00" % (24))
    else:
        time = ("%02d:00:00" % (32))

    return time

def qos_MEGAHIT(wildcards, input):
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size > 45:
        out_qos = "p71863_2048"
    else:
        out_qos = "p71863_0512"

    return out_qos

def part_MEGAHIT(wildcards, input):    
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size > 45:
        out_part = "zen3_2048"
    else:
        out_part = "zen3_0512"
    
    return out_part


# def node_HIPMER(wildcards, input):
#     in_size = input.size_mb/1024
#     in_size = max(1, int(in_size))

#     if in_size < 15:
#         node = 1
#     elif in_size < 40:
#         node = 2
#     elif in_size < 120:
#         node = 1
#     else:
#         node = 2

#     return node

def node_HIPMER(wildcards, input):
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size < 20:
        node = 1
    elif in_size < 40:
        node = 2
    elif in_size < 100:
        node = 3
    elif in_size < 140:
        node = 4
    elif in_size < 170:
        node = 5
    else:
        node = 6

    return node

def time_HIPMER(wildcards, input):
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size < 15:
        time = ("00:%02d:00" % (59))
    elif in_size < 40:
        time = ("00:%02d:00" % (59))
    elif in_size < 120:
        time = ("%02d:00:00" % (4))
    else:
        time = ("%02d:00:00" % (12))

    return time

def qos_HIPMER(wildcards, input):
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size < 15:
        qos = "p71863_0512"
    elif in_size < 40:
        qos = "p71863_0512"
    elif in_size < 120:
        qos = "p71863_2048"
    else:
        qos = "p71863_2048"

    return qos

def part_HIPMER(wildcards, input):
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))
    
    if in_size < 15:
        part = "zen3_0512"
    elif in_size < 40:
        part = "zen3_0512"
    elif in_size < 120:
        part = "zen3_2048"
    else:
        part = "zen3_2048"

    return part
## ---










### deprecated:


## extract match from list
def match_from_list(mylist, pattern):
    from itertools import compress
    
    bool_filter = [pattern in i for i in mylist]
    # x = list(compress(mylist, bool_filter))
    x = list(compress(mylist, bool_filter))
    
    return x
## ---


## set time resource of a rule - VSC4
def get_time_spades(wildcards, input):
    import math

    time_coef = 7
    run_time = input.size_mb/1024 * time_coef
    hours = max(2, math.ceil(run_time)) # to ensure min time of 2h
    if hours > 72:
        return -1 # kill job if it takes more than 3d
    time = ("%02d:00:00" % (hours))
    return time

def get_time_megahit(wildcards, input):
    import math
    
    time_coef = 0.5
    run_time = input.size_mb/1024 * time_coef
    hours = max(1, math.ceil(run_time)+1) # to ensure min time of 1h
    time = ("%02d:00:00" % (hours))
    return time

def get_qos_vsc4(wildcards, input): # VSC4
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size > 30:
        out_qos = "mem_0384"
    else:
        out_qos = "mem_0096"

    return out_qos
## ---


## print system date and time
def print_datetime():
    from datetime import datetime

    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    dt = "date and time ="+ dt_string
    
    return dt	
## ---


## create table with all genomes and relative paths
def make_gnm_table(config):
    import pandas as pd
    
    gnm_table = pd.DataFrame()
    
    if config['ref_gnm'] and config['ref_dir']:
        if isinstance(config['ref_gnm'], list):
            for gnm in config['ref_gnm']:
                gnm_table = pd.concat(
                    [gnm_table, pd.DataFrame(dict({'dir' : config['ref_dir'].strip("/"),
                                                   'file' : gnm}), index = [0])], 
                    ignore_index = True, axis = 0)
       
        else:
            gnm_table = pd.concat(
                [gnm_table, pd.DataFrame(dict({'dir' : config['ref_dir'].strip("/"),
                                               'file' : config['ref_gnm']}), index = [0])], 
                ignore_index = True, axis = 0)
    
    if config['query_gnm'] and config['query_dir']:
        if isinstance(config['query_gnm'], list):
            for gnm in config['query_gnm']:
                gnm_table = pd.concat(
                    [gnm_table, pd.DataFrame(dict({'dir' : config['query_dir'].strip("/"),
                                                   'file' : gnm}), index = [0])],
                ignore_index = True, axis = 0)
        
        else:
            gnm_table = pd.concat(
                [gnm_table, pd.DataFrame(dict({'dir' : config['query_dir'].strip("/"),
                                               'file' : config['query_gnm']}), index = [0])], 
                ignore_index = True, axis = 0)   
            
    return gnm_table
## ---



## create list of genomes with related path from make_gnm_table output
def path_all(gnm_table):    
    import os 
    
    all_list = gnm_table.apply(lambda x : os.path.join(*x), axis=1).values.tolist()
    
    return all_list

def path_ref(gnm_table):
    import os 

    # gnm_table_f = gnm_table[gnm_table.dir.str.match(config["ref_dir"].strip("/"))]
    gnm_table_f = gnm_table[gnm_table.dir == config["ref_dir"].strip("/")]
    ref_list = gnm_table_f.apply(lambda x : os.path.join(*x), axis=1).values.tolist()
    
    return ref_list

def path_query(gnm_table):
    import os 

    # gnm_table_f = gnm_table[gnm_table.dir.str.match(config["query_dir"].strip("/"))]
    gnm_table_f = gnm_table[gnm_table.dir == config["query_dir"].strip("/")]
    query_list = gnm_table_f.apply(lambda x : os.path.join(*x), axis=1).values.tolist()
    
    return query_list
## ---


## create file paths by rowise concatenating a table
def paths_from_rows(in_table, sel_file=False, sel_rows=False, sel_cols=False, extension=False):
    import os
    import pandas as pd
    
    in_table_filt = in_table
    
    if sel_file != False:
        if 'file' in in_table_filt:
            if sel_file in in_table_filt.file:
                in_table_filt = in_table_filt.loc[in_table_filt.file == sel_file, ]                            
            else:
                print("Required file is not in the table.\n") 
                exit() 
        else:
            print("You asked to select for a specific file but there is no 'file' column in the table.\n") 
            exit()             
        
    if sel_rows != False:
        if isinstance(sel_rows, list):
            sel_rows = [sel_rows]
        if all([isinstance(i, int) for i in sel_rows]):
            sel_rows = [x - 1 for x in sel_rows] # Assume people start counting from 1 and not 0 (as Python does).
            if max(sel_rows) > (in_table_filt.shape[0]): 
                print("Not enough rows for your selection.\n") 
                exit() 
            in_table_filt = in_table_filt.iloc[sel_rows, :]
        else: 
            print("Unsupported row selection. Function accepts only numerical indexes.\n") 
            exit() 
            
    if sel_cols != False:
        if sel_cols == "last_module":
            sel_cols = [(in_table_filt.shape[1]-1), 0]
            in_table_filt = in_table_filt.iloc[:, sel_cols]
        if sel_cols == "secondlast_module":
            sel_cols = [(in_table_filt.shape[1]-2), 0]
            in_table_filt = in_table_filt.iloc[:, sel_cols]
        else:
            if isinstance(sel_rows, list):
                sel_rows = [sel_rows]
            if all([isinstance(i, int) for i in sel_cols]):
                sel_cols = [x - 1 for x in sel_cols] # Assume people start counting from 1 and not 0 (as Python does).
                if max(sel_cols) > (in_table_filt.shape[1]): 
                    print("Not enough columns for your selection. Remember Python starts to count from 0.\n") 
                    exit() 
                in_table_filt = in_table_filt.iloc[:, sel_cols]
            elif all([isinstance(i, str) for i in sel_cols]):
                in_table_filt = in_table_filt.loc[:, sel_cols]
            else: 
                print("Unsupported column selection. Use either column names or column indexes.\n") 
                exit() 

    if extension != False:
        if 'file' in in_table_filt:
            pd.options.mode.chained_assignment = None  # default='warn'
            # in_table_filt.file = in_table_filt['file'].astype(str) + extension
            # in_table_filt.loc[:, 'file'] = in_table_filt.loc[:, 'file'].astype(str) + extension
            in_table_filt.file = [x + extension for x in in_table_filt.file]
        else:
            print("You asked to append an 'extension' but did NOT select the 'file' column.\n") 
            exit() 
    
    paths = in_table_filt.apply(lambda x: os.path.join(*x), axis=1)
    paths = paths.values.flatten().tolist()
    
    return paths
## ---