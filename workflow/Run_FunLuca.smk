#!/usr/bin/env python
# -*- coding: utf-8

"""
name: Trait-base annotation of prokaryotes
description: Snakemake pipeline for functional annotation of prokaryote genomes.
Based on the annotation pipeline published in https://doi.org/10.1038/s42003-022-03184-4
An annotation scaffold is built using prokka and then is enriched using:
KOfamscan, AntiSMASH, TCdb, CAZy, and specific dbs for phytohormones, vibrioferrin and DMSP
author: Luca Zoccarato (luca.zoccarato@boku.ac.at)
dependencies:
    - folder with genome fasta files
    - run_kofamscan, run_antismash need reference DBs which can be downloaded using the workflow "get_DBs.smk" (!needs FTP access)
TODO:
    - fix issue in "KM_reconstruction_wrapper", function "KMdiagram_fetcher", that says "Received HTTP/0.9 when not allowed" 
    - in antismash parsing, add gene_ann info and standardize last 2 cols of parsed table to be gene_ann & trait
    - in CAZy parsing, try to implement info on CAZy clusters (CGCs) annotation e.g. PULs
"""


#! load config files
configfile: "configs/FunLuca_config.yaml"



#! import python modules
import os
from scripts.myPy_FUN import find_inputs
from scripts.myPy_FUN import col2list
from scripts.myPy_FUN import comb_lists
from scripts.myPy_FUN import which_drep



#! setting vars
genome_dir = "input_genomes"
genome_ext = ".fna"
sample_tab = find_inputs(genome_dir, genome_ext, 
                      double_ext=True, uniq_PE=True)
sample_names = col2list(sample_tab, "uniq_PE")



#! genome analysis
coverm_list = expand(os.path.join(genome_anal_dir, 'coverm', '{sample}'),
                     sample=sample_names)
profile_list = expand(os.path.join(genome_anal_dir, 'instrain', '{sample}'),
                     sample=sample_names)
run_basic = [
    os.path.join(genome_anal_dir, 'checkm'),
    os.path.join(genome_anal_dir, 'gtdbtk'),
]
run_coverm = os.path.join(genome_anal_dir, 'coverm', "merged_TPM.tsv")



#! genome annotation
gnm_dir = os.path.join(genome_anal_dir, 'drep_fasta')
gtdbtk_dir = os.path.join(genome_anal_dir, 'gtdbtk')
trait_hm = os.path.join(genome_plot, 'hm_MASTERtraits_jacc.html')
run_FunLuca = trait_hm



#! declare outputs
rule all:
    input:
        run_FunLuca



#! load pipelines
include: "../rules/FunLuca.smk"
include: "../rules/get_reports.smk"



#! useful bash commands
# git https://github.com/lucaz88/FunLuca.git
# conda activate snakemake
# snakemake -s workflow/Run_FunLuca.smk --profile configs/snakemake/local -pn
# snakemake -s workflow/Run_FunLuca.smk --forceall --rulegraph | dot -Tpdf > dag.pdf