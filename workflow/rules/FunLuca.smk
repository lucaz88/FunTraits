#!/usr/bin/env python
# -*- coding: utf-8

"""
name: Annotation rules
author: Luca Zoccarato (luca.zoccarato@boku.ac.at)
TODO:
    - fix issue in "KM_reconstruction_wrapper", function "KMdiagram_fetcher", that says "Received HTTP/0.9 when not allowed" 
    - in antismash parsing, add gene_ann info and standardize last 2 cols of parsed table to be gene_ann & trait
    - in CAZy parsing, try to implement info on CAZy clusters (CGCs) annotation e.g. PULs
"""





##### summary table                

rule make_master_table:
    """
    Combine all annotation tables into a MASTER table.
    """
    input:
        ann_tables = expand(os.path.join(output_dir, "{module}_tab.tsv"), module = config["ann_modules"]),
    output:
        MASTER_table = os.path.join(output_dir, "MASTER_table.tsv")
    params:
        KM_long_names = config["KM_long_names"],
    conda:
        "../envs/R_data_parsing.yaml"
    script:
        "../scripts/make_master_table.R"
 




##### annotation rules   

rule parse_dbCAN:
    """
    Filter and merge CAZy annotations from dbCAN into a unique table.
    Edit R script to change the settings for hits filtering.
    """
    input:
        CAZy_folder = expand(os.path.join(output_dir, "dbCAN_CAZy", "{genome}"), 
                             genome=genome_names),
    output:
        CAZy_table = os.path.join(output_dir, "dbCAN_CAZy_tab.tsv")
    params:
        CAZy_subs = config["CAZy_subs"],
    conda:
        "../envs/R_data_parsing.yaml"
    script:
        "../scripts/parse_dbCAN.R"


rule run_dbCAN:
    """
    Run the automated Carbohydrate-active enzyme ANnotation.
    """
    input:
        i_gnm_dir = os.path.join(output_dir, "prokka", "{genome}"),
        ref_db = config["dbCAN_db"],
    output:
        i_cazy = directory(os.path.join(output_dir, "dbCAN_CAZy", "{genome}")),
        i_no_seq_gff = temp(os.path.join(output_dir, "dbCAN_CAZy", "{genome}", "no_seq.gff")),
    params:
        core = config["nCORE"],
    conda:
        "../envs/dbCAN.yaml",
    log:
        command = "_logs/run_dbCAN_{genome}.command",
    shell:
        '''
        cmd="
        sed '/##FASTA/Q' {input.i_gnm_dir}/*.gff > {output.i_no_seq_gff};
        run_dbcan
        {input.i_gnm_dir}/*.faa
        protein
        -c {output.i_no_seq_gff}
        --out_dir {output.i_cazy}
        --db_dir {input.ref_db}
        --dia_cpu {params.core}
        --hmm_cpu {params.core}
        --eCAMI_jobs {params.core};
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''


rule parse_DMSP:
    """
    Filter and merge DMSP annotations from blast search into a unique table.
    Edit R script to change the settings for hits filtering.
    """
    input:
        DMSP_files = expand(os.path.join(output_dir, "blast_DMSP", "{genome}.tsv"), 
                               genome=genome_names),
    output:
        blast_DMSP_tab = os.path.join(output_dir, "blast_DMSP_tab.tsv")
    conda:
        "../envs/R_data_parsing.yaml"
    script:
        "../scripts/parse_DMSP.R"
        
        
rule blast_DMSP:
    """
    Run blast against custom database of DMSP degradation genes.
    """
    input:
        i_gnm_dir = os.path.join(output_dir, "prokka", "{genome}"),
        ref_token = config["DMSP_db"]+"_token", 
    output:
        i_phytohormone = os.path.join(output_dir, "blast_DMSP", "{genome}.tsv"),
    params:
        core = config["nCORE"],
        ref_db = config["DMSP_db"],
    conda:
        "../envs/blast.yaml",
    log:
        command = "_logs/blast_DMSP_{genome}.command",
    shell:
        '''
        cmd="
        blastp
        -query {input.i_gnm_dir}/*.faa
        -db {params.ref_db}
        -out {output.i_phytohormone}
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'
        -evalue 1e-5
        -max_target_seqs 20
        -num_threads {params.core};
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''


rule makeblastdb_DMSP:
    """
    Build reference database for annotation of DMSP-degradation genes.
    """
    input:
        ref_fasta = config["DMSP_db"],
    output:
        ref_token = config["DMSP_db"]+"_token",
    params:
        core = config["nCORE"],
    conda:
        "../envs/blast.yaml",
    log:
        command = "_logs/makeblastdb_DMSP.command",
    shell:
        '''
        cmd="
        makeblastdb 
        -in {input.ref_fasta} 
        -dbtype prot;
        touch {output.ref_token}
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''


rule parse_vibrioferrin:
    """
    Filter and merge vibrioferrin annotations from blast search into a unique table.
    Edit R script to change the settings for hits filtering.
    """
    input:
        vibferr_files = expand(os.path.join(output_dir, "blast_vibrioferrin", "{genome}.tsv"), 
                               genome=genome_names),
    output:
        blast_vibrioferrin_tab = os.path.join(output_dir, "blast_vibrioferrin_tab.tsv")
    conda:
        "../envs/R_data_parsing.yaml"
    script:
        "../scripts/parse_vibrioferrin.R"
        
        
rule blast_vibrioferrin:
    """
    Run blast against custom database of vibrioferrin production & uptake genes.
    """
    input:
        i_gnm_dir = os.path.join(output_dir, "prokka", "{genome}"),
        ref_token = config["vibrioferrin_db"]+"_token", 
    output:
        i_vibrioferrin = os.path.join(output_dir, "blast_vibrioferrin", "{genome}.tsv"),
    params:
        core = config["nCORE"],
        ref_db = config["vibrioferrin_db"],
    conda:
        "../envs/blast.yaml",
    log:
        command = "_logs/blast_vibrioferrin_{genome}.command",
    shell:
        '''
        cmd="
        blastp
        -query {input.i_gnm_dir}/*.faa
        -db {params.ref_db}
        -out {output.i_vibrioferrin}
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'
        -evalue 1e-5
        -max_target_seqs 20
        -num_threads {params.core};
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''


rule makeblastdb_vibrioferrin:
    """
    Build reference database for annotation of vibrioferrin production and uptake genes.
    """
    input:
        ref_fasta = config["vibrioferrin_db"],
    output:
        ref_token = config["vibrioferrin_db"]+"_token",
    params:
        core = config["nCORE"],
    conda:
        "../envs/blast.yaml",
    log:
        command = "_logs/makeblastdb_vibrioferrin.command",
    shell:
        '''
        cmd="
        makeblastdb 
        -in {input.ref_fasta} 
        -dbtype prot;
        touch {output.ref_token}
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''


rule parse_phytohormones:
    """
    Filter and merge phytohormone annotations from blast search into a unique table.
    Edit R script to change the settings for hits filtering.
    """
    input:
        phyhorm_files = expand(os.path.join(output_dir, "blast_phytohormones", "{genome}.tsv"), 
                               genome=genome_names),
        suppl_tab = config["phyhorm_suppl_tab"],
    output:
        blast_phytohormones_tab = os.path.join(output_dir, "blast_phytohormones_tab.tsv")
    conda:
        "../envs/R_data_parsing.yaml"
    script:
        "../scripts/parse_phytohormones.R"
        
        
rule blast_phytohormones:
    """
    Run blast against custom database of phytohormone production genes.
    """
    input:
        i_gnm_dir = os.path.join(output_dir, "prokka", "{genome}"),
        ref_token = config["phytohormones_db"]+"_token", 
    output:
        i_phytohormone = os.path.join(output_dir, "blast_phytohormones", "{genome}.tsv"),
    params:
        core = config["nCORE"],
        ref_db = config["phytohormones_db"],
    conda:
        "../envs/blast.yaml",
    log:
        command = "_logs/blast_phytohormones_{genome}.command",
    shell:
        '''
        cmd="
        blastp
        -query {input.i_gnm_dir}/*.faa
        -db {params.ref_db}
        -out {output.i_phytohormone}
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'
        -evalue 1e-5
        -max_target_seqs 20
        -num_threads {params.core};
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''


rule makeblastdb_phytohormones:
    """
    Build reference database for annotation of phytohormone production genes.
    """
    input:
        ref_fasta = config["phytohormones_db"],
    output:
        ref_token = config["phytohormones_db"]+"_token",
    params:
        core = config["nCORE"],
    conda:
        "../envs/blast.yaml",
    log:
        command = "_logs/makeblastdb_phytohormones.command",
    shell:
        '''
        cmd="
        makeblastdb 
        -in {input.ref_fasta} 
        -dbtype prot;
        touch {output.ref_token}
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''


rule parse_gblast:
    """
    Filter and merge transporter annotations from gblast-BioV suite into a unique table.
    Filters for specific transporter sets are available in the R code.
    Edit R script to change the settings for hits filtering.
    """
    input:
        transp_folder = expand(os.path.join(output_dir, "BioV_transp", "{genome}"), 
                               genome=genome_names),
    output:
        BioV_transp_tab = os.path.join(output_dir, "BioV_transp_tab.tsv")
    params:
        transp_set = config["transp_set"]
    conda:
        "../envs/R_data_parsing.yaml"
    script:
        "../scripts/parse_gblast.R"


rule run_gblast:
    """
    Run gblast script of the BioV suit to annotate transporters using TCdb.
    """
    input:
        i_gnm_dir = os.path.join(output_dir, "prokka", "{genome}"),
    output:
        i_transp = directory(os.path.join(output_dir, "BioV_transp", "{genome}")),
    container:
        "docker://lucaz88/biovx"
    # conda:
    #     "../envs/BioVx.yaml",
    log:
        command = "_logs/run_gblast_{genome}.command",
    shell:
        '''
        HOME=$(pwd)"/databases/tcdb"
        cmd="
        gblast3.py
        -i {input.i_gnm_dir}/*.faa
        -o {output.i_transp};
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''
        # '''
        # cmd="
        # if [[ ! -L $CONDA_PREFIX/bin/hmmtop ]]; then ln -sf $HMMTOP_PATH/hmmtop $CONDA_PREFIX/bin; fi;
        # if [[ ! -L $CONDA_PREFIX/bin/gblast3.py ]]; then ln -sf $BIOVX_PATH/*.py $CONDA_PREFIX/bin; fi;
        # python2 $CONDA_PREFIX/bin/gblast3.py
        # -i {input.i_gnm_dir}/*.faa
        # -o {output.i_transp};
        # ";
        # echo $cmd >> {log.command};
        # eval $cmd
        # '''


rule parse_antismash:
    """
    Merge all secondary metabolite annotations from antiSMASH into a unique table.
    """
    input:
        html_folder = expand(os.path.join(output_dir, "antiSMASH", "{genome}"), 
                             genome=genome_names),
        gff_folder = expand(os.path.join(output_dir, "prokka", "{genome}"), 
                            genome=genome_names),
    output:
        antiSMASH_tab = os.path.join(output_dir, "antiSMASH_tab.tsv")
    conda:
        "../envs/R_data_parsing.yaml"
    script:
        "../scripts/parse_antismash.R"


rule run_antismash:
    """
    Run antiSMASH to annotate secondary metaolites.
    """
    input:
        i_gnm_dir = os.path.join(output_dir, "prokka", "{genome}"),
        antismash_db = config["antim_db"],
    output:
        i_KO = directory(os.path.join(output_dir, "antiSMASH", "{genome}")),
    params:
        core = config["nCORE"], 
    conda:
        "../envs/antiSMASH.yaml",
    log:
        command = "_logs/run_antismash_{genome}.command",
    shell:
        '''
        cmd="
        antismash
        -c {params.core}
        --taxon bacteria 
        --fullhmmer
        --clusterhmmer
        --tigrfam
        --smcog-trees
        --tta-threshold 0.65
        --cb-general
        --cb-subclusters
        --cb-knownclusters
        --asf
        --pfam2go
        --rre
        --cc-mibig
        --genefinding-tool prodigal-m
        --output-dir {output.i_KO}
        {input.i_gnm_dir}/*.gbk;
        #--cassis
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''


rule parse_manual_fromKOs:
    """
    Add manual trait annotations using parsed KO table.
    Currently implemented trait: DHPS and taurine utilization.
    """
    input:
        KEGG_KO_tab = os.path.join(output_dir, "KEGG_KO_tab.tsv"),
    output:
        KEGG_manual_tab = os.path.join(output_dir, "KEGG_manual_tab.tsv")
    conda:
        "../envs/R_data_parsing.yaml"
    script:
        "../scripts/parse_manual_fromKOs.R"


rule KM_reconstruction_wrapper:
    """
    Recombine KO annotations into a list of complete KEGG modules and add this info to the kofamscan table.
    """
    input:
        KEGG_KO_tab = os.path.join(output_dir, "KEGG_KO_tab.tsv"),
    output:
        KEGG_KM_tab = os.path.join(output_dir, "KEGG_KM_tab.tsv"),
    params:
        KM_str = config["KM_str"],
        ncore = config["nCORE"], 
    conda:
        "../envs/R_data_parsing.yaml"
    script:
        "../scripts/KM_reconstruction_wrapper.R"


rule parse_kofamscan:
    """
    Filter and merge KO annotations from kofamscan into a unique table.
    Only signif annotation are retained, see https://doi.org/10.1093/bioinformatics/btz859).
    """
    input:
        ko_files = expand(os.path.join(output_dir, "KEGG_KO", "{genome}_ko.txt"), 
                          genome=genome_names),
    output:
        KEGG_KO_tab = os.path.join(output_dir, "KEGG_KO_tab.tsv")
    conda:
        "../envs/R_data_parsing.yaml"
    script:
        "../scripts/parse_kofamscan.R"


rule run_kofamscan:
    """
    Run kofamscan to annotate KOs.
    """
    input:
        i_gnm_dir = os.path.join(output_dir, "prokka", "{genome}"),
        kegg_profiles = os.path.join(config["kegg_db"], "profiles"),
        kegg_ko_list = os.path.join(config["kegg_db"], "ko_list"),
    output:
        i_KO = os.path.join(output_dir, "KEGG_KO", "{genome}_ko.txt"),
        i_tmp = directory(os.path.join(output_dir, "KEGG_KO", "{genome}")),
    params:
        core = 100 #config["nCORE"], 
    conda:
        "../envs/kofamscan.yaml",
    log:
        command = "_logs/run_kofamscan_{genome}.command",
    shell:
        '''
        cmd="
        exec_annotation
        -o {output.i_KO}
        --profile={input.kegg_profiles}
        --ko-list={input.kegg_ko_list}
        --cpu={params.core}
        --tmp-dir={output.i_tmp}
        {input.i_gnm_dir}/*.faa;
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''


rule parse_prokka:
    """
    Merge all GFF files from prokka into a unique table.
    """
    input:
        gff_folder = expand(os.path.join(output_dir, "prokka", "{genome}"), 
                            genome=genome_names),
    output:
        GFF_table = os.path.join(output_dir, "prokka_tab.tsv")
    conda:
        "../envs/R_data_parsing.yaml"
    script:
        "../scripts/parse_prokka.R"
        
        
rule run_prokka:
    """
    Run prokka annotation on fasta files.
    """
    input:
        i_gnm = lambda wildcards: os.path.join(genome_dir, 
                                  genome_tab.file[genome_tab.file_noext == wildcards.genome].squeeze()),
    output:
        i_gnm_dir = directory(os.path.join(output_dir, "prokka", "{genome}")),
    params:
        core = config["nCORE"], 
    conda:
        "../envs/prokka.yaml",
    log:
        command = "_logs/run_prokka_{genome}.command",
    shell:
        '''
        cmd="
        prokka
        --outdir {output.i_gnm_dir}
        {input.i_gnm}
        --cpus {params.core}
        --compliant
        --rnammer
        --force
        #--norrna
        #--notrna;
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''