#!/usr/bin/env python
# -*- coding: utf-8

"""
name: ###
description: ###
author: Luca Zoccarato (luca.zoccarato@boku.ac.at)
dependencies:
    - parsed reads from 'parse_PacBio.smk'
TODO:
    - add inStrain profiling https://drep.readthedocs.io/en/latest/overview.html
    - implement sample metadata annotation in HM of merge_coverm
    - check log for commands using "2> {log.command}"
"""





# ##### Mapping genomes to fastq files

# if is_long_reads == "true":
#     db_ext = ".mmi"
# else:
#     db_ext = ".bt2"


# rule merge_plot_coverm:
#     """
#     Merge CoverM output of different samples.
#     Creates heatmaps showing patterns of MAGs TPM and completness across samples.
#     """
#     input:
#         coverm_out = coverm_list,
#         gtdbtk_dir = os.path.join(output_dir, 'gtdbtk'),
#     output:
#         merged_tpm = os.path.join(output_dir, 'coverm', "merged_TPM.tsv"),
#         merged_cover = os.path.join(output_dir, 'coverm', "merged_coverfrac.tsv"),
#         hm_tpm = os.path.join(output_dir, 'coverm', "hm_MAGs_TPM.html"),
#         hm_cover = os.path.join(output_dir, 'coverm', "hm_MAGs_coverfrac.html"),
#     params:
#         mapper_meth = config["coverm_methods"],
#         taxa_col = config["taxa_col"],
#         # smpl_col = config["smpl_col"],
#         dist_mt = "bray",
#         aggl_mt = "complete",
#     conda:
#         "../envs/R_plots.yaml"
#     script:
#         "../workflow/scripts/merge&plot_CoverM.R"


# if is_long_reads == "true":
#     rule mapping_MAGs_minimap2:
#         """
#         Map raw/cleaned reads agaist MAGs and assess MAG metrics with CoverM.
#         Choose with PacBio and Nanopore reads (uses minimap2).
#         """
#         input:
#             drep_mags_db = os.path.join(output_dir, 'mapping_mags', 'ALL_mags'+db_ext),
#             # parsed_reads = lambda wildcards: match_from_list(parsing_list, wildcards.sample),
#             parsed_reads = os.path.join(parsing_dir, 'decontam', '{sample}.fq.gz'),
#             MAG_dir = os.path.join(output_dir, 'drep_fasta'),
#         output:
#             map_sam = temp(os.path.join(output_dir, 'mapping_mags', '{sample}.sam')),
#             map_bam = os.path.join(output_dir, 'mapping_mags', '{sample}.bam'),
#             map_stats = os.path.join(output_dir, 'mapping_mags', '{sample}.stats'),
#             coverm_dir = directory(os.path.join(output_dir, 'coverm', '{sample}')),
#             coverm_out = os.path.join(output_dir, 'coverm', '{sample}', 'MAGs_mapping.txt'),
#         params:
#             mapping_algorithm = config["mapping_algorithm"],
#             file_ext = ".fa",
#             mapper_meth = config["coverm_methods"],
#             min_algn_len = 40,
#             min_identity = 97,
#             min_overlap = 80,
#             ncore = config["nCORE"], 
#         resources:
#             time = "12:00:00",
#             # partition = "mem_0384",
#             # qos = "mem_0384",
#         conda:
#             "../envs/mapping_reads.yaml"
#         log:
#             command = "_logs/mapping_MAGs_minimap2_{sample}.command"
#         shell:
#             '''
#             cmd="
#             minimap2
#             -ax {params.mapping_algorithm}
#             -o {output.map_sam}
#             -t {params.ncore}
#             {input.drep_mags_db}
#             {input.parsed_reads}
#             ;
#             samtools sort -@ {params.ncore} {output.map_sam} -o {output.map_bam}
#             ;
#             samtools index {output.map_bam}
#             ;
#             samtools idxstats {output.map_bam} > {output.map_stats}
#             ;
#             coverm 
#             genome
#             --bam-files {output.map_bam}
#             --genome-fasta-directory {input.MAG_dir}
#             --genome-fasta-extension {params.file_ext}
#             --methods {params.mapper_meth}
#             --min-covered-fraction 0 
#             --min-read-aligned-length {params.min_algn_len}
#             --min-read-percent-identity {params.min_identity}
#             --min-read-aligned-percent {params.min_overlap}
#             --output-file {output.coverm_out}
#             --threads {params.ncore}
#             ";
#             echo $cmd >> {log.command};
#             eval $cmd
#             '''


# else:
#     rule mapping_MAGs_bowtie2:
#         """
#         Map raw/cleaned reads agaist MAGs and assess MAG metrics with CoverM.
#         Choose with Illumina reads (uses bowtie2).
#         """
#         input:
#             drep_mags_db = os.path.join(output_dir, 'mapping_mags', 'ALL_mags'+db_ext),
#             parsed_R1 = os.path.join(parsing_dir, '3_remove_EUK', '{sample}_R1.fq.gz'),
#             parsed_R2 = os.path.join(parsing_dir, '3_remove_EUK', '{sample}_R2.fq.gz'),
#             MAG_dir = os.path.join(output_dir, 'drep_fasta'),
#         output:
#             map_sam = temp(os.path.join(output_dir, 'mapping_mags', '{sample}.sam')),
#             map_bam = os.path.join(output_dir, 'mapping_mags', '{sample}.bam'),
#             map_stats = os.path.join(output_dir, 'mapping_mags', '{sample}.stats'),
#             coverm_dir = directory(os.path.join(output_dir, 'coverm', '{sample}')),
#             coverm_out = os.path.join(output_dir, 'coverm', '{sample}', 'MAGs_mapping.txt'),
#         params:
#             file_ext = ".fa",
#             mapper_meth = config["coverm_methods"],
#             min_algn_len = 40,
#             min_identity = 97,
#             min_overlap = 80,
#             ncore = config["nCORE"], 
#         resources:
#             time = "12:00:00",
#             # partition = "mem_0384",
#             # qos = "mem_0384",
#         conda:
#             "../envs/mapping_reads.yaml"
#         log:
#             command = "_logs/mapping_MAGs_minimap2_{sample}.command"
#         shell:
#             '''
#             cmd="
#             bowtie2
#             -x {input.drep_mags_db}
#             -1 {input.parsed_R1}
#             -2 {input.parsed_R2}
#             -S {output.map_sam}
#             -p {params.ncore}
#             ;
#             samtools sort -@ {params.ncore} {output.map_sam} -o {output.map_bam}
#             ;
#             samtools index {output.map_bam}
#             ;
#             samtools idxstats {output.map_bam} > {output.map_stats}
#             ;
#             coverm 
#             genome
#             --bam-files {output.map_bam}
#             --genome-fasta-directory {input.MAG_dir}
#             --genome-fasta-extension {params.file_ext}
#             --methods {params.mapper_meth}
#             --min-covered-fraction 0
#             --min-read-aligned-length {params.min_algn_len}
#             --min-read-percent-identity {params.min_identity}
#             --min-read-aligned-percent {params.min_overlap}
#             --output-file {output.coverm_out}
#             --threads {params.ncore}
#             ";
#             echo $cmd >> {log.command};
#             eval $cmd
#             '''


# rule make_mags_db:
#     """
#     Prepare MAGs DB to screen reads.
#     """
#     input:
#         MAG_dir = os.path.join(output_dir, 'drep_fasta'),
#     output:
#         drep_mags_fa = os.path.join(output_dir, 'mapping_mags', 'ALL_mags.fa'),          
#         drep_mags_db = os.path.join(output_dir, 'mapping_mags', 'ALL_mags'+db_ext),
#     params:
#         is_long = is_long_reads,
#         file_ext = ".fa",
#         ncore = config["nCORE"],       
#     conda:
#         "../envs/mapping_reads.yaml"
#     log:
#         command = "_logs/make_mags_db.command"
#     shell:
#         '''
#         if [[ {params.is_long} = true ]]
#         then
#         cmd="
#         cat {input.MAG_dir}/*{params.file_ext} > {output.drep_mags_fa}
#         ;
#         minimap2 -d {output.drep_mags_db} {output.drep_mags_fa}
#         ";
#         else
#         cmd="
#         cat {input.MAG_dir}/*{params.file_ext} > {output.drep_mags_fa}
#         ;
#         bowtie2-build --threads {params.ncore} {output.drep_mags_fa} {output.drep_mags_db}
#         ;
#         touch {output.drep_mags_db}
#         ";
#         fi;
#         echo $cmd >> {log.command};
#         eval $cmd
#         '''





##### MAGs stats rules
rule run_gtdbtk:
    """
    Taxonomic classification of dereplicated MAGs with GTDB-Tk.
    """
    input:
        gnm_dir = genome_dir,
        GTDBTk_path = directory(config["GTDBTk_path"]),
    output:
        gtdbtk_dir = directory(os.path.join(output_dir, 'gtdbtk')),
    params:
        genome_ext = genome_ext,
        ncore = config["nCORE"],
        pplacer_cpus = 64, # issue with pplacer (it never ends when using more CPUs)
    conda:
        "../envs/GTDB_Tk.yaml"
    resources:
        time = "03:59:00",
        # partition = "mem_0384",
        # qos = "mem_0384",
    log:
        command = "_logs/run_gtdbtk.command"
    shell:
        '''
        GTDBTk_path=$(pwd)"/"{input.GTDBTk_path};
        GTDBTk_path=$(find $GTDBTk_path -type d -name "release*");
        export GTDBTK_DATA_PATH=$GTDBTk_path;
        cmd="
        gtdbtk
        classify_wf
        --genome_dir {input.gnm_dir}
        --out_dir {output.gtdbtk_dir}
        --extension {params.genome_ext}
        --cpus {params.ncore}
        --pplacer_cpus {params.pplacer_cpus}
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''
   

rule run_checkm:
    """
    Compute quality metrics on genomes.
    """
    input:
        gnm_dir = genome_dir,
    output:
        checkm_dir = directory(os.path.join(output_dir, 'checkm')),
    params:
        ext = re.sub("^\\.", "", genome_ext),
        ncore = 64, # issue with pplacer (it never ends when using more CPUs)
        # ncore = config["nCORE"],
    conda:
        "../envs/CheckM.yaml"
    # resources:
    #     time = "03:59:00",
    log:
        command = "_logs/run_checkm.command"
    shell:
        '''
        cmd="
        checkm
        lineage_wf
        -t {params.ncore}
        -x {params.ext}
        {input.gnm_dir}
        {output.checkm_dir}
        ;
        checkm 
        qa 
        -o 2 
        -t {params.ncore} 
        {output.checkm_dir}/lineage.ms 
        {output.checkm_dir} 
        > {output.checkm_dir}/results.tsv
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''