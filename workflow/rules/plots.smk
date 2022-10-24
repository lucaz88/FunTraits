#!/usr/bin/env python
# -*- coding: utf-8

"""
name: Plot rules for FunLuca
author: Luca Zoccarato (luca.zoccarato@boku.ac.at)
TODO:
    -
"""


rule heatmap_master_table:
    """
    Generate heat-map with traits in MASTER table.
    """
    input:
        MASTER_table = os.path.join(output_dir, "MASTER_table.tsv"),
        gtdbtk_dir = os.path.join(output_dir, 'gtdbtk'),
    output:
        outfile = os.path.join(plot_dir, 'hm_MASTERtraits_jacc.html'),
    params:
        ann_modules = config["ann_modules"],
        min_trait_occur = config["min_trait_occur"],
        dist_mt = config["dist_mt"],
        aggl_mt = config["aggl_mt"],
        taxa_col = config["taxa_col"],
        # smpl_col = config["smpl_col"],
    conda:
        "../envs/R_plots.yaml"
    script:
        "../workflow/scripts/heatmap_master_table.R"


# rule plot_FunLuca_MAGs:
#     """
#     Generate heat-map with traits mapped across MAGs.
#     """
#     input:
#         merged_tpm = os.path.join(analysis_dir, 'coverm', "merged_TPM.tsv"),
#         MASTER_table = os.path.join(annotation_dir, "MASTER_table.tsv"),
#     output:
#         out_dir = os.path.join(plot_dir, 'hm_MASTERtraits_mapped.html'),
#     params:
#         ann_modules = config["ann_modules"],
#         taxa_col = config["taxa_col"],
#         # smpl_col = config["smpl_col"],
#         min_trait_occur = config["min_trait_occur"],
#         min_tpm = config["min_tpm"],
#         dist_mt = config["dist_mt"],
#         aggl_mt = config["aggl_mt"],
#     conda:
#         "../envs/R_plots.yaml"
#     script:
#         "../workflow/scripts/plot_FunLuca_MAGs.R"