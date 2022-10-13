# Snakemake workflow for genetic trait annotation

This repo contains a Snakmeake implementation of the annotation pipeline published in  <https://www.nature.com/articles/s42003-022-03184-4> . It allows for and easy implementation and deployment of several annotation tools:

* prokka
* KEGG Orthology (including insight analyses of phytohormones production, DHPS and taurine utilization)
* BioVx (for membrane transporters)
* AntiSMASH (for secondary metabolites)
* Vibrioferrin biosynthesis and transport (blastp against UniProt)
* DMSP degradation pathways (blastp against UniProt)
* 

Thanks to the Snakemake architecture, it offers an easy scalability from local server to cluster computer. It's also easily expandable by new annotation tools.

Provided outputs:

* `MASTER_table.tsv` containing the comprehensive summary of all annotations. It's based on a GFF format in which each row represent a gene of a genomes in the input directory, while the columns holds the information of gene location and functional annotations.
*  An interactive HTML heatmap `hm_MASTERtraits_jacc.html` representing the annotated traits across genomes
* !!!ToDo add full genome report with checkM, GTKDB

## How to run:

#### Set up the environment

You need to have Snakemake installed properly in your system (<https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>).

Open a terminal and fetch the code from GitHub:

`git https://github.com/lucaz88/FunLuca.git`

#### Prepare input files

You can either test the WF with the two provided genomes in the folder `input_genomes`, or you can delete them and place your desired genomes in that folder.

Genomes are expected to be placed in that specific folder (i.e. `input_genomes`) and have the extension `.fna`. However you can change these two parameters by editing the variables `genome_dir` and `genome_ext` in the file `workflow/Run_FunLuca.smk.`

#### Edit the Snakemake config file

\###

#### Run the WF

In the terminal, prompt the following command to check that everything works properly:

`snakemake -s workflow/Run_FunLuca.smk --profile configs/snakemake/local -pn`

and then run it (`-n` flag execute a dry run that just output the WF rules without running them)

`snakemake -s workflow/Run_FunLuca.smk --profile configs/snakemake/local -p`

If you wanna check the rule scheme that would be executed by Snakemake and save it as PDF run

`snakemake -s workflow/Run_FunLuca.smk --forceall --rulegraph | dot -Tpdf > dag.pdf`
