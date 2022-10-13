# Snakemake workflow for genetic trait annotation

This repo contains a Snakmeake implementation of the annotation pipeline published in  <https://www.nature.com/articles/s42003-022-03184-4> . It allows for and easy implementation and deployment of several annotation tools:

* prokka
* KEGG Orthology (including insight analyses of phytohormones production, DHPS and taurine utilization)
* BioVx (for membrane transporters)
* AntiSMASH (for secondary metabolites)
* Vibrioferrin biosynthesis and transport (blastp against UniProt)
* DMSP degradation pathways (blastp against UniProt)
* dbCAN (for automated carbohydrate-active enzyme annotation)

Thanks to the Snakemake architecture, it offers an easy scalability from local server to cluster computer. It's also easily expandable by new annotation tools.

Outputs that will be generated in the `results` and `plots` folders are:

* `MASTER_table.tsv` containing the comprehensive summary of all annotations. It's based on a GFF format in which each row represent a gene of a genomes in the input directory, while the columns holds the information of gene location and functional annotations.
*  An interactive HTML heatmap `hm_MASTERtraits_jacc.html` representing the annotated traits across genomes

#### ToDo list

* add full genome report with checkM, GTKDB
* add rules for mapping genomes to fastq files and make plots

## How to run

#### Set up the environment

You need to have Snakemake installed properly in your system. You can usually achieve this by running the following commands but refer the official documentation for any issue or doubts (<https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>). In the terminal type:

    wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh
    chmod +x Miniconda3-py38_4.12.0-Linux-x86_64.sh
    bash Miniconda3-py38_4.12.0-Linux-x86_64.sh
    # close and reopen the terminal
    conda update conda
    conda install mamba -n base -c conda-forge
    mamba create -c conda-forge -c bioconda -n snakemake snakemake

Then fetch the code repository from GitHub:

`git https://github.com/lucaz88/FunLuca.git`

#### Prepare input files

You can either test the WF with the two provided genomes in the folder `input_genomes`, or you can delete them and place your desired genomes in that folder.

Genomes are expected to be placed in that specific folder (i.e. `input_genomes`), have the extension `.fna` and not be compressed. However you can change these parameters by editing the variables `genome_dir`, `genome_ext` and `sfds` in the file `workflow/Run_FunLuca.smk`.

#### Edit the Snakemake config file

This files controls how Snakemake use the available resources to run the WF. Two files are already provided to run on local (`configs/snakemake/local/config.yaml`) and on a computing cluster (`configs/snakemake/sbatch/config.yaml`), however check the "Profiles" section in the Snakemake documentation (https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) to learn how to properly tune such file.

#### Run the WF

In the terminal, prompt the following command to check that everything works properly:

`snakemake -s workflow/Run_FunLuca.smk --profile configs/snakemake/local -pn`

and then run it (`-n` flag execute a dry run that just output the WF rules without running them)

`snakemake -s workflow/Run_FunLuca.smk --profile configs/snakemake/local -p`

If you wanna check the rule scheme that would be executed by Snakemake and save it as PDF run

`snakemake -s workflow/Run_FunLuca.smk --forceall --rulegraph | dot -Tpdf > dag.pdf`
