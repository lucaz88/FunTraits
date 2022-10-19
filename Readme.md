# Snakemake workflow for genetic trait annotation

This repo contains a Snakmeake implementation of the annotation pipeline published in  <https://www.nature.com/articles/s42003-022-03184-4> . It allows for and easy implementation and deployment of several **annotation tools**:

* prokka
* KEGG Orthology (including insight analyses of phytohormones production, DHPS and taurine utilization)
* BioVx (for membrane transporters)
* AntiSMASH (for secondary metabolites)
* Vibrioferrin biosynthesis and transport (blastp against UniProt)
* DMSP degradation pathways (blastp against UniProt)
* dbCAN (for automated carbohydrate-active enzyme annotation)

Thanks to the Snakemake architecture, it offers an easy scalability from local server to cluster computer. It's also easily expandable by new annotation tools.

The **generated outputs** are:

* `_results/MASTER_table.tsv` containing the comprehensive summary of all annotations. It's based on a GFF format in which each row represent a gene of a genomes in the input directory, while the columns holds the information of gene location and functional annotations.
* Genome stats about completness (using CheckM) and taxonomy (using GTDB-Tk) in the folders `_results/checkm` and `_results/gtdbtk`
* An interactive HTML heatmap `_plot/hm_MASTERtraits_jacc.html` representing the annotated traits across genomes
* Logs of all executed rules in the directory `_logs`

### TODO list

* add feauture for mapping genomes to fastq files and make plots (code ready, just need implementation)

## How to run

#### Set up the environment

You need to have Snakemake, Mamba (recomended, not essential) and Singularity installed properly in your system. You can usually achieve this by running the following commands but refer the official documentations for any issue or doubts (<https://snakemake.readthedocs.io/en/stable/getting_started/installation.html> and <https://docs.sylabs.io/guides/3.0/user-guide/installation.html#installation>). In the terminal type:

    wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh
    chmod +x Miniconda3-py38_4.12.0-Linux-x86_64.sh
    bash Miniconda3-py38_4.12.0-Linux-x86_64.sh
    # close and reopen the terminal
    conda update conda
    conda install mamba -n base -c conda-forge
    mamba create -c conda-forge -c bioconda -n snakemake snakemake
    conda activate snakemake
    conda install -c conda-forge singularity

Then fetch the code repository from GitHub:

    git clone https://github.com/lucaz88/FunLuca.git
    cd FunLuca
    ls

#### Prepare input files

You can either try the WF with the test genomes provided in the folder `input_genomes`, or you can delete them and place your desired genomes in that folder.

Genomes are expected to be placed in that specific folder (i.e. `input_genomes`), have the extension `.fna` and not be compressed. However you can change these parameters by editing the variables `genome_dir`, `genome_ext` and `sfds` in the file `workflow/Run_FunLuca.smk`.

#### Edit the Snakemake config file

This files controls how Snakemake use the available resources to run the WF. Two files are already provided to run on local (`configs/snakemake/local/config.yaml`) and on a computing cluster (`configs/snakemake/sbatch/config.yaml`), however check the "Profiles" section in the Snakemake documentation (https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) to learn how to properly tune such file.
Maybe something you wanna set are the values `jobs` and `cores`, which control the degree of parallelization (i.e. the number of jobs to execute in parallel and how many cores EACH job can take), and adjust them to the resources available in your system. E.g. if your computer has 32 cores you could run 6 jobs in parallel with 5 cores each.

#### Run the WF

In the terminal, activate Snakemake (if you haven't done already) with `conda activate snakemake` and prompt the following command to check that everything works properly:

`snakemake -s workflow/Run_FunLuca.smk --profile configs/snakemake/local -pn`

and then run it (`-n` flag execute a dry run that just output the WF rules without running them)

`snakemake -s workflow/Run_FunLuca.smk --profile configs/snakemake/local -p`

If you wanna check the rule scheme that would be executed by Snakemake and save it as PDF run

`snakemake -s workflow/Run_FunLuca.smk --forceall --rulegraph | dot -Tpdf > dag.pdf`

!!! the first execution will take quite some time as it has to fetches the DBs of KEGG, CAZy and (especially) GTDB-Tk
