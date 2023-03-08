# Snakemake workflow for genetic trait annotation

This repo contains a Snakmeake implementation of the annotation pipeline published in  <https://www.nature.com/articles/s42003-022-03184-4> . The workflow (WF) will execute several **annotation tools**:

* prokka
* KEGG Orthology (including insight analyses of phytohormones production, DHPS and taurine utilization)
* BioVx (for membrane transporters)
* AntiSMASH (for secondary metabolites)
* Vibrioferrin biosynthesis and transport (blastp against UniProt)
* DMSP degradation pathways (blastp against UniProt)
* dbCAN (for automated carbohydrate-active enzyme annotation)

Thanks to the Snakemake architecture, the WF offers an easy scalability from local server to cluster computer. It can also be expanded to include new annotation tools.

The **generated outputs** are:

* `_results/MASTER_table.tsv` containing the comprehensive summary of all annotations. It's based on a GFF format in which each row represent a gene of a genomes in the input folder, while the columns holds the information of gene location and functional annotations.
* Genome stats about completness (using CheckM) and taxonomy (using GTDB-Tk) in the folders `_results/checkm` and `_results/gtdbtk`
* An interactive HTML heatmap `_plot/hm_MASTERtraits_jacc.html` representing the annotated traits across genomes. The file opens correclty with Google Chrome, but has issues with Firefox.
* Logs of all executed rules in the folder `_logs`
* All conda environments and docker used in the WF will be stored in the folder `_tools`
* Fetched databases (i.e. GTDBTk, KEGG pfam and dbCAN) will be stored in the folder `databases`. [!only if you know what you are doing: you can edit the parameter `ann_dbs_path` in the file `workflow/Run_FunLuca.smk` to point to an alternative path where you already have the correct databases saved.]

### TODO list

* add feauture for mapping genomes to fastq files and make plots (code ready, just need implementation)

## How to run

#### System requirements

In order to run smoothly it requeires a system with at least 55Gb of RAM as the WF runs GTDB-Tk to assign genomes' taxonomy.
Nodes/server used to run the WF needs to have FTP access (required to fetch annotation databases).

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

You can try the WF with the test genomes provided in the folder `input_genomes`.
If you wanna run the WF on your desired genomes you can either:

* copy their fasta files into the folder `input_genomes` and delete the test genomes. Fasta are expected to have the extension `.fna` and not be compressed.
* in the file `workflow/Run_FunLuca.smk`, edit the parameters `genome_dir` to specify a custom folder for the fasta files, `genome_ext` for alternative extension (e.g. `.fa`) and/or set `genome_compressed` to `True` if fasta are compresed

#### Edit the Snakemake config file

This files controls how Snakemake use the available resources to run the WF. Two files are already provided to run on local (`configs/snakemake/local/config.yaml`) and on a computing cluster (`configs/snakemake/sbatch/config.yaml`), however check the "Profiles" section in the Snakemake documentation (https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) to learn how to properly tune such file.
Maybe something you wanna set are the values `jobs` and `cores`, which control the degree of parallelization (i.e. the number of jobs to execute in parallel and how many cores EACH job can take), and adjust them to the resources available in your system. E.g. if your computer has 32 cores you could run 6 jobs in parallel with 5 cores each.

#### Run the WF

In the terminal, activate Snakemake and run the WF by selecting the relevant config file

* for local server execution:

    ```bash
    conda activate snakemake
    snakemake -s workflow/Run_FunLuca.smk --profile configs/snakemake/local -p
    ```

* for computer cluster execution (! remember to adjust the config file to your cluster environment):

    ```bash
    conda activate snakemake
    snakemake -s workflow/Run_FunLuca.smk --profile configs/snakemake/sbatch -p
    ```

If you just wanna check which rules will be exectuted in the WF without actually run it, you can perform a dry run by setting the `-n` flag, e.g.:
`snakemake -s workflow/Run_FunLuca.smk --profile configs/snakemake/local -pn`

If you wanna check the rule scheme of the WF and save it as PDF run

`snakemake -s workflow/Run_FunLuca.smk --forceall --rulegraph | dot -Tpdf > dag.pdf`

!!! Be aware that the first execution could take quite some time as snakemake will need to install all the required conda envromnents and singularity, as well as to fetche the DBs (especially KEGG ~2 Gb and GTDB-Tk ~66 Gb).
