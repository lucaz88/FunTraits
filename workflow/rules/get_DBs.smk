#!/usr/bin/env python
# -*- coding: utf-8

"""
name: Fetches DBs used for the trait-base annotation workflow.
author: Luca Zoccarato (luca.zoccarato@boku.ac.at)
dependencies:
    - reomve the token inside a DB folder to trigger the update
TODO:
    - implement token "trick" of GTDBTK for all DBs
    - implement 'create_cazy_db -h' in update_dbCAN_DB
"""


rule get_tcdb:
    """
    Get tcdb, database to gblast.
    Database is stored locally in a tmp folder as container env seems to mount paths outside pwd as read-only.
    '/app/test.faa' is already in the container.
    """
    output:
        TCDB_db = temp(directory("tmp_tcdb")),
        fake_out = temp(directory("tmp_tcdb/fake_out")),
    container:
        "docker://lucaz88/biovx"
    log:
        command = "_logs/get_tcdb.command",
    shell:
        '''
        export HOME={output.TCDB_db};
        cmd="
        mkdir -p {output.TCDB_db}
        ;
        gblast3.py
        -i /app/test.faa
        -o {output.fake_out};
        ";
        echo $cmd >> {log.command};
        eval $cmd
        '''


rule get_GTDBTk_DB:
    output:
        GTDBTk_path = directory(os.path.join(ann_dbs_path, config["GTDBTk_path"])),
    # conda:
    #     "../envs/GTDB_Tk.yaml",
    log:
        command = "_logs/get_GTDBTk_DB.command",
    shell:
        '''
        cmd="
        mkdir -p {output.GTDBTk_path}
        ;
        cd {output.GTDBTk_path}
        ;
        wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz
        ;
        tar xvzf gtdbtk_v2_data.tar.gz
        ;
        rm gtdbtk_v2_data.tar.gz
        ";
        eval $cmd
        '''

        
rule update_dbCAN_DB:
    output:
        dbCAN_dir = directory(os.path.join(ann_dbs_path, config["dbCAN_db"])),
    conda:
        "../envs/dbCAN.yaml",
    log:
        command = "_logs/update_dbCAN_DB.command",
    shell:
        '''
        cmd="
        mkdir -p {output.dbCAN_dir}
        ;
        cd {output.dbCAN_dir}
        ;
        wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.09242021.fa && diamond makedb --in CAZyDB.09242021.fa -d CAZy
        ;
        wget https://bcb.unl.edu/dbCAN2/download/Databases/V10/dbCAN-HMMdb-V10.txt && mv dbCAN-HMMdb-V10.txt dbCAN.txt && hmmpress dbCAN.txt
        ;
        wget http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb
        ;
        wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm && hmmpress tf-1.hmm
        ;
        wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm && hmmpress tf-2.hmm
        ;
        wget http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm && hmmpress stp.hmm
        ;
        touch $(date +'%Y.%m.%d')
        ";
        eval $cmd
        '''


rule get_antismash_DB:
    output:
        fake_out = os.path.join(ann_dbs_path, config["antim_db"]),
    conda:
        "../envs/antiSMASH.yaml",
    log:
        command = "_logs/get_antismash_DB.command",
    shell:
        '''
        cmd="
        download-antismash-databases
        ;
        touch {output.fake_out}
        ";
        eval $cmd
        '''


rule get_kofamscan_DB:
    output:
        kegg_db = directory(os.path.join(ann_dbs_path, config["kegg_db"])),
        kegg_profiles = directory(os.path.join(ann_dbs_path, config["kegg_db"], "profiles")),
        kegg_ko_list = os.path.join(ann_dbs_path, config["kegg_db"], "ko_list"),
    log:
        command = "_logs/get_kofamscan_DB.command",
    shell:
        '''
        cmd="
        mkdir -p {output.kegg_db}
        ;
        wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz -O {output.kegg_db}/profiles.tar.gz
        ;
        tar zxf {output.kegg_db}/profiles.tar.gz -C {output.kegg_db}
        ;
        rm {output.kegg_db}/profiles.tar.gz
        ;
        wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz -O {output.kegg_db}/ko_list.gz
        ;
        gzip -d {output.kegg_db}/ko_list.gz
        ;
        touch {output.kegg_profiles} {output.kegg_ko_list}
        ;
        touch $(date +'%Y.%m.%d')
        ";
        echo $cmd >> {log.command}; 
        eval $cmd
        '''