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


rule update_GTDBTk_db:
    output:
        GTDBTk_dir = directory(config["GTDBTk_dir"]),
        GTDBTk_token = config["GTDBTk_db"],
    conda:
        "../envs/GTDB_Tk.yaml",
    shell:
        '''
        cmd="
        test -d {output.GTDBTk_dir} || mkdir {output.GTDBTk_dir};
        cd {output.GTDBTk_dir};
        wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz;
        tar xvzf gtdbtk_v2_data.tar.gz;
        touch {output.GTDBTk_token};
        ";
        eval $cmd
        '''

        
rule update_dbCAN_DB:
    output:
        dbCAN_dir = directory(config["dbCAN_db"]),
    conda:
        "../envs/dbCAN.yaml",
    shell:
        '''
        cmd="
        test -d {output.dbCAN_dir} || mkdir {output.dbCAN_dir};
        cd {output.dbCAN_dir};
        wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.09242021.fa && diamond makedb --in CAZyDB.09242021.fa -d CAZy;
        wget https://bcb.unl.edu/dbCAN2/download/Databases/V10/dbCAN-HMMdb-V10.txt && mv dbCAN-HMMdb-V10.txt dbCAN.txt && hmmpress dbCAN.txt;
        wget http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb;
        wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm && hmmpress tf-1.hmm;
        wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm && hmmpress tf-2.hmm;
        wget http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm && hmmpress stp.hmm;
        ";
        eval $cmd
        '''


rule update_antismash_DB:
    output:
        fake_out = config["antim_db"],
    conda:
        "../envs/antiSMASH.yaml",
    shell:
        '''
        cmd="
        download-antismash-databases;
        touch {output.fake_out};
        ";
        eval $cmd
        '''


rule get_kofamscan_DB:
    output:
        kegg_db = directory(config["kegg_db"]),
        kegg_profiles = directory(os.path.join(config["kegg_db"], "profiles")),
        kegg_ko_list = os.path.join(config["kegg_db"], "ko_list"),
    log:
        command = "_logs/get_kofamscan_DB.command",
    shell:
        '''
        cmd="
        mkdir -p {output.kegg_db};
        wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz -O {output.kegg_db}/profiles.tar.gz;
        tar zxf {output.kegg_db}/profiles.tar.gz -C {output.kegg_db};
        rm {output.kegg_db}/profiles.tar.gz;
        wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz -O {output.kegg_db}/ko_list.gz;
        gzip -d {output.kegg_db}/ko_list.gz;
        touch {output.kegg_profiles} {output.kegg_ko_list};
        ";
        echo $cmd >> {log.command}; 
        eval $cmd
        '''