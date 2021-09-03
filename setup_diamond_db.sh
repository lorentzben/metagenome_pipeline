#!/usr/bin/env bash

mkdir ncbi_nr
cd ncbi_nr
mkdir logs

#first step is to download all of the nr files from ncbi

update_blastdb --decompress nr [*] > logs/blast_db_download.log

#prep ncbi db for use with diamond prepdb command
diamond prepdb --db nr -o diamond_db > logs/diamond_db_prep.log

#check to make sure that diamond is functional
diamond test > logs/diamond_test.log

#check to make sure the database is good to go
diamond dbinfo --db nr > logs/diamond_db_check.info