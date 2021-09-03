Metagenome Pipeline
-------------------------------------------------
This project's aim is focused on automating processes of a shotgun metagenomic analysis and learn the pipeline development tool Nextflow. The ultimate aim is to make the process more accessible to a wider audience. 

## Prerequisities
* Linux based system
* nextflow

## Install

#TODO update the main.nf call here.

```shell
$ git clone https://github.com/lorentzben/metagenome_pipeline.git
$ cd metagenome_pipeline
$ nextflow run main.nf 
```

## Setup diamond database
```shell
docker run lorentzb/diamond bash setup_diamond_db.sh

singularity run docker://lorentzb/diamond bash setup_diamond_db.sh


```