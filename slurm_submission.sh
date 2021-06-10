#!/bin/bash

#SBATCH --partition=batch
#SBATCH --job-name=Metagenomic_Pipeline
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
#SBATCH --mem=64gb

#Replace this with your UGA email to get notified on completion
#SBATCH --mail-user="bjl34716@uga.edu"
#SBATCH --mail-type=BEGIN,END,FAIL

module load Nextflow/20.04.1

nextflow run main.nf --input seqs --mapping mapping.tsv --outdir paired -resume