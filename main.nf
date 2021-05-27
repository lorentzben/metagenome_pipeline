#!/usr/bin/env nextflow

def helpMessage(){
	log.info"""
	Usage:

    This script should be called as follows:
    srun --pty  -p inter_p  --mem=40G --nodes=2 --ntasks-per-node=8 --time=12:00:00 --job-name=qlogin /bin/bash -l
    TODO update this command 
    nextflow run main.nf 

    Main arguments:
        --input [path/to/folder]      Folder containing reads from this analysis
        --outdir [file]               The output directory where the results will be saved 

    """.stripIndent()
}


if(params.input){
    Channel
        .fromPath(params.input)
        .ifEmpty {exit 1, log.info "Cannot find path file ${input}"}
        .set{ ch_seqs }
}

//show help message 
if (params.help){
    helpMessage()
    exit 0 
}
