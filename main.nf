#!/usr/bin/env nextflow

def helpMessage(){
	log.info"""
	Usage:

    This script should be called as follows:
    srun --pty  -p inter_p  --mem=40G --nodes=2 --ntasks-per-node=8 --time=12:00:00 --job-name=qlogin /bin/bash -l
    TODO update this command 
    nextflow run main.nf --input seqs --mapping mapping.tsv --outdir paired 
    This program expects paired end reads, a single end read pipeline can be developed if required. 

    Main arguments:
        --input [path/to/folder]      Folder containing reads from this analysis
        --mapping [file]              Mapping file with, sample-id \t forward-read \t reverse-read
        --outdir [file]               The output directory where the results will be saved 

    """.stripIndent()
}


if(params.input){
    Channel
        .fromPath(params.input)
        .ifEmpty {exit 1, log.info "Cannot find path file ${input}"}
        .set{ ch_seqs }
}

if(params.mapping){
    Channel
        .fromPath(params.mapping)
        .ifEmpty {exit 1, log.info "Cannot find path file ${mapping}"}
        .into{ ch_mapping_file ; ch_mapping_file_sam ; ch_mapping_file_assembly ; ch_mapping_coverage ; ch_mapping_remapping ; ch_mapping_file_assembly_two ; ch_mapping_file_screen_and_combo; ch_mapping_file_find_orf ; ch_mapping_file_screen_orf; ch_mapping_file_cluster; ch_mapping_gene_lib; ch_mapping_supported_genes ; ch_mapping_derep_gene_lib ; ch_mapping_gene_bam_to_fasta; ch_mapping_file_rename; ch_mapping_lib_stat}
}

Channel
    .fromPath("${baseDir}/Gallus_gallus")
    .ifEmpty{ exit 1, log.info "Cannot find ref database, please collect: https://support.illumina.com/sequencing/sequencing_software/igenome.html"}
    .set{ ch_chicken_ref }

Channel 
    .fromPath("${baseDir}/screen_and_combine.sh")
    .ifEmpty{ exit 1, log.info "Cannot find screen and combine script" }
    .set{ ch_screening_script }

Channel
    .fromPath("${baseDir}/screen_orf_over_100.sh")
    .ifEmpty{ exit 1, log.info "Cannot find orf screen script"}
    .set { ch_orf_screen_script }

Channel
    .fromPath("${baseDir}/rename_post_orf.py")
    .ifEmpty{ exit 1, log.info "Cannot find rename script"}
    .set{ ch_rename_script }

Channel
    .fromPath("${baseDir}/de_duplicate.sh")
    .ifEmpty{ exit 1, log.info "Cannot find script "}
    .set{ ch_de_dup_script }

//show help message 
if (params.help){
    helpMessage()
    exit 0 
}

process QualityTrim{

    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    container "docker://lorentzb/trim-galore"

    input:
    path seq_dir from ch_seqs

    output:
    file "done.txt" into ch_done
    path "clean_reads" into ch_cleaned
    path "clean_reads" into ch_cleaned_contam
    

    shell:
    '''
    #!/usr/bin/env bash

    eval "$(conda shell.bash hook)"
    conda activate metagenome  

    trim_galore -q 20 --fastqc --paired --cores 4 -o clean_reads !{seq_dir}/*.fastq

    echo "done" > done.txt

    '''
}

process multiQC{
    publishDir "${params.outdir}/multi_qc", mode: 'copy'

    container "docker://lorentzb/multi_qc"

    input:
    path "clean_reads" from ch_cleaned

    output:
    file "multiqc_report.html" into ch_multi_qc_report

    script:
    """
    #!/usr/bin/env bash

    multiqc clean_reads/    
    
    """
    
}

process CheckForContamination{
    publishDir "${params.outdir}/decontaminate", mode: 'copy'

    container "docker://lorentzb/bowtie"

    input:
    file mapping from ch_mapping_file
    path "clean_reads" from ch_cleaned_contam
    path "Gallus_gallus" from ch_chicken_ref

    output:
    path "decontam" into ch_decontam

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    subprocess.run(['mkdir decontam'],shell=True)
    import pandas as pd
    samples = pd.read_csv("${mapping}",sep='\t')

    for index, row in samples.iterrows():
        forward = row['forward-read']
        reverse = row['reverse-read']
        stub = row['sequence-id']
        command = "bowtie2 -p 4 -x Gallus_gallus/genome -1 clean_reads/" +forward[:-6]+"_val_1.fq -2 clean_reads/"+reverse[:-6]+ "_val_2.fq -S decontam/"+stub+".sam"
        result = subprocess.run([command], shell=True)
    
  
    """
}

process ConvertSamtoBamandModify{
    publishDir "${params.outdir}/decontaminate", mode: 'copy'

    container "docker://lorentzb/samtools"

    input:
    file mapping from ch_mapping_file_sam
    path "decontam" from ch_decontam
    
    output:
    path "no_host" into ch_no_host_seqs
    path "no_host" into ch_no_host_mapping
    path "no_host" into ch_no_host_gene_lib
    path "no_host" into ch_no_host_supported_genes

    script:
    """
    #!/usr/bin/env python3
    import glob
    import subprocess

    sample_list = glob.glob('./decontam/*.sam')
    subprocess.run(['mkdir no_host'],shell=True)

    for sample in sample_list:
        sample = sample.split('/')[2]
        #convert sam to bam 
        sam_conv = 'samtools view -bS decontam/'+sample+' > '+sample[:-4]+'.bam'
        subprocess.run([sam_conv], shell=True)

        #filter the unmapped reads
        filter_command = 'samtools view -b -f 12 -F 256 '+sample[:-4]+'.bam > '+sample[:-4]+'_bothUnmapped.bam'
        subprocess.run([filter_command],shell=True)

        #split reads into fastq files
        sort_command = 'samtools sort -n -m 5G -@ 2 '+sample[:-4]+'_bothUnmapped.bam -o '+sample[:-4]+'_sorted.bam'
        subprocess.run([sort_command], shell=True)
        split_command = 'samtools fastq -@ 8 '+sample[:-4]+'_sorted.bam -1 no_host/'+sample[:-4]+'_R1.fastq.gz -2 no_host/'+sample[:-4]+'_R2.fastq.gz -0 /dev/null -s /dev/null -n'
        subprocess.run([split_command],shell=True)
    """
}

process RoundOneAssemble{

    publishDir "${params.outdir}/initial_assembly", mode: 'copy'

    container "docker://lorentzb/megahit"

    input:
    file mapping from ch_mapping_file_assembly
    path "no_host" from ch_no_host_seqs
    
    output:
    path "first_contigs" into ch_first_contigs
    path "first_contigs" into ch_remapping
    path "first_contigs" into ch_remapping_pull
    path "first_contigs" into ch_first_assembly
     
    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd

    samples = pd.read_csv("${mapping}",sep='\t')
    subprocess.run(['mkdir first_contigs'], shell=True)

    for index, row in samples.iterrows():
        stub = row['sequence-id']
    
        megahit_command = "megahit -1 no_host/"+stub+"_R1.fastq.gz -2 no_host/"+stub+"_R2.fastq.gz --presets meta-large -o first_contigs/"+stub+"_assembly"
        subprocess.run([megahit_command],shell=True)
    """
}

process FindUnmappedReads{
    publishDir "${params.outdir}/unmapped_reads", mode: 'copy'

    container "docker://lorentzb/bowtie"

    input:
    file mapping from ch_mapping_coverage
    path "first_contigs" from ch_remapping
    path "no_host" from ch_no_host_mapping
    
    output:
    path "sams" into ch_sams_to_be_separated
     

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd

    samples = pd.read_csv("${mapping}",sep='\t')
    subprocess.run(['mkdir sams'], shell=True)


    for index, row in samples.iterrows():
        stub = row['sequence-id']

        build_index = "bowtie2-build first_contigs/"+stub+"_assembly/final.contigs.fa final.contigs"
        subprocess.run([build_index], shell=True)

        map_reads = "bowtie2 -x final.contigs -1 no_host/"+stub+"_R1.fastq.gz -2 no_host/"+stub+"_R2.fastq.gz -S sams/"+stub+".sam" 
        subprocess.run([map_reads], shell=True)
    """
}

process PullUnmappedOut{

    publishDir "${params.outdir}/unmapped_reads", mode: 'copy'

    container "docker://lorentzb/samtools"

    input:
    path "first_contigs" from ch_remapping_pull
    file mapping from ch_mapping_remapping
    path "sams" from ch_sams_to_be_separated
    
    output:
    path "round_two_reads" into ch_reads_to_be_reassembled

    script:
    """
    #!/usr/bin/env python3
    
    import subprocess
    import pandas as pd

    samples = pd.read_csv("${mapping}",sep='\t')
    
    subprocess.run(['mkdir round_two_reads'],shell=True)

    for index, row in samples.iterrows():
        stub = row['sequence-id']

        #convert sam to bam 
        sam_conv = 'samtools view -bS sams/'+stub+'.sam > '+stub+'.bam'
        subprocess.run([sam_conv], shell=True)

        #filter the unmapped reads
        filter_command = 'samtools view -b -f 12 '+stub+'.bam > '+stub+'_bothUnmapped.bam'
        subprocess.run([filter_command],shell=True)

        #split reads into fastq files
        sort_command = 'samtools sort -n -m 5G -@ 2 '+stub+'_bothUnmapped.bam -o '+stub+'_sorted.bam'
        subprocess.run([sort_command], shell=True)

        split_command = 'samtools fastq -@ 8 '+stub+'_sorted.bam -1 round_two_reads/'+stub+'_R1.fastq.gz -2 round_two_reads/'+stub+'_R2.fastq.gz -0 /dev/null -s /dev/null -n'
        subprocess.run([split_command],shell=True)
    """
}

process RoundTwoAssembly{

    publishDir "${params.outdir}/second_assembly", mode: 'copy'

    container "docker://lorentzb/megahit"

    input:
    file mapping from ch_mapping_file_assembly_two
    path "round_two_reads" from ch_reads_to_be_reassembled
    
    output:
    path "second_contigs" into ch_second_contigs
    path "second_contigs" into ch_second_contigs_screen
    
    
    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd

    samples = pd.read_csv('${mapping}',sep='\t')
    subprocess.run(['mkdir second_contigs'], shell=True)

    for index, row in samples.iterrows():
        stub = row['sequence-id']
    
        megahit_command = "megahit -1 round_two_reads/"+stub+"_R1.fastq.gz -2 round_two_reads/"+stub+"_R2.fastq.gz --presets meta-large -o second_contigs/"+stub+"_assembly"
        subprocess.run([megahit_command],shell=True)
    """
}

process ScreenAndCombine{

    publishDir "${params.outdir}", mode: 'copy'

    container "docker://lorentzb/bioawk"

    input:
    file mapping from ch_mapping_file_screen_and_combo
    path "first_contigs" from ch_first_assembly
    path "second_contigs" from ch_second_contigs_screen
    file "screen_and_combine.sh" from ch_screening_script
    
    output:
    path "final_contigs" into ch_final_contigs_eval
    path "final_contigs" into ch_final_contigs_orfs
    
    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd

    samples = pd.read_csv('${mapping}',sep='\t')
    subprocess.run(['mkdir final_contigs'], shell=True)

    for index, row in samples.iterrows():
        stub = row['sequence-id']

        screen_and_comb_command = "echo " +stub+" | bash screen_and_combine.sh"
        subprocess.run([screen_and_comb_command],shell=True)
        

    """
}

process EvaluateAssembiles{

    publishDir "${params.outdir}", mode: 'copy'

    container "docker://lorentzb/assembly_stats"

    input:
    path "final_contigs" from ch_final_contigs_eval
    
    output:
    file "final_contigs_consolidated.csv" into ch_contig_stats
    
    
    shell:
    '''
    #!/usr/bin/env bash
    cd final_contigs
    
    for i in *.fasta; do
        python3 /opt/fasta_metadata_parser-0.0.16/fasta_meta_data_parser.py "$i"
    done

    cd ..

    python3 /opt/fasta_metadata_parser-0.0.16/consolidate_assembly_records.py final_contigs

    '''
}

process FindORF{

    publishDir "${params.outdir}/prodigal", mode: 'copy'

    container "docker://lorentzb/prodigal"

    input:
    file mapping from ch_mapping_file_find_orf
    path "final_contigs" from ch_final_contigs_orfs
    
    output:
    path "prodigal_inital" into ch_inital_orfs
    
    
    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd

    samples = pd.read_csv('${mapping}',sep='\t')
    subprocess.run(['mkdir prodigal_inital'], shell=True)

    for index, row in samples.iterrows():
        stub = row['sequence-id']

        prod_command = "prodigal.linux -i final_contigs/"+stub+"_final.fasta -o prodigal_inital/"+stub+".gbk  -d prodigal_inital/"+stub+".fna -p meta"
        subprocess.run([prod_command], shell=True)
    """
}

process RenameORFHeaders{

     publishDir "${params.outdir}/renamedFasta", mode: 'copy'

    container "docker://lorentzb/rename"

    input:
    file mapping from ch_mapping_file_rename
    path "prodigal_inital" from ch_inital_orfs
    file "rename_post_orf.py" from ch_rename_script
    
    output:
    path "prodigal_renamed" into ch_inital_orfs_renamed
    
    
    
    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd

    samples = pd.read_csv('${mapping}',sep='\t')
    subprocess.run(['mkdir prodigal_renamed'], shell=True)

    for index, row in samples.iterrows():
        stub = row['sequence-id']

        command = "python3 rename_post_orf.py -i prodigal_inital/"+stub+".fna -o prodigal_renamed/"+stub+"_renamed.fna"
        subprocess.run([command], shell=True)
    """

}

process ScreenORFover100{
    publishDir "${params.outdir}", mode: 'copy'

    container "docker://lorentzb/bioawk"

    input:
    file mapping from ch_mapping_file_screen_orf
    path "prodigal_renamed" from ch_inital_orfs_renamed
    file "screen_orf_over_100.sh" from ch_orf_screen_script
    
    output:
    path "orf_over_100" into ch_orfs_over_100
    
    
    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd

    samples = pd.read_csv('${mapping}',sep='\t')
    subprocess.run(['mkdir orf_over_100'], shell=True)

    for index, row in samples.iterrows():
        stub = row['sequence-id']

        screen_and_comb_command ="echo "+stub+" | bash screen_orf_over_100.sh"
        subprocess.run([screen_and_comb_command],shell=True)
    """


}

process RunCdHit{
    publishDir "${params.outdir}", mode: 'copy'

    container "docker://lorentzb/cdhit"

    input:
    file mapping from ch_mapping_file_cluster
    path "orf_over_100" from ch_orfs_over_100
    
    
    output:
    path "clustered_orf" into ch_clustered_orf
    
    
    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd

    samples = pd.read_csv('${mapping}',sep='\t')
    subprocess.run(['mkdir clustered_orf'], shell=True)

    for index, row in samples.iterrows():
        stub = row['sequence-id']

        cdhit_command = "cd-hit -i orf_over_100/"+stub+"_over_100.fasta -o clustered_orf/"+stub+" -c 0.95 -G 0 -aS 0.9 -g 1 -d 0 -T 8 -M 0"
        subprocess.run([cdhit_command],shell=True)
    """
}

process IndexAndAlignGeneLib{
    publishDir "${params.outdir}", mode: 'copy'

    container "docker://lorentzb/bowtie"

    input:
    
    file mapping from ch_mapping_gene_lib
    path "no_host" from ch_no_host_gene_lib
    path "clustered_orf" from ch_clustered_orf
    
    output:
    path "sams" into ch_sams_from_gene_lib
     

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd

    samples = pd.read_csv("${mapping}",sep='\t')
    subprocess.run(['mkdir sams'], shell=True)


    for index, row in samples.iterrows():
        stub = row['sequence-id']

        build_index = 'bowtie2-build clustered_orf/'+stub+' '+stub
        subprocess.run([build_index], shell=True)

        map_reads = "bowtie2 -x "+stub+" -1 no_host/"+stub+"_R1.fastq.gz -2 no_host/"+stub+"_R2.fastq.gz -S sams/"+stub+".sam" 
        subprocess.run([map_reads], shell=True)
    """
}

process ConvertSamsToBams{

    publishDir "${params.outdir}/gene_lib", mode: 'copy'

    container "docker://lorentzb/samtools"

    input:
    
    file mapping from ch_mapping_supported_genes
    path "sams" from ch_sams_from_gene_lib
    
    output:
    path "bams" into ch_bams_from_gene_lib

    script:
    """
    #!/usr/bin/env python3
    
    import subprocess
    import pandas as pd

    samples = pd.read_csv("${mapping}",sep='\t')
    
    subprocess.run(['mkdir bams'],shell=True)

    for index, row in samples.iterrows():
        stub = row['sequence-id']

        #convert sam to bam 
        sam_conv = 'samtools view -bS sams/'+stub+'.sam > bams/'+stub+'.bam'
        subprocess.run([sam_conv], shell=True)
    """
}

process FilterSupportedGenes{
    publishDir "${params.outdir}/gene_lib", mode: 'copy'

    container "docker://lorentzb/samtools"

    input:
    
    file mapping from ch_mapping_derep_gene_lib
    path "bams" from ch_bams_from_gene_lib
    
    output:
    path "derepped_bams" into ch_bams_dereplicated

    script:
    """
    #!/usr/bin/env python3
    
    import subprocess
    import pandas as pd

    samples = pd.read_csv("${mapping}",sep='\t')
    
    subprocess.run(['mkdir derepped_bams'],shell=True)

    for index, row in samples.iterrows():
        stub = row['sequence-id']

        #steps to sort and mark duplicate alignments

        collate_command = "samtools collate -o " +stub+"_collate.bam bams/"+stub+".bam"
        subprocess.run([collate_command], shell=True)

        fixmate_command = "samtools fixmate -m "+stub+"_collate.bam " +stub+"_fixmate.bam"
        subprocess.run([fixmate_command], shell=True)

        sort_command = "samtools sort -o " +stub+ "_positionsort.bam " +stub+"_fixmate.bam"
        subprocess.run([sort_command], shell=True)

        markdup_command = "samtools markdup -s "+stub+"_positionsort.bam " +stub+"_markdup.bam"
        subprocess.run([markdup_command], shell=True)

        #pull duplicate alignments (supported evidence) into separate bam file

        separate_alignments_command = "samtools view -F 0X400 -b "+stub+"_markdup.bam > derepped_bams/" +stub+"_only_dups.bam"
        subprocess.run([separate_alignments_command], shell=True)
               
    """

}

process ConvertBamsToFasta{ 
    publishDir "${params.outdir}/gene_lib", mode: 'copy'

    container "docker://lorentzb/samtools"

    input:
    
    file mapping from ch_mapping_gene_bam_to_fasta
    path "derepped_bams" from ch_bams_dereplicated
    file "de_duplicate.sh" from ch_de_dup_script
    
    output:
    path "geneLibrary" into ch_gene_library

    script:

    """
    #!/usr/bin/env python3
    
    import subprocess
    import pandas as pd

    samples = pd.read_csv("${mapping}",sep='\t')
    
    subprocess.run(['mkdir geneLibrary'],shell=True)

    for index, row in samples.iterrows():
        stub = row['sequence-id']

        #split reads into fasta files
        sort_command = 'samtools sort -n -m 5G -@ 2 derepped_bams/'+stub+'_only_dups.bam -o '+stub+'_sorted.bam'
        subprocess.run([sort_command], shell=True)

        split_command = 'samtools fasta -@ 8 '+stub+'_sorted.bam > '+stub+'_library.fasta'
        subprocess.run([split_command],shell=True)

        de_duplicate_fasta = 'echo '+stub+' | bash de_duplicate.sh'
        subprocess.run([de_duplicate_fasta],shell=True)
    """
}

process IntegrityandGeneLibStats{

    publishDir "${params.outdir}/gene_lib_stat", mode: 'copy'

    container "docker://lorentzb/integstat"

    input:
    
    file mapping from ch_mapping_lib_stat
    path "geneLibrary" from ch_gene_library
    
    
    output:
    path "stats" into ch_gene_library_stat

    script:

    """
    #!/usr/bin/env python3
    
    import subprocess
    import pandas as pd

    samples = pd.read_csv("${mapping}",sep='\t')
    subprocess.run(['mkdir stats'],shell=True)

    for index, row in samples.iterrows():
        stub = row['sequence-id']

        gtool_command = "python3 /opt/gtool-0.1.0/gtool.py -sg geneLibrary/"+stub+".fasta > "+stub+".log"
        subprocess.run([gtool_command],shell=True)

    """
}