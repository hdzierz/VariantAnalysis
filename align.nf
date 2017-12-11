#!/usr/bin/env nextflow


/*
 * The reference genome file
  */
genome_file = file(params.genome)

/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs }
  
/*
 * Indexing target genome. 
 * Config in nextflow.config:
 *  tool: bwa
 *  input: File genome_file
 *  output: Channel genome_index
 *  params.build_index: true
 */

process build_index {
    module = params.build_index_module

    input:
    file genome_file from genome_file
      
    output:
    file 'genome.index*' into genome_index

    script:
    """
        bwa index -a bwtsw ${genome_file} -p genome.index
    """
}


label_align = 'align'

/*
 * Aligning raw reads against reference genome
 * Config in nextflow.config:
 *  tool: bwa
 *  input: File genoem_file, Channel read_pairs
 *  output: file bam
*/ 

process align{
    module = params.align_module  
    publishDir params.pub_align, mode: "move", overwrite: true, pattern: "*.bam"

    input:
    file genome_file from genome_file
    file genome_index from genome_index
    set pair_id, file(reads) from read_pairs

    output:
    set pair_id, file("${label_align}_${pair_id}.bam") into bam_align

    script:
    """
        bwa mem -t 8 -M genome.index ${reads} | samtools view -u - | samtools sort -O bam -T genome_file -o ${label_align}_${pair_id}.bam
    """
}

