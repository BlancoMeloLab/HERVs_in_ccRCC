#!/usr/bin/env nextflow
// for paired-end samples

nextflow.enable.dsl = 2

// Parameters
params.reads   = "./raw/*_{R1,R2}_001.fastq.gz"   // paired-end reads: two files per sample
params.gtf     = "/fh/fast/blancomelo_d/grp/telescope_pipeline/ref/GRCh38_GENCODE_rmsk_TE.gtf"
params.starin  = "/fh/fast/blancomelo_d/grp/telescope_pipeline/ref/hg38_gc_star_index_020926"
params.outdir  = "./results"
params.trim_len = 101   // shortest library in the cohort; longer reads are cut to this length

//
// Processes
//
process fastp_trim {
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.{json,html}"

    input:
        tuple val(sample_id), path(read_files)

    output:
        tuple val(sample_id), path("${sample_id}_trim_R?.fastq.gz"), emit: reads
        path "${sample_id}.fastp.json", emit: json
        path "${sample_id}.fastp.html", emit: html

    script:
    """
    fastp \
        --in1 ${read_files[0]} --in2 ${read_files[1]} \
        --out1 ${sample_id}_trim_R1.fastq.gz \
        --out2 ${sample_id}_trim_R2.fastq.gz \
        --detect_adapter_for_pe \
        --max_len1 ${params.trim_len} --max_len2 ${params.trim_len} \
        --length_required 50 \
        --thread 8 \
        --json ${sample_id}.fastp.json \
        --html ${sample_id}.fastp.html
    """
}

process star_map {
    publishDir "${params.outdir}/star_aln", mode: 'copy', pattern: "star_aln/${sample_id}.sorted.bam*"
    publishDir "${params.outdir}/star_counts", mode: 'copy', pattern: "star_aln/${sample_id}.ReadsPerGene.out.tab"

    input:
        path starin
        tuple val(sample_id), path(read_files)

    output:
        tuple val(sample_id), path("star_aln/${sample_id}.sorted.bam"), path("star_aln/${sample_id}.sorted.bam.bai"), emit: bam
        tuple val(sample_id), path("star_aln/${sample_id}.ReadsPerGene.out.tab"), emit: gene_counts

    script:
    """
    mkdir -p star_aln

    STAR --genomeDir ${starin} \
         --runThreadN 24 \
         --readFilesCommand gunzip -c \
         --readFilesIn ${read_files[0]} ${read_files[1]} \
         --outFileNamePrefix star_aln/tmp_ \
         --quantMode GeneCounts \
         --outSAMattributes NH HI NM MD AS XS \
         --outFilterMultimapNmax 500 \
         --winAnchorMultimapNmax 500 \
         --limitOutSAMoneReadBytes 1000000 \
         --outSAMtype BAM SortedByCoordinate

    # Rename to include sample ID
    mv star_aln/tmp_Aligned.sortedByCoord.out.bam star_aln/${sample_id}.sorted.bam
    mv star_aln/tmp_ReadsPerGene.out.tab star_aln/${sample_id}.ReadsPerGene.out.tab
    
    # Index for Telescope
    samtools index star_aln/${sample_id}.sorted.bam
    """
}

process telescope {
    publishDir "${params.outdir}/telescope", mode: 'copy'

    input:
        tuple val(sample_id), path(bamfile), path(bai)
        path gtf

    output:
        tuple val(sample_id), path("tel_out/${sample_id}/${sample_id}-TE_counts.tsv"), emit: te_counts
        path "tel_out/${sample_id}/*", emit: all_outputs

    script:
    """
    mkdir -p tel_out/${sample_id}
    telescope assign \
        --exp_tag ${sample_id} \
        --attribute gene_id \
        --reassign_mode choose \
        --theta_prior 200000 \
        --max_iter 200 \
        --updated_sam \
        --outdir tel_out/${sample_id} \
        ${bamfile} \
        ${gtf}
    """
}

process concatenate_counts {
    publishDir "${params.outdir}/combined_counts", mode: 'copy'

    input:
        tuple val(sample_id), path(gene_counts), path(te_counts)
    
    output:
        path "${sample_id}.combined_counts.tsv"

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Create header
    echo -e "feature_id\t${sample_id}" > ${sample_id}.combined_counts.tsv
    
    # Add STAR gene counts (column 2 is unstranded count, column 1 is gene ID, col 4 is reverse strand)
    # Skip first 4 header lines from STAR output
    tail -n +5 ${gene_counts} | awk 'BEGIN{OFS="\t"}{print \$1, \$4}' >> ${sample_id}.combined_counts.tsv
    
    # Add Telescope TE counts
    # Assuming telescope report has header, skip it and add counts
    # Adjust column numbers based on actual telescope output format
    tail -n +2 ${te_counts} | awk 'BEGIN{OFS="\t"}{print \$1, \$2}' >> ${sample_id}.combined_counts.tsv
    
    echo "Combined counts created for ${sample_id}"
    echo "Total features: \$(wc -l < ${sample_id}.combined_counts.tsv)"
    """
}

//
// Workflow definition
//
workflow {
    main:
        // Input channels for paired-end reads
        def read_files = channel
            .fromFilePairs(params.reads, flat: true)
            .map { sample_id, read1, read2 ->
                // Remove common suffixes from sample ID if needed
                def clean_id = sample_id.replaceAll(/_001$/, '')
                tuple(clean_id, [read1, read2])
            }

        def gtf = file(params.gtf)
        def starin = file(params.starin)

        // Trim adapters and cap read length
        fastp_trim(read_files)

        // Map with STAR using trimmed reads
        star_map(starin, fastp_trim.out.reads)
        
        // Map with Telescope (needs BAM, BAI, and GTF)
        telescope(star_map.out.bam, gtf)
        
        // Combine the gene counts and TE counts for each sample
        // Join the two channels on sample_id
        def combined_input = star_map.out.gene_counts.join(telescope.out.te_counts)
        
        // Concatenate counts
        concatenate_counts(combined_input)
}
