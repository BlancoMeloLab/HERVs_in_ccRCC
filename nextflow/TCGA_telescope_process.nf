#!/usr/bin/env nextflow
/*
Pipeline for downloading and processing TCGA gene expression files with telescope
SK 2025-08-11
Modified: 2026-02-03 - Added concatenation of STAR gene counts and Telescope TE counts

check that manifest and token are good; should put in a work dir on fast

run nextflow like:
-- from dir -- need relative path to nf
-- cd into work dir
-- may need to create tempdir depending on cluster set up

module load STAR
module load Telescope
module load gdc-client/1.6.1-GCCcore-10.2.0 
module load Python/3.8.2-GCCcore-9.3.0 nextflow
module load SAMtools/1.17-GCC-12.2.0

To run nextflow:
nextflow run -resume TCGA_telescope_process.nf --manifest gdc_manifest_test.txt --token <token>

Pipeline steps:

1. Download BAM files from TCGA
    1a. GDC Data Transfer Tool Client -- command line tool to access genomic data commons and download BAM files; will need authentication token
2. Convert to fastq using samtools 
3. Map with STAR
    3a. high multimapping tolerance
4. Mapping w/ Telescope
    4a. hg19 ref gtf + 8k HERVs + LINEs 
5. Concatenate STAR gene counts with Telescope TE counts

*/

// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
// Parameter values can be provided by the user in a params.json or using flags e.g. --manifest /path/to/manifest.txt
params.manifest = false
params.token = false
params.gtf = "" // gtf
params.starin = "" // star index
params.outdir = "" // outdir


// Processes
// 1. Download BAM
process gdc_download {
    input: 
        path manifest
        path token
    output:
        path "bams/*.bam"

    """
    gdc-client download -m "${manifest}" -t "${token}"
    mkdir bams
    find -type f -iname "*.bam" -exec mv --backup=numbered -t bams/ {} +
    """
}

// 2. Revert to fastq, updated to hpc temp
process bam_2_fastq {
    input:
        path X
    output:
        tuple env(OUTPUT_NAME), path("*.fq.gz")
    script:
        """#!/bin/bash
        set -euo pipefail

        OUTPUT_NAME="${X.name.replaceAll(/.bam$/, '')}"
        OUTPUT_NAME_R1="\${OUTPUT_NAME}_R1.fq.gz"
        OUTPUT_NAME_R2="\${OUTPUT_NAME}_R2.fq.gz"

        echo "OUTPUT_NAME: \$OUTPUT_NAME"
        echo "OUTPUT_NAME_R1: \$OUTPUT_NAME_R1"
        echo "OUTPUT_NAME_R2: \$OUTPUT_NAME_R2"

        samtools collate -u -T /tempdir/temp -O "$X" | samtools fastq -1 "\$OUTPUT_NAME_R1" -2 "\$OUTPUT_NAME_R2" -T /tempdir/temp -n -
        """
}

// 3. Map with STAR - MODIFIED to output ReadsPerGene file
process star_map {
    publishDir "${params.outdir}/star_aln", mode: 'copy', pattern: "star_aln/${sample_id}.sorted.bam*"
    publishDir "${params.outdir}/star_counts", mode: 'copy', pattern: "star_aln/${sample_id}.ReadsPerGene.out.tab"

    input:
        path starin
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("star_aln/${sample_id}.sorted.bam"), path("star_aln/${sample_id}.sorted.bam.bai"), emit: bam
        tuple val(sample_id), path("star_aln/${sample_id}.ReadsPerGene.out.tab"), emit: gene_counts

    script:
    """
    mkdir -p star_aln

    STAR --genomeDir ${starin} \
         --runThreadN 24 \
         --readFilesCommand gunzip -c \
         --readFilesIn *_R1.fq.gz *_R2.fq.gz \
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

// 4. Map with Telescope - MODIFIED to output te_counts.tsv specifically
process telescope {
    publishDir "${params.outdir}/telescope", mode: 'copy'

    input:
        tuple val(sample_id), path(bamfile), path(bai), path(gtf)
    output:
        tuple val(sample_id), path("tel_out/${sample_id}/${sample_id}-TE_counts.tsv"), emit: te_counts
        path "tel_out/${sample_id}/*", emit: all_outputs

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

// 5. NEW PROCESS: Concatenate STAR gene counts with Telescope TE counts
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
    
    # Add STAR gene counts (column 2 is unstranded count, column 1 is gene ID)
    # Skip first 4 header lines from STAR output
    tail -n +5 ${gene_counts} | awk 'BEGIN{OFS="\t"}{print \$1, \$2}' >> ${sample_id}.combined_counts.tsv
    
    # Add Telescope TE counts
    # Assuming telescope report has header, skip it and add counts
    # Adjust column numbers based on actual telescope output format
    tail -n +2 ${te_counts} | awk 'BEGIN{OFS="\t"}{print \$1, \$2}' >> ${sample_id}.combined_counts.tsv
    
    echo "Combined counts created for ${sample_id}"
    echo "Total features: \$(wc -l < ${sample_id}.combined_counts.tsv)"
    """
}

// Workflow
workflow {
    manifest = file(params.manifest, checkIfExists: true)
    token = file(params.token, checkIfExists: true)
    gtf = file(params.gtf, checkIfExists: true)
    starin = file(params.starin, checkIfExists: true)

    // Download BAMs
    gdc_download(manifest, token)
    
    // Convert to FASTQ
    bam_2_fastq(gdc_download.out.flatten())
    
    // Map with STAR (now outputs both BAM and gene counts)
    star_map(starin, bam_2_fastq.out)
    
    // Map with Telescope (needs BAM, BAI, and GTF)
    telescope(star_map.out.bam.map { sample_id, bam, bai -> tuple(sample_id, bam, bai, gtf) })
    
    // Combine the gene counts and TE counts for each sample
    // Join the two channels on sample_id
    combined_input = star_map.out.gene_counts.join(telescope.out.te_counts)
    
    // Concatenate counts
    concatenate_counts(combined_input)
}

