#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false
if (params.help) {
    log.info"""
    -----------------------------------------------------------------------------
    pdx-nf: a bulk RNAseq workflow for PDX samples
    =============================================================================
    Required arguments:
    -------------------
    --fastq         Full path to directory with reads. Default: ./fastq
   
    Description:
    ------------
    Aligns bulk RNAseq samples to a combined human and mouse reference.
    Assumes only one fastq pair per sample.
    Keep only those reads that align to the human genome and count genes.
    All reference files can be set via params in nextflow.config.

    -----------------------------------------------------------------------------
    """.stripIndent()
    exit 0
}

// Logging
log.info("\n")
log.info("Fastq directory       (--fastq)                 :${params.fastq}")
log.info("STAR Index            (--star_index)            :${params.star_index}")
log.info("Results               (--outdir)                :${params.outdir}")

// Read pair channel: tuple of pair_id and paired fastqs
Channel
    .fromFilePairs( params.fastq )
    .ifEmpty { error "Cannot find any reads matching: ${params.fastq}" }
    .set { read_pairs_ch }

// Remove optical duplicates from fastq
process dedup {
    tag "$pair_id"

    input:
      tuple val(pair_id), path(reads)

    output:
      tuple val(pair_id), path("${pair_id}.read*.fastq.gz"), emit: deduped_reads
      path("${pair_id}.clumpify.out.txt")

    script:
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/BBmap/v38.34/clumpify.sh \
      in1=${pair_id}_R1_001.fastq.gz \
      in2=${pair_id}_R2_001.fastq.gz \
      out1=${pair_id}.read1.fastq.gz \
      out2=${pair_id}.read2.fastq.gz \
      dupedist=12000 \
      dedupe=t \
      optical=t 2> ${pair_id}.clumpify.out.txt
      """
}

// Trim adapters
process trim {
    module 'cutadapt/2.10'
    tag "$pair_id"

    publishDir "${params.outdir}/trimmed_reads", mode:"copy"
    
    input:
      tuple val(pair_id), path(deduped_reads)

    output:
      tuple val(pair_id), path("${pair_id}.trim*.fastq.gz"), emit: trimmed_reads
      path("${pair_id}.cutadapt.out.txt")

    script:
      """
      cutadapt \
      -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      -o ${pair_id}.trim1.fastq.gz \
      -p ${pair_id}.trim2.fastq.gz \
      -j 8 \
      -m 20 \
      ${pair_id}.read1.fastq.gz ${pair_id}.read2.fastq.gz > ${pair_id}.cutadapt.out.txt
      """
}

// Maps each read-pair with STAR
process star {
    tag "$pair_id"
      
    input:
      path(star_index)
      tuple val(pair_id), path(reads)
  
    output:
      tuple val(pair_id), path("${pair_id}_unsorted.bam"), emit: bam

    script:
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/STAR/2.7.9a/STAR --runMode alignReads \
      --genomeDir $star_index \
      --twopassMode Basic \
      --readFilesIn $reads \
      --readFilesCommand zcat \
      --runThreadN ${task.cpus} \
      --outSAMtype BAM Unsorted \
      --quantMode TranscriptomeSAM 
      mv Aligned.out.bam ${pair_id}_unsorted.bam
      """
}

// Filter alignments for human reads
process filter_alignment {
    tag "${pair_id}"
 
    publishDir "${params.outdir}/bams", mode:"copy"

    input:
      tuple val(pair_id), path(bam)

    output:
      tuple val(pair_id), path("${pair_id}_human_sorted.bam"), path("${pair_id}_human_sorted.bam.bai"), emit: indexed_bam
      path("${pair_id}.idxstats")

    script:
      """
        /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.10/samtools view -q 255 -h ${pair_id}_unsorted.bam | grep human | sed -e 's/human_//' > ${pair_id}_human.sam
        /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.10/samtools view -S -b ${pair_id}_human.sam > ${pair_id}_human.bam
        /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.10/samtools sort ${pair_id}_human.bam -o ${pair_id}_human_sorted.bam
        /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.10/samtools index ${pair_id}_human_sorted.bam
        /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.10/samtools idxstats ${pair_id}_human_sorted.bam | sort -V > ${pair_id}.idxstats
      """
}

// FeatureCounts
process feature_counts {
    tag "${pair_id}"
 
    publishDir "${params.outdir}/feature_counts", mode:"copy"

    input:
      tuple val(pair_id), path(bam), path(bai)
      path(gtf)

    output:
      path("${pair_id}.counts")

    script:
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/Subread/1.6.3/bin/featureCounts \
      -T 8 \
      -p \
      -C \
      -s 2 \
      --largestOverlap \
      -a $gtf \
      -o ${pair_id}.counts \
      $bam
      """
}

// FastqScreen
process fastqscreen {
    tag "${pair_id}"
 
    publishDir "${params.outdir}/qc_metrics", mode:"copy"

    input:
      tuple val(pair_id), path(reads)

    output:
    path("${pair_id}.trim1_screen.html")
    path("${pair_id}.trim1_screen.txt")

    script:
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/fastq_screen/v0.14.0/fastq_screen \
      --conf /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/FastQ_Screen_Genomes/redwood_fastq_screen.conf \
      --subset 1000000 ${pair_id}.trim1.fastq.gz -outdir ./
      """
}

// FastQC
process fastqc {
    tag "${pair_id}"
 
    publishDir "${params.outdir}/qc_metrics", mode:"copy"

    input:
      tuple val(pair_id), path(reads)
      
    output:
      path("${pair_id}.trim1_fastqc.html")
      path("${pair_id}.trim2_fastqc.html")
      path("${pair_id}.trim1_fastqc.zip")
      path("${pair_id}.trim2_fastqc.zip")

    script:
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/FastQC/0.11.5/fastqc -T ${task.cpus} -f fastq ${pair_id}.trim1.fastq.gz
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/FastQC/0.11.5/fastqc -T ${task.cpus} -f fastq ${pair_id}.trim2.fastq.gz
      """
}

// Picard RNAseq Metrics
process rnaseq_metrics {
    tag "${pair_id}"
 
    publishDir "${params.outdir}/qc_metrics", mode:"copy"

    input:
      tuple val(pair_id), path(bam), path(bai)

    output:
      path("${pair_id}.rna_metrics")

    script:
      """
      /usr/bin/java -Xmx20G -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/picard/2.9.0/picard.jar CollectRnaSeqMetrics  REF_FLAT=/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/release104/Homo_sapiens.GRCh38.104.refflat \
      STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/release104/Homo_sapiens.GRCh38.104.rRNA.interval \
      I=${pair_id}_human_sorted.bam  O=${pair_id}.rna_metrics
      """
}

workflow {
    dedup(read_pairs_ch)
    trim(dedup.out.deduped_reads)
    star(params.star_index, trim.out.trimmed_reads)
    filter_alignment(star.out.bam)
    feature_counts(filter_alignment.out.indexed_bam, params.gtf)
    fastqscreen(trim.out.trimmed_reads)
    fastqc(trim.out.trimmed_reads)
    rnaseq_metrics(filter_alignment.out.indexed_bam)
}
