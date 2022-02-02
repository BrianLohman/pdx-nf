# pdx-nf
`pdx-nf` is a A Nextflow workflow to generate gene counts from PDX bulk RNAseq data  

# Usage
`nextflow run pdx-nf --fastq [path to fastqs]`

# Arguments:
+ `--fastq` 
    + Path to fastqs, default is `./fastq`

# Outline
This workflow will:
+ Make channel from paired end reads in `--fastq`
+ Remove optical duplicates with Clumpify
+ Trim adapters with Cutadapt (hard coded adapter sequence)
+ Map reads with STAR
+ Filter BAMs for alignments that map to the human reference. Index BAMS and calculate idxstats with Samtools. Also published to oudir/bams
+ Count genes with featureCounts. Published to `outdir/feature_counts`
+ RNAseq metrics with Picard. Published to `outdir/qc_metrics`
+ FastqScreen to look for contamination. Published to `outdir/qc_metrics`
+ FastQC. Published to `outdir/qc_metrics`

# Results
Workflow results are published to `--outdir`, whose default is `./results`  
These include:
+ BAMs from STAR filtered for human reads to `./results/bams`
+ Samtools idxstats to `./results/bams`
+ featureCounts to `./results/feature_counts`
+ RNAseq quality metrics (FastqScreen, FastQC, Picard RNAseqMetrics) to `./results/qc_metrics`

# Defaults
Default values for paths to reference files are in `nextflow.config`.  
Default STAR Index is a combination of Human and Mouse for 150 base pair reads  
Human genome is GRCh38 and Ensembl v96  
Mouse genome is GRCm38 and Ensembl v96  