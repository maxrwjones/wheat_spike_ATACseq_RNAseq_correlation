# Correlating developmental series of ATAC-seq and RNA-seq data

BASH and R scripts for Chapter 4 of my PhD thesis "Exploring genetic loci for the improvement of inflorescence traits in the polyploid cereals wheat and tef".

This code was used to explore whether, across wheat spike development, correlations could be identified between the accessibility of accessible chromatin regions (ACRs) and the expression of neighbouring genes.

_Note:
- I use the term CACRs (consensus ACRs) here because I am interested in ACRs which were detected across all developmental stages analysed
- I also call CACRs 'peaks' as a shorthand in some scripts
- BASH scripts are set up to run using the SLURM job handler. SLURM lines can be deleted if running locally or replaced as needed. _

### Initial set-up
The code can be executed linearly from script 1 to script 20 and requires the user to set up a UNIX directory with the following sub-directories and files:
- scripts (containing scripts 1-20 from this pipeline)
- logs
- outputs
- inputs

The 'inputs' sub-directory should contain:
- A genome annotation file. Here I have used a file called "IWGSC_v1.1_HC_20170706.sorted.gff3".
- A set of .narrowPeak files containing ACRs called from ATAC-seq data. One file per developmental stage. These should be placed inside a subfolder called 'MACS2_peaks'.
- A set of .bedGraph files of the trimmed, mapped, and cleaned ATAC-seq data which have been normalised for final library depth. These should be placed inside a subfolder called 'scaled_read_coverage'.
- A set of .tsv files from the pseudomapper Kallisto. Each file represents one rep from one developmental stage and should give count and tpm data for each gene in the annotation. These should be placed inside a subfolder called 'RNAseq_abundance_files'.
- A file containing unmethylated region (UMR) calls. Here I have used a file called "CS_UMRs_cov5_unsplit_coords.bed"


### Description of operations
1) Streamline original peak files
2) Identify consensus ACRs (CACRs) and retrieve fold-change and q-value stats for CACRs
3) Make a version of the consensus peak file which can be viewed using the Integrative Genomics Viewer
4) Filter GFF3 to keep only genes then generate extended coordinates for each gene (except in certain edge cases as explained in script)
5) Make single gene BEDs
6) Extract maximum read coverage values under each CACR for each replicate of each stage
7) Make single "consensus_peaks_with_max_cov.bed" from stagewise max cov files
8) Intersect CACRs with single gene BEDs to assign them to genes
9) Test genes for differential expression across developmental stages (ImpulseDE2)
10) Test CACRs for differential chromatin accessibility across developmental stages (ImpulseDE2)
11) Identify which peaks intersect UMRs and then test these for differential chromatin accessibility (ImpulseDE2)
12) Enrichment analysis of DEGs for CACRs and dCACRs, then test simulation of random CACR-gene pairs
13) Enrichment analysis of DEGs for UMR-supported CACRs and dCACRs, then test simulation of random CACR-gene pairs
14) Normalise CACRs and gene trajectories to the same scale to enable correlation
15) Calculate sums of squares difference between CACRs and genes to find **enhancer** candidates
16) Calculate sums of squares difference between CACRs and genes to find **silencer** candidates
17) Simulate **enhancer** candidates
18) Simulate **silencer** candidates
19) Compare simulated **enhancer** candidates with 'real' candidates. Is there any capacity to detect real signal?
20) Compare simulated **silencer** candidates with 'real' candidates. Is there any capacity to detect real signal?
