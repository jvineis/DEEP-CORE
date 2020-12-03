# DEEP-CORE
Analysis of Metagenomic Data from four salt marsh sediment samples.
## Lets start with relative abundance for each of the MAGs.  

### 1. Map the short reads from each metagenomic sample to the collection of scaffolds containing all scaffolds for all MAGs. This is the slurm script that I use, but you can see that bbmap was the mapper of choice here. These are reads derived from JGI, which are interleaved fastq files. The ref file is the fasta file containing all scaffolds for all MAGs.

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --mem=100Gb
    #SBATCH --partition=short
    #SBATCH --array=1-55


    READS=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_file-names-of-reads-to-map.txt)
    bbmap.sh nodisk=true interleaved=true ambiguous=random in=/scratch/vineis.j/DEEP_CORE/QUALITY_READS/${READS}.filter-METAGENOME.fastq ref=x_ALL_MAGS_CONCATENATED.fa    out=${READS}.bam bamscript=to_bam.sh covstats=${READS}.covstats.txt scafstats=${READS}.scafstats.txt

### 2. This command will result in a very useful file with an ending of *.scafstats.txt.* This file is key to calculating the relative abundance of each MAG across all samples. The first few lines looks like this. 

    #ID	Avg_fold	Length	Ref_GC	Covered_percent	Covered_bases	Plus_reads	Minus_reads	Read_GC	Median_fold	Std_Dev
    s_SalMarSW160110MG_3_8930491_Bin_10_1_contigs_c000000003170	0.0651	9210	0.4617	5.7220	527	2	2	0.4800	0	0.28
    s_SalMarSW160110MG_3_8930491_Bin_10_1_contigs_c000000015068	0.2270	3964	0.4435	22.0484	874	3	3	0.4322	0	0.43
    s_SalMarSW160110MG_3_8930491_Bin_10_1_contigs_c000000014165	0.3667	4091	0.4371	34.2459	1401	5	5	0.4340	0	0.53
    s_SalMarSW160110MG_3_8930491_Bin_10_1_contigs_c000000013048	0.2800	4285	0.4583	20.4667	877	4	4	0.4367	0	0.64
    s_SalMarSW160110MG_3_8930491_Bin_10_1_contigs_c000000018623	0.1703	3524	0.5477	14.7276	519	2	2	0.5400	0	0.43

### 3. We will use this file to collect the average fold coverage for each MAG using the script "estimate-MAG-coverage-from-bbmap-covstats-v2.py"
