# DEEP-CORE
Analysis of Metagenomic Data from four salt marsh sediment samples.
JGI PROJECT:  Combining high resolution organic matter characterization and microbial meta-omics to assess the effects of nutrient loading on salt marsh carbon sequestration
JGI LINK: https://genome.jgi.doe.gov/portal/Comhiguestration/Comhiguestration.info.html
JGI PROJECT ID: 503576

## Get the assemblies for each of the samples using the curl commands found in this file deep-core-get-assemblies.shx.. GLOBUS download would be a better alternative to all those curl scripts and you can access this using the JGI LINK above.
## Get the short reads for each of the samples from JGI. I would recommend usin GLOBUS
## Binning
### 1. fix the assemblies so that we are dealing with contig names that makes downstream programs happy and filter contigs that are below 2500 nt.
 
     #!/bin/bash
     #SBATCH --nodes=1
     #SBATCH --time=01:00:00
     #SBATCH --partition=express
     #SBATCH --array=1-55
     ASSEMBLY=$(sed -n "$SLURM_ARRAY_TASK_ID"p samples.txt)
     
     anvi-script-reformat-fasta ${ASSEMBLY}-final.contigs.fasta --simplify-names -o ${ASSEMBLY}_filter_contigs.fa -l 2500
     

### 2. Map the short reads using bowtie2 to reference assemblies ending in "filter_contigs.fa" created in step 1. This step was conducted so that only reads from nearby depths were mapped back to the assembly.  We organized the samples into eight different groups based on depth range. Below is a table of the group different groups

     110_120_130
     140_150_170
     190_200_210
     220_230_240
     270_290_310
     320_330_340
     370_390_410
     420_430_440
     
#### 2a. This is what the two SLURM mapping scripts look like. The first bash script (00_mapping_master.shx) is referenced in the second script.

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --time=6:00:00
    #SBATCH --mem=10Gb
    #SBATCH --partition=general

    bowtie2 -x ASSEMBLIES/${ASSEMBLY} -U QUALITY_READS/${READS} -q -S MAPPING/${ASSEMBLY}_${READS}.sam
    samtools view -bS -F 4 MAPPING/${ASSEMBLY}_${READS}.sam -o MAPPING/${ASSEMBLY}_${READS}.bam

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #!bin/bash
    #This is how I run the mapping for each of the asssemblies, but limited to reads that are derived from similar depths.

    #I made separate "assembly" and "reads" files that contain subsets of three adjacent depths.  This is not perfect, but should allow for increased stability during binning

    # here are the for loops i used to execute the mapping commands

    for ASSEMBLY in `cat assemblies_110_120_130.txt`; do for READS in `cat reads_110_120_130.txt`; do echo "${ASSEMBLY}, ${READS}"; export ASSEMBLY READS; sbatch 00_mapping_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_140_150_170.txt`; do for READS in `cat reads_140_150_170.txt`; do echo "${ASSEMBLY}, ${READS}"; export ASSEMBLY READS; sbatch 00_mapping_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_190_200_210.txt`; do for READS in `cat reads_190_200_210.txt`; do echo "${ASSEMBLY}, ${READS}"; export ASSEMBLY READS; sbatch 00_mapping_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_220_230_240.txt`; do for READS in `cat reads_220_230_240.txt`; do echo "${ASSEMBLY}, ${READS}"; export ASSEMBLY READS; sbatch 00_mapping_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_270_290_310.txt`; do for READS in `cat reads_270_290_310.txt`; do echo "${ASSEMBLY}, ${READS}"; export ASSEMBLY READS; sbatch 00_mapping_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_320_330_340.txt`; do for READS in `cat reads_320_330_340.txt`; do echo "${ASSEMBLY}, ${READS}"; export ASSEMBLY READS; sbatch 00_mapping_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_370_390_410.txt`; do for READS in `cat reads_370_390_410.txt`; do echo "${ASSEMBLY}, ${READS}"; export ASSEMBLY READS; sbatch 00_mapping_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_420_430_440.txt`; do for READS in `cat reads_420_430_440.txt`; do echo "${ASSEMBLY}, ${READS}"; export ASSEMBLY READS; sbatch 00_mapping_master.shx; sleep 1; done; done

### 3. The next step is to construct a contigs database for each of the assemblies.

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=06:00:00
    #SBATCH --mem=250Gb
    #SBATCH --partition=general
    ##SBATCH --array=1-55

    ASSEMBLY=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_assemblies.txt)
    anvi-gen-contigs-database -f ASSEMBLIES/${ASSEMBLY}_filter_contigs.fa -o ASSEMBLIES/${ASSEMBLY}_filter_contigs.db
    anvi-run-hmms -c ASSEMBLIES/${ASSEMBLY}_filter_contigs.db

### 4. Profile each of the mapping files using a similar approach to the mapping.  I created a list of the bam files for each of the depth categories to facilitate this process.
  
    #!/bin/bash
    #
    #SBATCH --nodes=10
    #SBATCH --tasks-per-node=20
    #SBATCH --time=24:00:00
    #SBATCH --mem=80Gb
    #SBATCH --partition=short

    anvi-profile -i MAPPING/${MAPPING_BAM} -c ASSEMBLIES/${ASSEMBLY}_filter_contigs.db -M 3000 -o MAPPING/${MAPPING_BAM}-PROFILE -T 40
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!bin/bash
     for ASSEMBLY in `cat assemblies_110_120_130.txt`; do for MAPPING_BAM in `cat bams_110_120_130.txt`; do echo "${ASSEMBLY}, ${MAPPING_BAM}"; export ASSEMBLY READS; sbatch 03_profiling_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_140_150_170.txt`; do for MAPPING_BAM in `cat bams_140_150_170.txt`; do echo "${ASSEMBLY}, ${MAPPING_BAM}"; export ASSEMBLY READS; sbatch 03_profiling_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_190_200_210.txt`; do for MAPPING_BAM in `cat bams_190_200_210.txt`; do echo "${ASSEMBLY}, ${MAPPING_BAM}"; export ASSEMBLY READS; sbatch 03_profiling_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_220_230_240.txt`; do for MAPPING_BAM in `cat bams_220_230_240.txt`; do echo "${ASSEMBLY}, ${MAPPING_BAM}"; export ASSEMBLY READS; sbatch 03_profiling_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_270_290_310.txt`; do for MAPPING_BAM in `cat bams_270_290_310.txt`; do echo "${ASSEMBLY}, ${MAPPING_BAM}"; export ASSEMBLY READS; sbatch 03_profiling_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_320_330_340.txt`; do for MAPPING_BAM in `cat bams_320_330_340.txt`; do echo "${ASSEMBLY}, ${MAPPING_BAM}"; export ASSEMBLY READS; sbatch 03_profiling_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_370_390_410.txt`; do for MAPPING_BAM in `cat bams_370_390_410.txt`; do echo "${ASSEMBLY}, ${MAPPING_BAM}"; export ASSEMBLY READS; sbatch 03_profiling_master.shx; sleep 1; done; done
    for ASSEMBLY in `cat assemblies_420_430_440.txt`; do for MAPPING_BAM in `cat bams_420_430_440.txt`; do echo "${ASSEMBLY}, ${MAPPING_BAM}"; export ASSEMBLY READS; sbatch 03_profiling_master.shx; sleep 1; done; done

### 5. Merge the profiles for each of the assemblies.

## Relative abundance for each of the MAGs.  

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
    
#### The above is a bit older... here is an update to the way I mapped each sample to the collection of contigs for all MAGs reconstructed in this study.  I'm running it from my scratch directory... /scratch/vineis.j/DEEP_CORE/ALL-MAGS

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --mem=100Gb
    #SBATCH --partition=short
    #SBATCH --array=1-55


    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p sample-names.txt)
    fastq_file=$( echo "/work/jennifer.bowen/JOE/DEEP-CORE-GLOBUS-        DOWNLOAD/${SAMPLE}/QC_Filtered_Raw_Data/"*"fastq")
    bbmap.sh nodisk=true interleaved=true ambiguous=random in=${fastq_file} ref=x_ALL_MAGS_CONCATENATED.fa out=MAPPING/${SAMPLE}.bam bamscript=to_bam.sh covstats=/MAPPING/covstats/${SAMPLE}.covstats.txt scafstats=/MAPPING/scafstats/${SAMPLE}.scafstats.txt

### 2. This command will result in a very useful file with an ending of *.scafstats.txt.* This file is key to calculating the relative abundance of each MAG across all samples. The first few lines looks like this. 

    #ID	Avg_fold	Length	Ref_GC	Covered_percent	Covered_bases	Plus_reads	Minus_reads	Read_GC	Median_fold	Std_Dev
    s_SalMarSW160110MG_3_8930491_Bin_10_1_contigs_c000000003170	0.0651	9210	0.4617	5.7220	527	2	2	0.4800	0	0.28
    s_SalMarSW160110MG_3_8930491_Bin_10_1_contigs_c000000015068	0.2270	3964	0.4435	22.0484	874	3	3	0.4322	0	0.43
    s_SalMarSW160110MG_3_8930491_Bin_10_1_contigs_c000000014165	0.3667	4091	0.4371	34.2459	1401	5	5	0.4340	0	0.53
    s_SalMarSW160110MG_3_8930491_Bin_10_1_contigs_c000000013048	0.2800	4285	0.4583	20.4667	877	4	4	0.4367	0	0.64
    s_SalMarSW160110MG_3_8930491_Bin_10_1_contigs_c000000018623	0.1703	3524	0.5477	14.7276	519	2	2	0.5400	0	0.43

### 3. We will use this file to collect the average fold coverage for each MAG using the script "estimate-MAG-coverage-from-bbmap-covstats-v2.py". Running this script requires a the scafstats.txt file specified above and a list of the scaffold ids in the fasta file (x_ALL_MAGS_CONCATENATED.fa). I run the script like this. What happens next: The script collects the values in the Avg_fold column and returns the average of this list for each MAG.   

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=100Gb
    #SBATCH --partition=express
    #SBATCH --array=1-55


    READS=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_file-names-of-reads-to-map.txt)
    python ~/scripts/estimate-MAG-coverage-from-bbmap-covstats-v2.py --mappingfile ${READS}.covstats.txt --out ${READS}-MAG-bbmap-Avg-fold.txt --ids x_all-MAG-ids.txt
    
##### This is what the first few lines of the resulting *MAG-bbmap-Avg-fold.txt* file looks like.

    a_s_SalMarSW160110MG_3_12718.2.279671.TCGGTTAC-GTAACCGA.covstats.txt	stat
    s_SalMarSW160110MG_3_8930491_Bin_10_1	17.9314
    s_SalMarSW160110MG_3_8930491_Bin_12_1	4.7454
    s_SalMarSW160110MG_3_8930491_Bin_12_2	5.8875
    s_SalMarSW160110MG_3_8930491_Bin_14_contigs	11.6489
    s_SalMarSW160110MG_3_8930491_Bin_16_contigs	14.0543
    
### 4. Therefore, we have an average of averages reported for each of our MAGs.  We need to correct this for the number of reads available in order to be able to compare the relative of a bin across samples. If we divide this number by the Mbp within each dataset, we have an estimate the *relative average fold coverage/Mbp*. The Mbp is equivalent to the "millions of base pairs" in a sample.  So if we have 2,000,000 nucleotides (I don't mean reads here, I mean nucleotides!) in our short read data, then the Mbp would be 2. Once I past together each of the "MAG-bbmap-Avg-fold.txt" file, I run this script to tabulate the relative fold coverage per sample. 

## Improved relative abundance for each of the MAGs. I created a script that will calculate the reads per Kbp per Mbp of sequenicng data (RPKM) for each MAG across all samples specified. Here is what the sbatch looks like. [reads mapped to MAG]/[total MAG length (/1000)] / [total nt in metagenomic sample/1000000]. This calculation is corrected for the side of the genome and the number of reads in the dataset allowing for comparison of MAGs within and across all samples.


#### Working from this directory /scratch/vineis.j/DEEP_CORE/ALL-MAGS 

    #!/bin/bash
    #
    #SBATCH --nodes=2
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=100Gb
    #SBATCH --partition=short
    #SBATCH --time=23:00:00
    python ~/scripts/estimate-MAG-coverage-from-bbmap-covstats-v3.py --mappingfiles /scratch/vineis.j/DEEP_CORE/ALL-MAGS/ --out x_ALL-MAG-RPKM-matrix-short-partition.txt --ids x_all-MAG-ids.txt --nts x_nt-per-sample.txt
    
## PHYLOGENOMICS: Here we outline the steps to generate the phylogenetic tree for the manuscript.  Step one is to create a list of the ribosomal proteins that we will use to create our phylogeny. Ribosomal proteins are routinely used to estimate phylogenetic relationships and here is a list of the genes that we use in a file called "x_gene-names.txt" and another file to get started in the external-genomes.txt file (both are contained in this repository . NOTE: I was working here /work/jennifer.bowen/DEEP-CORE/ALL-MAGS during the creation of this tutorial. This analysis (steps below) have been replaced to include the GEM genomes from intertidal and estuary sediments. The details of that analysis are found here. /work/jennifer.bowen/GEM-and-DEEP-CORE-for-estuary-paper

### 1. Start with running the command to pull out concatenated genes for each of the MAGs and then create the tree.  All this hapens using this amazing script.

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    SBATCH --time=10:00:00
    #SBATCH --mem=200Gb
    #SBATCH --partition=short

    #ANVIO 6.2

    anvi-get-sequences-for-hmm-hits --external-genomes -o x_concatenated-ribosomal-proteins-15.fa --hmm-source Campbell_et_al --gene-names x_gene-names.txt --return-best-hit --get-aa-sequences --concatenate --max-num-genes-missing-from-bin 15

    anvi-gen-phylogenetic-tree -f x_concatenated-ribosomal-proteins-15.fa -o x_phylogenomic-ribosomal-tree.tre
    

## ANI at the phylum level. Because the number of MAGs is so large, I partitioned the ANI analysis by phylum. 

### 1. I started by making text files that contain the names of the taxa for each of 15 major phyla identified in the DEEP-CORE based on the taxonomic assigment of the gtdbtk. I did not include the phyla that contained a single genome... Obviously!
    
    (base) [vineis.j@login-01 ALL-MAGS]$ cut -f 2 x_GTDBtk-OUTPUT-real-new-mags-included-multi-cpu/gtdbtk.ar122.summary.tsv | cut -f 2 -d ";" | sort | uniq -c
      1 classification
     26 p__Aenigmarchaeota
     51 p__Asgardarchaeota
     72 p__Crenarchaeota
     15 p__Euryarchaeota
      1 p__Halobacterota
     33 p__Nanoarchaeota
     12 p__Thermoplasmatota
    
    (base) [vineis.j@login-01 ALL-MAGS]$ cut -f 2 x_GTDBtk-OUTPUT-real-new-mags-included-multi-cpu/gtdbtk.bac120.summary.tsv | cut -f 2 -d ";" | sort | uniq -c
      1 classification
      4 p__AABM5-125-24
     65 p__Acidobacteriota
    127 p__Actinobacteriota
     34 p__Aerophobota
    178 p__Bacteroidota
     52 p__Bipolaricaulota
     17 p__Caldatribacteriota
     10 p__Calditrichota
      4 p__Campylobacterota
    524 p__Chloroflexota
     29 p__Chloroflexota_A
      5 p__Cloacimonadota
      1 p__Dadabacteria
      8 p__Deinococcota
      3 p__Dependentiae
    200 p__Desulfobacterota
      1 p__Desulfobacterota_A
      3 p__Eisenbacteria
      4 p__FCPU426
      1 p__Fermentibacterota
      1 p__Fibrobacterota
     15 p__Gemmatimonadota
      2 p__Hydrogenedentota
      3 p__Krumholzibacteriota
      4 p__KSB1
     18 p__Marinisomatota
      1 p__MBNT15
      1 p__Methylomirabilota
      6 p__Myxococcota
      3 p__Nitrospirota
      9 p__Omnitrophota
    147 p__Patescibacteria
    271 p__Planctomycetota
    116 p__Proteobacteria
      2 p__RBG-13-66-14
      1 p__SM23-31
     61 p__Spirochaetota
     24 p__TA06_A
     12 p__UBP3
      4 p__UBP7_A
      6 p__Verrucomicrobiota_A
     13 p__WOR-3
    112 p__WOR-3_B
     32 p__Zixibacteria
     
    cat x_BAC-PHYLA.txt x_ARC-PHYLA.txt > x_PHYLA-list.txt

#### I manually removed the phyla that contained a single MAG from the x_PHYLA-list.txt and also the lines that contained the word "classification"

### 2. Create a list for each of the phyla contained in the x_PHYLA-list.txt file.

    cat x_GTDBtk-OUTPUT-real-new-mags-included-multi-cpu/gtdbtk*12*.summary.tsv | sort > x_all_gtdbtk-output.txt
    
##### I manually removed the last two lines of the x_all_gtdbtk-output.txt file that contained header information.. then I created a separate external genomes file for each of the phyla.

    for i in `cat x_PHYLA-list.txt`; do grep $i x_all_gtdbtk-output.txt | cut -f 1 > x_$i-list.txt; done
    for i in `cat x_PHYLA-list.txt`; do python x_make-external-genomes-from-list.py x_$i-list.txt x_$i-external-genomes.txt; done
    
### 3. Run ANI for each of the external-genomes files.

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=10
    #SBATCH --time=24:00:00
    #SBATCH --mem=300Gb
    #SBATCH --partition=short
    #SBATCH --array=1-41

    MAGs=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_PHYLA-list.txt)

    anvi-compute-genome-similarity -e x_${MAGs}-external-genomes.txt -o x_ANI-${MAGs} -T 40
    anvi-dereplicate-genomes --ani-dir x_ANI-${MAGs}/ -o x_ANI_dereplication-${MAGs} --program fastANI --method ANIb --use-full-percent-identity --min-full-percent-identity 0.90 --similarity-threshold 0.95

