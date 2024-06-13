# tcircumcincta_genome - genome-wide analyses

### Stephen Doyle

## Overall aims
- to identifiy regions of the genome containing genes associated with drug treatment response
- analysing several different datasets, each with different sets fo samples that response differently to drug treatment
- analysing all in a similar way will hopefully reveal drug-treatment specific QTLs


## Sample sets
- "Choi"
    - reanalysis of published data from a backcross between suceptible and multidrug resistant stains
- "Farms"
    - from Jenni's PhD, focused on sampling pre and post treatment of two UK farms
- "Strains"
    - archival strains Tci1, Tci5 and Tci7



```bash
# strain data sequenced at Sanger, so need to retrieve it from iRODs
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA

module load irods_extractor/v3.0.10

irods_extractor --studyid 7548 --runid 48426
```




## mapping
```bash
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ

module load mapping-helminth/v1.0.8

bsub.py --queue yesterday 1 mapping_poolseq "mapping-helminth --input sample_manifest.txt --reference teladorsagia_circumcincta_tci2_wsi3.0.genome.fa"

```

where "sample_manifest.txt" is: 
```bash
ID,R1,R2
F2_POST,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Farm2_Post_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Farm2_Post_2.fq.gz
F2_PRE,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Farm2_Pre_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Farm2_Pre_2.fq.gz
F3_POST_A,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Farm3_PostA_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Farm3_PostA_2.fq.gz
F3_POST_B,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Farm3_PostB_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Farm3_PostB_2.fq.gz
F3_PRE_A,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Farm3_PreA_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Farm3_PreA_2.fq.gz
F3_PRE_B,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Farm3_PreB_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Farm3_PreB_2.fq.gz
RS3,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/RS3_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/RS3_2.fq.gz
SINBRED,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Sinbred_1_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/Sinbred_1_2.fq.gz
MTci1_pool_adultMF_post-BZ,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#1_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#1_2.fastq.gz
MTci2_pool_L4,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#2_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#2_2.fastq.gz
MTci2_pool_adultMF,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#3_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#3_2.fastq.gz
MTci5_pool_L4_MOTRI,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#4_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#4_2.fastq.gz
MTci5_pool_adultMF_post-BZ,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#5_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#5_2.fastq.gz
MTci5_pool_adultM_post-IVM,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#6_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#6_2.fastq.gz
MTci7_pool_adultMF_MOX-R,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#7_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA/results/fastqs/48426_1#7_2.fastq.gz
```


```bash
# when mapping is completed, run multiqc

multiqc .

```



## Setup for running popgen analyses

```bash
mkdir /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/ANALYSIS

# go to mapping directory


# nice way to find all the bams, and them make symbolic links to a target directory
OUTPUT_DIR=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/ANALYSIS
find ~+ -type f -name '*.bam*' -exec ln -vs "{}" $OUTPUT_DIR/ ';'



cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/ANALYSIS

ln -s ../teladorsagia_circumcincta_tci2_wsi3.0.genome.fa


```




## Run between sample Fst comparisons using grenedalf
- decided to run each group of samples seperately, as they are quite distinct datasets both biologically and technically (sequencing platform, read lengths, amount of data, pool sizes etc), and so there are all sorts of biases that would not be controlled for by comparing them

```bash
module load grenedalf/0.3.0

# choi data
bsub.py 1 grendalf_choi_fst \
"grenedalf fst \
--reference-genome-fasta-file teladorsagia_circumcincta_tci2_wsi3.0.genome.fa \
--allow-file-overwriting \
--method unbiased-nei \
--file-prefix tc_choi_poolseq \
--sam-min-map-qual 30 \
--sam-min-base-qual 30 \
--filter-sample-min-count 2 \
--filter-sample-min-coverage 10 \
--filter-sample-max-coverage 100 \
--pool-sizes 800 \
--window-type sliding \
--window-sliding-width 100000 \
--write-pi-tables \
--separator-char tab \
--sam-path RS3.bam \
--sam-path SINBRED.bam"


# farm data
bsub.py 1 grendalf_farm_fst \
"grenedalf fst \
--reference-genome-fasta-file teladorsagia_circumcincta_tci2_wsi3.0.genome.fa \
--allow-file-overwriting \
--method unbiased-nei \
--file-prefix tc_jm_farm_poolseq \
--sam-min-map-qual 30 \
--sam-min-base-qual 30 \
--filter-sample-min-count 2 \
--filter-sample-min-coverage 10 \
--filter-sample-max-coverage 100 \
--pool-sizes 182 \
--window-type sliding \
--window-sliding-width 100000 \
--write-pi-tables \
--separator-char tab \
--sam-path F2_POST.bam \
--sam-path F2_PRE.bam \
--sam-path F3_POST_A.bam \
--sam-path F3_POST_B.bam \
--sam-path F3_PRE_A.bam \
--sam-path F3_PRE_B.bam"


# strain data
bsub.py 1 grendalf_strains_fst \
"grenedalf fst \
--reference-genome-fasta-file teladorsagia_circumcincta_tci2_wsi3.0.genome.fa \
--allow-file-overwriting \
--method unbiased-nei \
--file-prefix tc_strains_poolseq \
--sam-min-map-qual 30 \
--sam-min-base-qual 30 \
--filter-sample-min-count 2 \
--filter-sample-min-coverage 10 \
--filter-sample-max-coverage 100 \
--pool-sizes 200 \
--window-type sliding \
--window-sliding-width 100000 \
--write-pi-tables \
--separator-char tab \
--sam-path MTci1_pool_adultMF_post-BZ.bam \
--sam-path MTci2_pool_adultMF.bam \
--sam-path MTci2_pool_L4.bam \
--sam-path MTci5_pool_adultMF_post-BZ.bam \
--sam-path MTci5_pool_adultM_post-IVM.bam \
--sam-path MTci5_pool_L4_MOTRI.bam \
--sam-path MTci7_pool_adultMF_MOX-R.bam"
```


## 

- farm2 - pre v post
- farm3 - pre v post
- farm2 post v farm3 post
- choi RS3 v Sinbred

```R

library(tidyverse)
library(patchwork)

colours <- c("#d0e11c", "#73d056", "#2db27d", "#21918c", "#2e6e8e", "#46327e")

jm100k <- read.table("tc_jm_farm_poolseqfst.csv", header=T)
jm100k <- jm100k %>% filter(grepl("chr_[12345X]", chrom)) %>% arrange(chrom, start)


jm100k_F2 <- jm100k %>% select(chrom, start, end, snps, F2_POST.1.F2_PRE.1)
jm100k_F2$name <- "1.F2_POST_F2_PRE"
colnames(jm100k_F2) <- c("chrom", "start", "end", "snps", "fst", "name")

jm100k_F3 <- jm100k %>% select(chrom, start, end, snps, F3_POST_A.1.F3_PRE_A.1)
jm100k_F3$name <- "2.F3_POST_F3_PRE"
colnames(jm100k_F3) <- c("chrom", "start", "end", "snps", "fst", "name")

jm100k_postpost <- jm100k %>% select(chrom, start, end, snps, F2_POST.1.F3_POST_A.1)
jm100k_postpost$name <- "3.F2_POST_F3_POST_A"
colnames(jm100k_postpost) <- c("chrom", "start", "end", "snps", "fst", "name")

gw_sig <- jm100k %>% summarise(gw=mean(F3_PRE_A.1.F3_PRE_B.1) + 5*sd(F3_PRE_A.1.F3_PRE_B.1))

data <- bind_rows(jm100k_F2, jm100k_F3, jm100k_postpost)

plot_jm100k <- ggplot(data, aes((start+50000)/1e6, fst, col=chrom)) + 
    geom_point(size=0.1) + 
    labs(title="") +
    facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
    theme_bw() + 
    theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
    labs(title="A", x="Genomic position (Mb)", y="Fst") +
    scale_colour_manual(values=colours) +
    geom_hline(yintercept=gw_sig$gw, linetype='dashed')


plot_jm100k

choi100k <- read.table("tc_choi_poolseqfst.csv", header=T)
choi100k <- choi100k %>% filter(grepl("chr_[12345X]", chrom)) %>% arrange(chrom, start)

choi100k$name <- "RS3.SINBRED"
colnames(choi100k) <- c("chrom", "start", "end", "snps", "fst", "name")


plot_choi_100k <- ggplot(choi100k, aes((start+50000)/1e6, fst, col=chrom)) + 
    geom_point(size=0.1) + 
    facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
    theme_bw() + 
    theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
    labs(title="B", x="Genomic position (Mb)", y="Fst")+
    scale_colour_manual(values=colours)



plot_jm100k + plot_choi_100k + plot_layout(ncol=1, height = c(3,1))




colours <- c("#d0e11c", "#73d056", "#2db27d", "#21918c", "#2e6e8e", "#46327e")


strains100k <- read.table("tc_strains_poolseqfst.csv", header=T)
strains100k <- strains100k %>% filter(grepl("chr_[12345X]", chrom)) %>% arrange(chrom, start)

#MTci1_pool_adultMF_post-BZ.1:MTci2_pool_adultMF.1
#MTci2_pool_L4.1:MTci2_pool_adultMF.1
#MTci2_pool_adultMF.1:MTci5_pool_L4_MOTRI.1
#MTci2_pool_adultMF.1:MTci5_pool_adultMF_post-BZ.1
#MTci2_pool_adultMF.1:MTci5_pool_adultM_post-IVM.1
#MTci2_pool_adultMF.1:MTci7_pool_adultMF_MOX-R.1

strains100k_t1_t2 <- strains100k  %>% select(chrom, start, end, snps, MTci1_pool_adultMF_post.BZ.1.MTci2_pool_adultMF.1)
strains100k_t1_t2$name <- "1.Tci1(BZ)_v_Tci2"
colnames(strains100k_t1_t2) <- c("chrom", "start", "end", "snps", "fst", "name")

strains100k_t5_t2 <- strains100k  %>% select(chrom, start, end, snps, MTci2_pool_adultMF.1.MTci5_pool_L4_MOTRI.1)
strains100k_t5_t2$name <- "2.Tci5_v_Tci2"
colnames(strains100k_t5_t2) <- c("chrom", "start", "end", "snps", "fst", "name")

strains100k_t5bz_t2 <- strains100k  %>% select(chrom, start, end, snps, MTci2_pool_adultMF.1.MTci5_pool_adultMF_post.BZ.1)
strains100k_t5bz_t2$name <- "3.Tci5(BZ)_v_Tci2"
colnames(strains100k_t5bz_t2) <- c("chrom", "start", "end", "snps", "fst", "name")

strains100k_t5ivm_t2 <- strains100k %>% select(chrom, start, end, snps, MTci2_pool_adultMF.1.MTci5_pool_adultM_post.IVM.1)
strains100k_t5ivm_t2$name <- "4.Tci5(IVM)_v_Tci2"
colnames(strains100k_t5ivm_t2) <- c("chrom", "start", "end", "snps", "fst", "name")

strains100k_t7_t2 <- strains100k %>% select(chrom, start, end, snps, MTci2_pool_adultMF.1.MTci7_pool_adultMF_MOX.R.1)
strains100k_t7_t2$name <- "5.Tci7(MOX)_v_Tci2"
colnames(strains100k_t7_t2) <- c("chrom", "start", "end", "snps", "fst", "name")

data <- bind_rows(strains100k_t1_t2, strains100k_t5_t2, strains100k_t5bz_t2, strains100k_t5ivm_t2, strains100k_t7_t2)


data_gws <- data %>%
    group_by(name) %>%
    summarise(GWS = mean(fst) + 3*sd(fst))

plot_strains100k <- ggplot(data, aes((start+50000)/1e6, fst, col=chrom)) + 
    geom_hline(data = data_gws,  aes(yintercept = GWS),  linetype = "dashed",  col = "black") +
    geom_point(size=0.1) + 
    labs(title="") +
    facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
    theme_bw() + 
    theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
    labs(title="C", x="Genomic position (Mb)", y="Fst") +
    scale_colour_manual(values=colours)




plot_jm100k + plot_choi_100k + plot_strains100k + plot_layout(ncol=1, height = c(3,1,5))
```

strains100k_t7_t5 <- strains100k %>% select(chrom, start, end, snps, MTci5_pool_adultM_post.IVM.1.MTci7_pool_adultMF_MOX.R.1)
strains100k_t7_t5$name <- "Tci7_v_Tci5"
colnames(strains100k_t7_t5) <- c("chrom", "start", "end", "snps", "fst", "name")



ggplot(strains100k_t7_t5, aes((start+50000)/1e6, fst, col=chrom)) + 
    geom_point(size=0.1) + 
    labs(title="") +
    facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
    theme_bw() + 
    theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
    labs(title="C", x="Genomic position (Mb)", y="Fst") +
    scale_colour_manual(values=colours)











```bash
module load grenedalf/0.3.0

# choi data
bsub.py 1 grendalf_choi_diversity \
"grenedalf diversity \
--reference-genome-fasta-file teladorsagia_circumcincta_tci2_wsi3.0.genome.fa \
--allow-file-overwriting \
--measure all \
--file-prefix tc_choi_poolseq \
--sam-min-map-qual 30 \
--sam-min-base-qual 30 \
--filter-sample-min-count 2 \
--filter-sample-min-coverage 10 \
--filter-sample-max-coverage 100 \
--pool-sizes 800 \
--window-type sliding \
--window-sliding-width 100000 \
--popoolation-format \
--sam-path RS3.bam \
--sam-path SINBRED.bam"


# farm data
bsub.py 1 grendalf_farm_diversity \
"grenedalf diversity \
--reference-genome-fasta-file teladorsagia_circumcincta_tci2_wsi3.0.genome.fa \
--allow-file-overwriting \
--measure all \
--file-prefix tc_jm_farm_poolseq \
--sam-min-map-qual 30 \
--sam-min-base-qual 30 \
--filter-sample-min-count 2 \
--filter-sample-min-coverage 10 \
--filter-sample-max-coverage 100 \
--pool-sizes 182 \
--window-type sliding \
--window-sliding-width 100000 \
--popoolation-format \
--sam-path F2_POST.bam \
--sam-path F2_PRE.bam \
--sam-path F3_POST_A.bam \
--sam-path F3_POST_B.bam \
--sam-path F3_PRE_A.bam \
--sam-path F3_PRE_B.bam"


# strain data
bsub.py 1 grendalf_strains_diversity \
"grenedalf diversity \
--reference-genome-fasta-file teladorsagia_circumcincta_tci2_wsi3.0.genome.fa \
--allow-file-overwriting \
--measure all \
--file-prefix tc_strains_poolseq \
--sam-min-map-qual 30 \
--sam-min-base-qual 30 \
--filter-sample-min-count 2 \
--filter-sample-min-coverage 10 \
--filter-sample-max-coverage 100 \
--pool-sizes 200 \
--window-type sliding \
--window-sliding-width 100000 \
--popoolation-format \
--sam-path MTci1_pool_adultMF_post-BZ.bam \
--sam-path MTci2_pool_adultMF.bam \
--sam-path MTci2_pool_L4.bam \
--sam-path MTci5_pool_adultMF_post-BZ.bam \
--sam-path MTci5_pool_adultM_post-IVM.bam \
--sam-path MTci5_pool_L4_MOTRI.bam \
--sam-path MTci7_pool_adultMF_MOX-R.bam"
```








### TESTING - CNV analysis
- https://github.com/andrewkern/poolDiffCNV

```bash

# step 1

samtools view RS3.bam | python2.7 findEvertedInserts.py > everted_inserts_poolA.tsv
samtools view SINBRED.bam | python2.7 findEvertedInserts.py > everted_inserts_poolB.tsv

samtools view RS3.bam | python2.7 findDistantInserts.py 100000 > distant_inserts_poolA.tsv &
samtools view SINBRED.bam | python2.7 findDistantInserts.py 100000 > distant_inserts_poolB.tsv &

# step 2

#--- duplications
cat everted_inserts_poolA.tsv | python2.7  clusterEvertedInserts.py 100000 > clustered_everted_inserts_poolA.tsv &
cat everted_inserts_poolB.tsv | python2.7  clusterEvertedInserts.py 100000 > clustered_everted_inserts_poolB.tsv &

#--- deletions
cat distant_inserts_poolA.tsv | python2.7 clusterDistantInserts.py 2000 5 100 > clustered_distant_inserts_poolA.tsv &
cat distant_inserts_poolB.tsv | python2.7 clusterDistantInserts.py 2000 5 100 > clustered_distant_inserts_poolB.tsv &



```






## Benzimdazole analyses
```bash
# working dir:
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/ANALYSIS

# file containing position of variants
> btubulin_variant_positions.txt
# edit this to include the three isotype 1 variants 

# extract allele count data for each farm
vcftools --gzvcf TCIRC.raw.vcf.gz \
    --positions btubulin_variant_positions.txt \
    --extract-FORMAT-info AD \
    --out btubulin_variant_positions

# convert allele count data to variant frequency
grep "^tci2" btubulin_variant_positions.AD.FORMAT |\
    awk -F '[\t, ]' '{print $1, $2, $4/($3+$4), $6/($5+$6), $8/($7+$8), $10/($9+$10), $12/($11+$12), $14/($13+$14), $16/($15+$16), $18/($17+$18), $20/($19+$20), $22/($21+$22), $24/($23+$24), $26/($25+$26), $28/($27+$28), $30/($29+$30), $32/($31+$32)}' OFS="\t" > btubulin_variant_positions.AD.freq

```

```R
# load libraries
library(tidyverse)
library(reshape2)

library(rstatix)

# reformat data
btub1 <- read.table("btubulin_variant_positions.AD.freq")

colnames(btub1) <- c("CHR", "POS", "Farm 1 (post)", "Farm 1 (pre)", "Farm 2 (post)", "F3_POST_B", "Farm 2 (pre)", "F3_PRE_B", "MTci1", 	"MTci2", "MTci2_pool_L4", "MTci5 (post-BZ)",	"MTci5 (post-IVM)",	"MTci5", "MTci7", "RS3", "SINBRED")

btub1 <- btub1 %>% 
            mutate(POS = str_replace(POS,  c("62291012", "62291019"),  c("E198L", "F200Y")))

btub1 <- melt(btub1,  id = c("CHR",  "POS"),  variable.name = "SAMPLE_ID")
colnames(btub1) <- c("CHR", "POS", "SAMPLE_ID", "ALLELE_FREQ")

groups <- rep(c("FARM", "STRAIN", "CHOI"), c(12,14,4))

groups <- melt(groups)
colnames(groups) <- "group"

data <- data.frame(btub1, groups)

data <- data %>% filter(!grepl("_B", SAMPLE_ID)) %>% filter(!grepl("_pool_L4", SAMPLE_ID))

# make the figure
ggplot(data, aes(x = SAMPLE_ID, y = ALLELE_FREQ, fill = factor(POS))) +
     geom_bar(position = "dodge",  stat = "identity") +
     labs(title = "A",  x="Sampling location",  y="Resistant allele frequency",  fill = "Variant") +
     theme_bw() + 
     theme(text = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust = 1)) + 
     facet_grid(. ~ group, space="free_x", scales="free_x") +
     ylim(0, 1)

# save it
ggsave("FigureSX_USfarm_btub1.pdf",  useDingbats=FALSE, width=170, height=100, units="mm")
ggsave("FigureSX_USfarm_btub1.png")

```







## make a plot comparing Tcirc and Haem Fst values
- use Tcirc fst data (above)
- remap haemonchus XQTL data
- generate 1-to-1 ortholog comparison between Tcirc and Haem




### Mapping Haemonchus data 
- datasets
XQTL_F3_L3_n200_IVM_pre_01	21395_2_1
XQTL_F3_L3_n200_IVM_post_01	21395_2_2
XQTL_F3_L3_n200_IVM_pre_02	23241_2_1
XQTL_F3_L3_n200_IVM_post_02	23241_2_2
XQTL_F3_L3_n200_IVM_pre_03	23241_7_1
XQTL_F3_L3_n200_IVM_post_03	23241_7_2


```bash
# Download the data from iRODs using irods_extractor
cd ~/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/DATA

irods_extractor --studyid 2679  --runid 21395 --laneid 2
irods_extractor --studyid 2679  --runid 23241 --laneid 2
irods_extractor --studyid 2679  --runid 23241 --laneid 7

# Mapping to Haemonchus reference genome
module load mapping-helminth/v1.0.9

mapping-helminth --input hc_xqtl_sample_manifect.txt --reference HAEM_V4_final.chr.fa

rm -rf work

# nice way to find all the bams, and them make symbolic links to a target directory
OUTPUT_DIR=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/ANALYSIS
find ~+ -type f -name 'XQTL*.bam*' -exec ln -vs "{}" $OUTPUT_DIR/ ';'

cd /lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/ANALYSIS

# get reference - same version as in WBPS18
ln -s ../../../../haemonchus_contortus/GENOME/REF/HAEM_V4_final.chr.fa

# run Fst analysis
#--- note - running in 10 kb, rather than 100kb in Tc, as Hc genome smaller and reflects what was done in the original analysis
#--- the IVM QTL is also much smaller in Hc, so would get lost at 100 kb window

module load grenedalf/0.3.0

bsub.py 1 grendalf_hc_xqtl_ivm_fst \
"grenedalf fst \
--reference-genome-fasta-file HAEM_V4_final.chr.fa \
--allow-file-overwriting \
--method unbiased-nei \
--file-prefix hc_xqtl_ivm_poolseq \
--sam-min-map-qual 30 \
--sam-min-base-qual 30 \
--filter-sample-min-count 2 \
--filter-sample-min-coverage 10 \
--filter-sample-max-coverage 200 \
--pool-sizes 400 \
--window-type sliding \
--window-sliding-width 10000 \
--write-pi-tables \
--separator-char tab \
--sam-path XQTL_F3_L3_n200_IVM_pre_01.bam \
--sam-path XQTL_F3_L3_n200_IVM_post_01.bam \
--sam-path XQTL_F3_L3_n200_IVM_pre_02.bam \
--sam-path XQTL_F3_L3_n200_IVM_post_02.bam \
--sam-path XQTL_F3_L3_n200_IVM_pre_03.bam \
--sam-path XQTL_F3_L3_n200_IVM_post_03.bam"
```







### Orthofinder analysis to ID 1-to-1 orthologs
- using the previous orthofinder run data used to compare genomes

```bash
cd ~/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/ANALYSIS/ORTHOFINDER_CHR5

ln -s /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/GENOME_COMPARISON/ORTHOFINDER/proteins/OrthoFinder/Results_Feb22/Orthologues/Orthologues_teladorsagia_circumcincta_tci2_wsi3.0.proteins.longest-isoform/teladorsagia_circumcincta_tci2_wsi3.0.proteins.longest-isoform__v__hc.unique.proteins.tsv

cat teladorsagia_circumcincta_tci2_wsi3.0.proteins.longest-isoform__v__hc.unique.proteins.tsv | awk '{if(NF==3) print $2,$3}' OFS="\t" > hc_v_tc_1-to-1-orthologs.list

wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3.gz
gunzip haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3

mv haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3 hc.gff3

ln -s /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/V3/teladorsagia_circumcincta_tci2_wsi3.0.annotation.gff3 tc.gff3

# some strange hidden characters need to be removed
sed -i $'s/[^[:print:]\t]//g' hc_v_tc_1-to-1-orthologs.list

# get gene coordinates
awk -F '[\t;]' '$3=="mRNA" {print $1, $4, $5, $7, $9}' OFS="\t" tc.gff3 | sed 's/ID=//g' > tc.gff.coords
awk -F '[\t;]' '$3=="mRNA" {print $1, $4, $5, $7, $9}' OFS="\t" hc.gff3 | sed 's/ID=transcript://g' > hc.gff.coords

# extract positional data for both Hc and Ce
cat hc_v_tc_1-to-1-orthologs.list | while read -r tcirc haem ; do 
    COORDS1=$(grep "${tcirc}" tc.gff.coords) ; 
    COORDS2=$(grep "${haem}" hc.gff.coords) ; 
    echo -e "${COORDS1}\t${COORDS2}" ; 
    done > hc_v_tc_1-to-1.coords


cat hc_v_tc_1-to-1.coords | grep "tci2_wsi3.0_chr_5" | grep "hcontortus_chr5_Celeg_TT_arrow_pilon" > hc_v_tc_1-to-1.chr5.coords

wc -l hc_v_tc_1-to-1.chr5.coords
#> 1253 hc_v_tc_1-to-1.chr5.coords
```




### Making the figure
working dir: ~/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/POOLSEQ/ANALYSIS

```R
library(tidyverse)

# strains
rawdata_strains <- read.table("tc_strains_poolseqfst.csv", header=T)

data_strains <- rawdata_strains %>% 
mutate(rawdata_strains, mean_fst  = rowMeans(select(.,MTci2_pool_adultMF.1.MTci5_pool_L4_MOTRI.1, MTci2_pool_adultMF.1.MTci5_pool_adultMF_post.BZ.1, MTci2_pool_adultMF.1.MTci5_pool_adultM_post.IVM.1, MTci2_pool_adultMF.1.MTci7_pool_adultMF_MOX.R.1)))

data_strains_chr5 <- data_strains %>% filter(grepl("chr_5", chrom)) %>% arrange(chrom, start)

data_strains_chr5$group <- "STRAIN"

data_strains_chr5 <- data_strains_chr5 %>% select(chrom, start, end, snps, mean_fst, group)

data_strains_chr5 <- data_strains_chr5 %>% mutate(., mean_fst_pc = mean_fst/max(mean_fst))

# choi
rawdata_choi <- read.table("tc_choi_poolseqfst.csv", header=T)

rawdata_choi <- rawdata_choi %>% 
mutate(rawdata_choi, mean_fst = rowMeans(select(., RS3.1.SINBRED.1), na.rm = TRUE)) 

data_choi_chr5 <- rawdata_choi %>% filter(grepl("chr_5", chrom)) %>% arrange(chrom, start)

data_choi_chr5$group <- "CHOI"

data_choi_chr5 <- data_choi_chr5 %>% select(chrom, start, end, snps, mean_fst, group)

data_choi_chr5 <- data_choi_chr5 %>% mutate(., mean_fst_pc = mean_fst/max(mean_fst))



# farms
rawdata_farms <- read.table("tc_jm_farm_poolseqfst.csv", header=T)

data_farms <- rawdata_farms %>% 
mutate(rawdata_farms, mean_fst = rowMeans(select(., F2_POST.1.F2_PRE.1, F2_POST.1.F3_POST_A.1, F3_POST_A.1.F3_PRE_A.1), na.rm = TRUE))

data_farms_chr5 <- data_farms %>% filter(grepl("chr_5", chrom)) %>% arrange(chrom, start)

data_farms_chr5$group <- "FARM"

data_farms_chr5 <- data_farms_chr5 %>% select(chrom, start, end, snps, mean_fst, group)

data_farms_chr5 <- data_farms_chr5 %>% mutate(., mean_fst_pc = mean_fst/max(mean_fst))


data <- bind_rows(data_strains_chr5, data_choi_chr5, data_farms_chr5)


# plot - top panel
plot_tc_fst_heatmap <- 
    ggplot(data, aes(start/1e6, group, fill=mean_fst_pc*100)) + 
    geom_tile() + 
    scale_fill_gradient(low = "white", high = "blue", na.value="white") +
    theme_minimal() +
    scale_x_continuous(n.breaks = 10, limits=c(0,91)) +
    labs(fill="% max Fst") +
    theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank()) 




# 1-to-1 orthologs on chromosome 5 - middle panel
ortho_data <- read.table("ORTHOFINDER_CHR5/hc_v_tc_1-to-1.chr5.coords", header=F)

ortho_data <- ortho_data %>% mutate(., col = ifelse(V7 > 37200000 & V7 < 37700000, "Hc_IVM_QTL", "Other"))


plot_orthologs <- 
    ggplot(ortho_data) +
    geom_segment(aes(x = V2/1e6, y = 1, xend = V7/1e6, yend = 0, col=col), alpha=0.2) +
    geom_segment(data = subset(ortho_data, col=="Hc_IVM_QTL"), aes(x = V2/1e6, y = 1, xend = V7/1e6, yend = 0, col=col)) +
    scale_color_manual(values = c("Other"="grey", "Hc_IVM_QTL"="red")) + 
    theme_minimal() +
    labs(col="1-to-1 othologs") +
    scale_x_continuous(n.breaks = 10, limits=c(0,91)) +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) 





# Haemonchus data - bottom panel

# for testing - used previosuly mapped data
# ln -s ~sd21/lustre_link/haemonchus_contortus/XQTL/05_ANALYSIS/IVM/XQTL_IVM.merged.fst


rawdata_hc <- read.table("hc_xqtl_ivm_poolseqfst.csv", header = T)

data_hc <- rawdata_hc %>% 
mutate(rawdata_hc, mean_fst  = rowMeans(select(.,XQTL_F3_L3_n200_IVM_post_01.1.XQTL_F3_L3_n200_IVM_pre_01.1, XQTL_F3_L3_n200_IVM_post_02.1.XQTL_F3_L3_n200_IVM_pre_02.1, XQTL_F3_L3_n200_IVM_post_03.1.XQTL_F3_L3_n200_IVM_pre_03.1)))

data_hc_chr5 <- data_hc %>% filter(grepl("chr5", chrom)) %>% arrange(chrom, start)

data_hc_chr5 <- data_hc_chr5 %>% mutate(., mean_fst_pc = mean_fst/max(mean_fst, na.rm=TRUE))

plot_hc_fst_heatmap <- 
    ggplot(data_hc_chr5, aes(start/1e6, 'Haem_XQTL', fill=mean_fst_pc*100)) +    
    geom_tile() + 
    scale_fill_gradient(low = "white", high = "blue", na.value="white") + 
    labs(x="Genome position (Mb)", fill="% max Fst") +
    theme_minimal() +
    scale_x_continuous(n.breaks = 10, limits=c(0,91)) +
    theme(axis.title.y=element_blank())




# bring it all together
plot_tc_fst_heatmap + 
    plot_orthologs + 
    plot_hc_fst_heatmap + 
    plot_layout(ncol=1, heights = c(3,3,1), guides = "collect")  
```