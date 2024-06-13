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
F2_POST,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Farm2_Post_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Farm2_Post_2.fq.gz
F2_PRE,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Farm2_Pre_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Farm2_Pre_2.fq.gz
F3_POST_A,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Farm3_PostA_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Farm3_PostA_2.fq.gz
F3_POST_B,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Farm3_PostB_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Farm3_PostB_2.fq.gz
F3_PRE_A,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Farm3_PreA_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Farm3_PreA_2.fq.gz
F3_PRE_B,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Farm3_PreB_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Farm3_PreB_2.fq.gz
RS3,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/RS3_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/RS3_2.fq.gz
SINBRED,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Sinbred_1_1.fq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/DATA/Sinbred_1_2.fq.gz
s
```




OUTPUT_DIR=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/ANALYSIS
find ~+ -type f -name '*.bam*' -exec ln -vs "{}" $OUTPUT_DIR/ ';'


cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ/ANALYSIS




bsub.py --queue long 10 grenedalf_div \
"grenedalf diversity \
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
--sam-path F3_PRE_B.bam \
--sam-path RS3.bam \
--sam-path SINBRED.bam"



RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-SINBRED:1-theta-watterson.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-SINBRED:1-theta-pi.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-SINBRED:1-tajimas-d.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-RS3:1-theta-watterson.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-RS3:1-theta-pi.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-RS3:1-tajimas-d.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_PRE_B:1-theta-watterson.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_PRE_B:1-theta-pi.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_PRE_B:1-tajimas-d.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_PRE_A:1-theta-watterson.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_PRE_A:1-theta-pi.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_PRE_A:1-tajimas-d.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_POST_B:1-theta-watterson.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_POST_B:1-theta-pi.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_POST_B:1-tajimas-d.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_POST_A:1-theta-watterson.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_POST_A:1-theta-pi.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_POST_A:1-tajimas-d.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F2_PRE:1-theta-watterson.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F2_PRE:1-theta-pi.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F2_PRE:1-tajimas-d.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F2_POST:1-theta-watterson.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F2_POST:1-theta-pi.csv
RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F2_POST:1-tajimas-d.csv



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









data_preA <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_PRE_A:1-theta-watterson.csv", header=F)
data_preA_chr <- data_preA %>% filter(grepl("chr", V1))
data_preA_chr_5 <- data_preA %>% filter(grepl("chr_5", V1))

data_postA <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_POST_A:1-theta-watterson.csv", header=F)
data_postA_chr <- data_postA %>% filter(grepl("chr", V1))
data_postA_chr_5 <- data_postA %>% filter(grepl("chr_5", V1))

plot_A_chr <- ggplot() + geom_point(aes(1:nrow(data_preA_chr), data_preA_chr$V5/data_postA_chr$V5, col=data_preA_chr$V1), size=0.5) + ylim(0,15)

plot_A_chr_5 <- ggplot() + geom_point(aes(data_preA_chr_5$V2, data_preA_chr_5$V5/data_postA_chr_5$V5, col=data_preA_chr_5$V1), size=0.5)



data_preB <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_PRE_B:1-theta-watterson.csv", header=F)
data_preB_chr <- data_preB %>% filter(grepl("chr", V1)) 
data_preB_chr_5 <- data_preB %>% filter(grepl("chr_5", V1))

data_postB <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_POST_B:1-theta-watterson.csv", header=F)
data_postB_chr <- data_postB %>% filter(grepl("chr", V1))
data_postB_chr_5 <- data_postB %>% filter(grepl("chr_5", V1))

plot_B_chr <- ggplot() + geom_point(aes(1:nrow(data_preB_chr), data_preB_chr$V5/data_postB_chr$V5, col=data_preB_chr$V1), size=0.5) + ylim(0,15)

plot_B_chr_5 <- ggplot() + geom_point(aes(data_preB_chr_5$V2, data_preB_chr_5$V5/data_postB_chr_5$V5, col=data_preB_chr_5$V1), size=0.5)

plot_A_chr + plot_B_chr + plot_layout(ncol=1, guides = "collect")

plot_A_chr_5 + plot_B_chr_5 + plot_layout(ncol=1, guides = "collect")






data_preA <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_PRE_A:1-tajimas-d.csv", header=F)
data_preA_chr <- data_preA %>% filter(grepl("chr", V1))
data_preA_chr_5 <- data_preA %>% filter(grepl("chr_5", V1))

data_postA <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_POST_A:1-tajimas-d.csv", header=F)
data_postA_chr <- data_postA %>% filter(grepl("chr", V1))
data_postA_chr_5 <- data_postA %>% filter(grepl("chr_5", V1))

plot_A_chr <- ggplot() + geom_point(aes(1:nrow(data_preA_chr), data_preA_chr$V5/data_postA_chr$V5, col=data_preA_chr$V1), size=0.5) + ylim(-5,5)

plot_A_chr_5 <- ggplot() + geom_point(aes(data_preA_chr_5$V2, data_preA_chr_5$V5/data_postA_chr_5$V5, col=data_preA_chr_5$V1), size=0.5) + ylim(-5,5)



data_preB <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_PRE_B:1-tajimas-d.csv", header=F)
data_preB_chr <- data_preB %>% filter(grepl("chr", V1)) 
data_preB_chr_5 <- data_preB %>% filter(grepl("chr_5", V1))

data_postB <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_POST_B:1-tajimas-d.csv", header=F)
data_postB_chr <- data_postB %>% filter(grepl("chr", V1))
data_postB_chr_5 <- data_postB %>% filter(grepl("chr_5", V1))

plot_B_chr <- ggplot() + geom_point(aes(1:nrow(data_preB_chr), data_preB_chr$V5/data_postB_chr$V5, col=data_preB_chr$V1), size=0.5) + ylim(-5,5)

plot_B_chr_5 <- ggplot() + geom_point(aes(data_preB_chr_5$V2, data_preB_chr_5$V5/data_postB_chr_5$V5, col=data_preB_chr_5$V1), size=0.5) + ylim(-5,5)

plot_A_chr + plot_B_chr + plot_layout(ncol=1, guides = "collect")

plot_A_chr_5 + plot_B_chr_5 + plot_layout(ncol=1, guides = "collect")





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

plot_A_chr_X <- ggplot() + geom_point(aes(data_preA_chr_X$V2, data_preA_chr_X$V5/data_postA_chr_X$V5, col=data_preA_chr_X$V1), size=0.5)