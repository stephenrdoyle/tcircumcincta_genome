

## mapping
```bash
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/POOLSEQ

module load mapping-helminth/v1.0.8

bsub.py --queue long 10 mapping "mapping-helminth --input sample_manifest.txt --reference teladorsagia_circumcincta_tci2_wsi2.4.fa"

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
--file-prefix RS3_v_Sinbred.bwa.freq0.1.cov.10.div \
--popoolation-format \
--sam-min-map-qual 30 \
--sam-min-base-qual 30 \
--pool-sizes 400 \
--window-type sliding \
--window-sliding-width 10000 \
--filter-sample-min-count 2 \
--reference-genome-fasta-file teladorsagia_circumcincta_tci2_wsi2.4.fa \
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

data_preA <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_PRE_A:1-theta-pi.csv", header=F)
data_preA_chr <- data_preA %>% filter(grepl("chr", V1))
data_preA_notchr <- data_preA %>% filter(!grepl("chr", V1))
data_preA_chr_5 <- data_preA %>% filter(grepl("chr_5", V1))

data_postA <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_POST_A:1-theta-pi.csv", header=F)
data_postA_chr <- data_postA %>% filter(grepl("chr", V1))
data_postA_notchr <- data_postA %>% filter(!grepl("chr", V1))
data_postA_chr_5 <- data_postA %>% filter(grepl("chr_5", V1))

plot_A_chr <- ggplot() + geom_point(aes(1:nrow(data_preA_chr), data_preA_chr$V5/data_postA_chr$V5, col=data_preA_chr$V1), size=0.5) + ylim(0,15)

plot_A_notchr <- ggplot() + geom_point(aes(1:nrow(data_preA_notchr), data_preA_notchr$V5/data_postA_notchr$V5, col=data_preA_notchr$V1), size=0.5) + ylim(0,15) + theme(legend.position="none")

plot_A_chr_5 <- ggplot() + geom_point(aes(data_preA_chr_5$V2, data_preA_chr_5$V5/data_postA_chr_5$V5, col=data_preA_chr_5$V1), size=0.5)



data_preB <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_PRE_B:1-theta-pi.csv", header=F)
data_preB_chr <- data_preB %>% filter(grepl("chr", V1)) 
data_preB_notchr <- data_preB %>% filter(!grepl("chr", V1)) 
data_preB_chr_5 <- data_preB %>% filter(grepl("chr_5", V1))

data_postB <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-F3_POST_B:1-theta-pi.csv", header=F)
data_postB_chr <- data_postB %>% filter(grepl("chr", V1))
data_postB_notchr <- data_postB %>% filter(!grepl("chr", V1)) 
data_postB_chr_5 <- data_postB %>% filter(grepl("chr_5", V1))

plot_B_chr <- ggplot() + geom_point(aes(1:nrow(data_preB_chr), data_preB_chr$V5/data_postB_chr$V5, col=data_preB_chr$V1), size=0.5) + ylim(0,15)

plot_B_notchr <- ggplot() + geom_point(aes(1:nrow(data_preB_notchr), data_preB_notchr$V5/data_postB_notchr$V5, col=data_preB_notchr$V1), size=0.5) + ylim(0,15) + theme(legend.position="none")

plot_B_chr_5 <- ggplot() + geom_point(aes(data_preB_chr_5$V2, data_preB_chr_5$V5/data_postB_chr_5$V5, col=data_preB_chr_5$V1), size=0.5)

plot_A_chr + plot_B_chr + plot_layout(ncol=1, guides = "collect")

plot_A_notchr + plot_B_notchr + plot_layout(ncol=1, guides = "collect")

plot_A_chr_5 + plot_B_chr_5 + plot_layout(ncol=1, guides = "collect")


plot_AB_pre_chr <- ggplot() + geom_point(aes(1:nrow(data_preA_chr), data_preA_chr$V5/data_preB_chr$V5, col=data_preA_chr$V1), size=0.5)

plot_AB_post_chr <- ggplot() + geom_point(aes(1:nrow(data_postA_chr), data_postA_chr$V5/data_postB_chr$V5, col=data_postA_chr$V1), size=0.5)


plot_AB_pre_chr + plot_AB_post_chr + plot_layout(ncol=1, guides = "collect")









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
library(patchwork)

data_preA <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-SINBRED:1-theta-pi.csv", header=F)
data_preA_chr <- data_preA %>% filter(grepl("chr", V1)) %>% arrange(V1,V2)
data_preA_notchr <- data_preA %>% filter(!grepl("chr", V1)) %>% arrange(V1,V2)
data_preA_chr_5 <- data_preA %>% filter(grepl("chr_5", V1)) %>% arrange(V1,V2)
data_preA_chr_X <- data_preA %>% filter(grepl("chr_X", V1)) %>% arrange(V1,V2)

data_postA <- read.table("RS3_v_Sinbred.bwa.freq0.1.cov.10.divdiversity-RS3:1-theta-pi.csv", header=F)
data_postA_chr <- data_postA %>% filter(grepl("chr", V1)) %>% arrange(V1,V2)
data_postA_notchr <- data_postA %>% filter(!grepl("chr", V1)) %>% arrange(V1,V2)
data_postA_chr_5 <- data_postA %>% filter(grepl("chr_5", V1)) %>% arrange(V1,V2)
data_postA_chr_X <- data_postA %>% filter(grepl("chr_X", V1)) %>% arrange(V1,V2)

plot_A_chr <- ggplot() + geom_point(aes(1:nrow(data_preA_chr), data_preA_chr$V5/data_postA_chr$V5, col=data_preA_chr$V1), size=0.5) + ylim(0,15)

plot_A_notchr <- ggplot() + geom_point(aes(1:nrow(data_preA_notchr), data_preA_notchr$V5/data_postA_notchr$V5, col=data_preA_notchr$V1), size=0.5) + ylim(0,15) + theme(legend.position="none")

plot_A_chr_5 <- ggplot() + geom_point(aes(data_preA_chr_5$V2, data_preA_chr_5$V5/data_postA_chr_5$V5, col=data_preA_chr_5$V1), size=0.5)

plot_A_chr_X <- ggplot() + geom_point(aes(data_preA_chr_X$V2, data_preA_chr_X$V5/data_postA_chr_X$V5, col=data_preA_chr_X$V1), size=0.5)