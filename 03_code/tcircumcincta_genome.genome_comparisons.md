# Genome comparisons



## Cumulative genome size plots
- want to compare the three teladorsagia genomes

```bash 
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/GENOME_COMPARISON

ln -s ../REFERENCES/teladorsagia_circumcincta_tci2_wsi2.4.fa tci2_wsi2.4.fa
ln -s ../REFERENCES/DNAZOO/Tcircumcincta.DNAzoo.fa DNAZOO.fa
ln -s ../REFERENCES/WASHU/teladorsagia_circumcincta.PRJNA72569.WBPS17.fa WASHU.fa



# make an index file to get data from
for i in *fa; do
    samtools faidx ${i};
    done


for i in WASHU.fa.fai DNAZOO.fa.fai tci2_wsi2.4.fa.fai; do
    name=$(echo ${i%.fa.fai})
    cat ${i} | sort -k2nr | awk -v name=${name} '{print name,$1,$2,total += $2; $0 = total}' OFS="\t";
done > cumulative_genomes_data.txt
```

```R
# libraries
library(tidyverse)
library(ggsci)
library(viridis)


data1 <- read.table("cumulative_genomes_data.txt", sep="\t", header=F)
data2 <- data1 %>% group_by(V1) %>% summarise(V3 = as.numeric(max(V3))) %>% mutate(V2="start") %>% mutate(V4=1) %>% relocate(V1,V2,V3,V4)
data <- bind_rows(data1,data2)
max <- data1 %>% group_by(V1) %>% summarise(V3 = as.numeric(max(V4)))

ggplot(data) + 
    geom_vline(data=max,aes(xintercept=V3/1e6, col=V1), linetype="dotted", linewidth=0.75) + 
    geom_line(aes(V4/1e6,V3/1e6,col=V1),linewidth=1) +
    geom_point(aes(V4/1e6,V3/1e6,col=V1),size=2) + 
    theme_bw() + theme(legend.position="bottom") +
    scale_color_npg() +
    ylim(0,100) +
    labs(x="Cumulative genome size (Mb)", y= "Scaffold/contig length (Mb)", colour="Genome assembly")

ggsave("cumulative_genome_comparison.pdf", height=5, width=5, units="in")
ggsave("cumulative_genome_comparison.png")
```
![](../04_analysis/cumulative_genome_comparison.png)


```R 
# libraries
library(tidyverse)
library(treemapify)

data1 <- read.table("cumulative_genomes_data.txt", sep="\t", header=F)
data2 <- data1 %>% group_by(V1) %>% summarise(V3 = as.numeric(max(V3))) %>% mutate(V2="start") %>% mutate(V4=1) %>% relocate(V1,V2,V3,V4)
data <- bind_rows(data1,data2)

ggplot(data, aes(area=V3, fill=V3)) + geom_treemap() + facet_grid(.~V1)


data <- read.table("cumulative_genomes_data.txt", sep="\t", header=F)
data <- data %>% group_by(V1) %>% top_n(n = 5000, wt = V3)
ggplot(data, aes(area=V3, fill=V3)) + geom_treemap() + facet_grid(.~V1)


data2 <- data %>% filter(V1=="tci2_wsi1.0")
ggplot(data2, aes(area=V3, fill=V3)) + geom_treemap() + facet_grid(.~V1)
```



## Comparison between tci2_2.4 and other genomes

```bash
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa.gz
gunzip haemonchus_contortus.PRJEB506.WBPS18.genomic.fa.gz


>hc_chr.fa
grep ">" haemonchus_contortus.PRJEB506.WBPS18.genomic.fa | grep "chr" | sed 's/>//'g | cut -f1 -d " " | while read -r NAME; do 
    samtools faidx haemonchus_contortus.PRJEB506.WBPS18.genomic.fa ${NAME} >> hc_chr.fa; 
    done



>tc_2.4_chr.fa
grep ">" tci2_wsi2.4.fa | grep "chr" | grep -v "mtDNA" | sed 's/>//'g | while read -r NAME; do 
    samtools faidx tci2_wsi2.4.fa ${NAME} >> tc_2.4_chr.fa; 
    done


>dnazoo_chr.fa
for NAME in HiC_scaffold_1 HiC_scaffold_2 HiC_scaffold_3 HiC_scaffold_4 HiC_scaffold_5 HiC_scaffold_6; do 
    samtools faidx DNAZOO.fa ${NAME} >> dnazoo_chr.fa; 
    done

# map tc_2.4 to DNAZOO
minimap2 -x asm20 tc_2.4_chr.fa dnazoo_chr.fa  > tc2.4_v_dnazoo.paf

# map tc_2.4 to Hcon
minimap2 -x asm20 tc_2.4_chr.fa hc_chr.fa  > tc2.4_v_haem.paf

# clean up the paf file for plotting
for i in *.paf; do
cat ${i}| /nfs/users/nfs_s/sd21/lustre_link/software/GENOME_ASSEMBLY/minimap/utils/bin/layout | cut -f1-19 > ${i%.paf}.layout;
done


samtools faidx tc_2.4_chr.fa
samtools faidx dnazoo_chr.fa
samtools faidx hc_chr.fa

cat <(sort tc_2.4_chr.fa.fai) <(sort dnazoo_chr.fa.fai) | awk '{print $1, "1", $2}' OFS="\t" > tc_2.4_dnazoo_chromosome_lengths.txt
cat <(sort tc_2.4_chr.fa.fai) <(sort hc_chr.fa.fai) | awk '{print $1, "1", $2}' OFS="\t"  > tc_2.4_haem_chromosome_lengths.txt


```

```R
library(tidyverse)
library(circlize)
library(viridis)

# tc2.4 vs dnazoo
chr <- read.table("tc_2.4_dnazoo_chromosome_lengths.txt", header=F)
colnames(chr) <- c("chr", "start", "end")
data<-read.table("tc2.4_v_dnazoo.layout")
data<-data %>% filter(V18>25000)

dnazoo_links <- data %>% select(V8, V10, V11, V18)
colnames(dnazoo_links) <- c("chr", "start", "end", "value")
tc_links <- data %>% select(V13, V15, V16, V18)
colnames(tc_links) <- c("chr", "start", "end", "value")



grid.col = c(tci2_wsi2.0_chr_1 = "#d0e11c", tci2_wsi2.0_chr_2 = "#73d056", tci2_wsi2.0_chr_3 = "#2db27d", tci2_wsi2.0_chr_4 = "#21918c", tci2_wsi2.0_chr_5 = "#2e6e8e", tci2_wsi2.0_chr_X = "#46327e", HiC_scaffold_1 = "grey", HiC_scaffold_2 = "grey", HiC_scaffold_3 = "grey", HiC_scaffold_4 = "grey", HiC_scaffold_5 = "grey", HiC_scaffold_6 = "grey")

names <- c("tci2_wsi2.0_chr_1", "tci2_wsi2.0_chr_2", "tci2_wsi2.0_chr_3", "tci2_wsi2.0_chr_4", "tci2_wsi2.0_chr_5", "tci2_wsi2.0_chr_X", "HiC_scaffold_1", "HiC_scaffold_2", "HiC_scaffold_3", "HiC_scaffold_4", "HiC_scaffold_5", "HiC_scaffold_6")

colours <- c("#d0e11c", "#73d056", "#2db27d", "#21918c", "#2e6e8e", "#46327e", "grey", "grey", "grey", "grey", "grey", "grey")

names_colours <- data.frame(names, colours)


tc_links_col <- left_join(tc_links, names_colours, by=c("chr" = "names"))

pdf("tc2.4_dnazoo_circulize.pdf")
circos.clear()

circos.genomicInitialize(chr)
circos.track(ylim = c(0, 1), bg.col = grid.col, bg.border = NA, track.height = 0.075)
circos.genomicLink(dnazoo_links, tc_links_col, border = NA, col = tc_links_col$colours)

circos.clear()
dev.off()





# tc2.4 vs haem
chr <- read.table("tc_2.4_haem_chromosome_lengths.txt", header=F)
colnames(chr) <- c("chr", "start", "end")
data<-read.table("tc2.4_v_haem.layout")
data<-data %>% filter(V18>5000)

haem_links <- data %>% select(V8, V10, V11, V18)
colnames(dnazoo_links) <- c("chr", "start", "end", "value")
tc_links <- data %>% select(V13, V15, V16, V18)
colnames(tc_links) <- c("chr", "start", "end", "value")




grid.col = c(tci2_wsi2.0_chr_1 = "#d0e11c", tci2_wsi2.0_chr_2 = "#73d056", tci2_wsi2.0_chr_3 = "#2db27d", tci2_wsi2.0_chr_4 = "#21918c", tci2_wsi2.0_chr_5 = "#2e6e8e", tci2_wsi2.0_chr_X = "#46327e", hcontortus_chr1_Celeg_TT_arrow_pilon = "grey", hcontortus_chr2_Celeg_TT_arrow_pilon = "grey", hcontortus_chr3_Celeg_TT_arrow_pilon = "grey", hcontortus_chr4_Celeg_TT_arrow_pilon = "grey", hcontortus_chr5_Celeg_TT_arrow_pilon = "grey", hcontortus_chrX_Celeg_TT_arrow_pilon = "grey")

names <- c("tci2_wsi2.0_chr_1", "tci2_wsi2.0_chr_2", "tci2_wsi2.0_chr_3", "tci2_wsi2.0_chr_4", "tci2_wsi2.0_chr_5", "tci2_wsi2.0_chr_X", "hcontortus_chr1_Celeg_TT_arrow_pilon", "hcontortus_chr2_Celeg_TT_arrow_pilon", "hcontortus_chr3_Celeg_TT_arrow_pilon", "hcontortus_chr4_Celeg_TT_arrow_pilon", "hcontortus_chr5_Celeg_TT_arrow_pilon", "hcontortus_chrX_Celeg_TT_arrow_pilon")

colours <- c("#d0e11c", "#73d056", "#2db27d", "#21918c", "#2e6e8e", "#46327e", "grey", "grey", "grey", "grey", "grey", "grey")

names_colours <- data.frame(names, colours)


tc_links_col <- left_join(tc_links, names_colours, by=c("chr" = "names"))


pdf("tc2.4_haem_circulize.pdf")
circos.clear()

circos.genomicInitialize(chr)
circos.track(ylim = c(0, 1), bg.col = grid.col, bg.border = NA, track.height = 0.075)
circos.genomicLink(haem_links, tc_links_col, border = NA, col = tc_links_col$colours)

circos.clear()
dev.off()







```

bsub.py 50 promer "promer --mum --prefix tc2.4_v_hc tc_2.4_chr.fa hc_chr.fa"

show-coords -lTH -L 1000 tc2.4_v_hc.delta | awk '{print $15,$11,$3,$4,"+",$14,$10,$1,$2,$5,$5,60}' OFS="\t" > tc2.4_v_hc.promer.pseudo-paf


# clean up the paf file for plotting
cat tc2.4_v_hc.promer.pseudo-paf | /nfs/users/nfs_s/sd21/lustre_link/software/GENOME_ASSEMBLY/minimap/utils/bin/layout | cut -f1-19 > tc2.4_v_hc.promer.layout



```R
# tc2.4 vs haem
chr <- read.table("tc_2.4_haem_chromosome_lengths.txt", header=F)
colnames(chr) <- c("chr", "start", "end")
data<-read.table("tc2.4_v_hc.promer.layout")
data<-data %>% filter(V18>2000)

haem_links <- data %>% select(V8, V10, V11, V18)
colnames(haem_links) <- c("chr", "start", "end", "value")
tc_links <- data %>% select(V13, V15, V16, V18)
colnames(tc_links) <- c("chr", "start", "end", "value")


grid.col = c(tci2_wsi2.0_chr_1 = "#d0e11c", tci2_wsi2.0_chr_2 = "#73d056", tci2_wsi2.0_chr_3 = "#2db27d", tci2_wsi2.0_chr_4 = "#21918c", tci2_wsi2.0_chr_5 = "#2e6e8e", tci2_wsi2.0_chr_X = "#46327e", hcontortus_chr1_Celeg_TT_arrow_pilon = "grey", hcontortus_chr2_Celeg_TT_arrow_pilon = "grey", hcontortus_chr3_Celeg_TT_arrow_pilon = "grey", hcontortus_chr4_Celeg_TT_arrow_pilon = "grey", hcontortus_chr5_Celeg_TT_arrow_pilon = "grey", hcontortus_chrX_Celeg_TT_arrow_pilon = "grey")

names <- c("tci2_wsi2.0_chr_1", "tci2_wsi2.0_chr_2", "tci2_wsi2.0_chr_3", "tci2_wsi2.0_chr_4", "tci2_wsi2.0_chr_5", "tci2_wsi2.0_chr_X", "hcontortus_chr1_Celeg_TT_arrow_pilon", "hcontortus_chr2_Celeg_TT_arrow_pilon", "hcontortus_chr3_Celeg_TT_arrow_pilon", "hcontortus_chr4_Celeg_TT_arrow_pilon", "hcontortus_chr5_Celeg_TT_arrow_pilon", "hcontortus_chrX_Celeg_TT_arrow_pilon")

colours <- c("#d0e11c", "#73d056", "#2db27d", "#21918c", "#2e6e8e", "#46327e", "grey", "grey", "grey", "grey", "grey", "grey")

names_colours <- data.frame(names, colours)


tc_links_col <- left_join(tc_links, names_colours, by=c("chr" = "names"))


pdf("tc2.4_haem_promer_circulize.pdf")
circos.clear()

circos.genomicInitialize(chr)
circos.track(ylim = c(0, 1), bg.col = grid.col, bg.border = NA, track.height = 0.075)
circos.genomicLink(haem_links, tc_links_col, border = NA, col = tc_links_col$colours)

circos.clear()
dev.off()
```




## Counting bases mapped in paf between Haemonchus and tc2.4
```bash 
cat tc2.4_v_haem.paf | awk '{print $1, $6, $11}' OFS="\t" | datamash --sort -g 1,2 sum 3
hcontortus_chr1_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_1	2962004
hcontortus_chr1_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_2	282667
hcontortus_chr1_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_3	290135
hcontortus_chr1_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_4	421295
hcontortus_chr1_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_5	317339
hcontortus_chr1_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_X	240613
hcontortus_chr2_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_1	418165
hcontortus_chr2_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_2	2169528
hcontortus_chr2_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_3	359760
hcontortus_chr2_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_4	404987
hcontortus_chr2_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_5	394526
hcontortus_chr2_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_X	235971
hcontortus_chr3_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_1	950880
hcontortus_chr3_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_2	270074
hcontortus_chr3_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_3	2217859
hcontortus_chr3_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_4	354561
hcontortus_chr3_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_5	348071
hcontortus_chr3_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_X	172252
hcontortus_chr4_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_1	332126
hcontortus_chr4_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_2	290435
hcontortus_chr4_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_3	376122
hcontortus_chr4_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_4	2604814
hcontortus_chr4_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_5	384548
hcontortus_chr4_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_X	213963
hcontortus_chr5_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_1	372399
hcontortus_chr5_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_2	303603
hcontortus_chr5_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_3	822999
hcontortus_chr5_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_4	399044
hcontortus_chr5_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_5	1800094
hcontortus_chr5_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_X	207964
hcontortus_chrX_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_1	687010
hcontortus_chrX_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_2	607586
hcontortus_chrX_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_3	624936
hcontortus_chrX_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_4	886592
hcontortus_chrX_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_5	665806
hcontortus_chrX_Celeg_TT_arrow_pilon	tci2_wsi2.0_chr_X	1844483

```

mkdir /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/GENOME_COMPARISON/ORTHOFINDER

cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/GENOME_COMPARISON/ORTHOFINDER

cp ../../ANNOTATION/BRAKER_V_AUGUSTUS/braker_augustus.apollo.2.4.clean.proteins.fa .
cat braker_augustus.apollo.2.4.clean.proteins.fa | cut -f1 -d " " > tc2.4.proteins.fa

bsub.py --queue long 5 get_unique "bash ./run_get_unique.sh"

# where "run_get_unique.sh" contains:
# cat tc2.4.proteins.fa | grep ">" | grep "\-00001" | sed 's/>//g' | cut -f1 -d " " | while read -r ID; do samtools faidx tc2.4.proteins.fa ${ID} >> tc2.4.unique.proteins.fa; done 

sed -i 's/\.//g' tc2.4.unique.proteins.fa

wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/teladorsagia_circumcincta/PRJNA72569/teladorsagia_circumcincta.PRJNA72569.WBPS18.protein.fa.gz    ## /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/GENOME_COMPARISON/ORTHOFINDER

gunzip teladorsagia_circumcincta.PRJNA72569.WBPS18.protein.fa.gz
mv teladorsagia_circumcincta.PRJNA72569.WBPS18.protein.fa proteins/

gunzip haemonchus_contortus.PRJEB506.WBPS18.protein.fa.gz 

# modified get unique for haem
bsub.py 10 get_unique_hc "bash run_get_unique.sh"

mv hc.unique.proteins.fa proteins/



module load orthofinder/2.3.3--0
bsub.py 50 --threads 20 --queue long orthofinder "orthofinder -f proteins/"


cd ~/lustre_link/teladorsagia_circumcincta/GENOME/GENOME_COMPARISON/ORTHOFINDER/proteins/OrthoFinder/Results_Nov10_1/Orthologues/Orthologues_tc2.4.unique.proteins

cat tc2.4.unique.proteins__v__hc.unique.proteins.tsv |\ 
    awk '{if(NF==3) print $2,$3}' OFS="\t" > hc_v_tc2.4_1-to-1-orthologs.list


/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/GENOME_COMPARISON/ORTHOFINDER/proteins/OrthoFinder/Results_Nov12/Orthologues/Orthologues_tc2.4.unique.proteins/hc_v_tc2.4_1-to-1-orthologs.list


wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3.gz
gunzip haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3

mv haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3 hc.gff3

ln -s ../../ANNOTATION/BRAKER_V_AUGUSTUS/braker_augustus.apollo.2.4.clean.gff3 tc2.4.gff3


ln -s /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/GENOME_COMPARISON/ORTHOFINDER/proteins/OrthoFinder/Results_Nov10_1/Orthologues/Orthologues_tc2.4.unique.proteins/hc_v_tc2.4_1-to-1-orthologs.list

# get gene coordinates
awk -F '[\t;]' '$3=="mRNA" {print $1, $4, $5, $7, $9}' OFS="\t" tc2.4.gff3 | sed 's/ID=//g' > tc2.4.gff.coords
awk -F '[\t;]' '$3=="mRNA" {print $1, $4, $5, $7, $9}' OFS="\t" hc.gff3 | sed 's/ID=transcript://g' > hc.gff.coords

# some strange hidden characters need to be removed
sed -i $'s/[^[:print:]\t]//g' hc_v_tc2.4_1-to-1-orthologs.list

# extract positional data for both Hc and Ce
cat hc_v_tc2.4_1-to-1-orthologs.list | while read -r tcirc haem ; do 
    COORDS1=$(grep "${tcirc}" tc2.4.gff.coords) ; 
    COORDS2=$(grep "${haem}" hc.gff.coords) ; 
    echo -e "${COORDS1}\t${COORDS2}"  ; 
    done > hc_v_tc2.4_1-to-1.hccoords


cat hc_v_tc2.4_1-to-1.hccoords | awk '{print $6, "100000000", $7,$8,$9, $1, "100000000", $2, $3, $3-$2 , "60"}' OFS="\t"  > hc_v_tc2.4_1-to-1.pseudo-paf

cat hc_v_tc2.4_1-to-1.pseudo-paf | /nfs/users/nfs_s/sd21/lustre_link/software/GENOME_ASSEMBLY/minimap/utils/bin/layout | cut -f1-19 > tc2.4_v_hc.promer.layout

cat tc2.4_v_hc.promer.layout | grep -v "scaf" | grep -v "_lg_" > tc2.4_v_hc.promer.chr.layout

```R
library(tidyverse)
library(circlize)
library(viridis)

# tc2.4 vs haem
chr <- read.table("tc_2.4_haem_chromosome_lengths.txt", header=F)
colnames(chr) <- c("chr", "start", "end")
data<-read.table("tc2.4_v_hc.promer.chr.layout")
#data<-data %>% filter(V18>2000)

haem_links <- data %>% select(V8, V10, V11, V18)
colnames(haem_links) <- c("chr", "start", "end", "value")
tc_links <- data %>% select(V13, V15, V16, V18)
colnames(tc_links) <- c("chr", "start", "end", "value")


grid.col = c(tci2_wsi2.0_chr_1 = "#d0e11c", tci2_wsi2.0_chr_2 = "#73d056", tci2_wsi2.0_chr_3 = "#2db27d", tci2_wsi2.0_chr_4 = "#21918c", tci2_wsi2.0_chr_5 = "#2e6e8e", tci2_wsi2.0_chr_X = "#46327e", hcontortus_chr1_Celeg_TT_arrow_pilon = "grey", hcontortus_chr2_Celeg_TT_arrow_pilon = "grey", hcontortus_chr3_Celeg_TT_arrow_pilon = "grey", hcontortus_chr4_Celeg_TT_arrow_pilon = "grey", hcontortus_chr5_Celeg_TT_arrow_pilon = "grey", hcontortus_chrX_Celeg_TT_arrow_pilon = "grey")

names <- c("tci2_wsi2.0_chr_1", "tci2_wsi2.0_chr_2", "tci2_wsi2.0_chr_3", "tci2_wsi2.0_chr_4", "tci2_wsi2.0_chr_5", "tci2_wsi2.0_chr_X", "hcontortus_chr1_Celeg_TT_arrow_pilon", "hcontortus_chr2_Celeg_TT_arrow_pilon", "hcontortus_chr3_Celeg_TT_arrow_pilon", "hcontortus_chr4_Celeg_TT_arrow_pilon", "hcontortus_chr5_Celeg_TT_arrow_pilon", "hcontortus_chrX_Celeg_TT_arrow_pilon")

colours <- c("#d0e11c", "#73d056", "#2db27d", "#21918c", "#2e6e8e", "#46327e", "grey", "grey", "grey", "grey", "grey", "grey")

names_colours <- data.frame(names, colours)


tc_links_col <- left_join(tc_links, names_colours, by=c("chr" = "names"))


pdf("tc2.4_haem_ortholog_circulize.pdf")
circos.clear()

circos.genomicInitialize(chr)
circos.track(ylim = c(0, 1), bg.col = grid.col, bg.border = NA, track.height = 0.075)
circos.genomicLink(haem_links, tc_links_col, border = NA, col = tc_links_col$colours)

circos.clear()
dev.off()
```

 sort -k1,1 -k2,2n hc_v_tc2.4_1-to-1.hccoords | awk '{if($2+0==$2 && $7+0==$7) print $0}' | grep -v "scaf" | grep -v "_lg_" > hc_v_tc2.4_1-to-1.hccoords.chrsorted

 ```R
library(ggplot2)
library(patchwork)
library(dplyr) 
 
data <- read.table("hc_v_tc2.4_1-to-1.hccoords.chrsorted", header=F)
plot <- ggplot(data, aes(V2/1e6, V7/1e6, group=V1))+
	geom_point(aes(col=V6), size=0.5)+
	facet_grid(.~V1)+
	labs(x="T. circumcincta chromosome position (bp)", y="H. contortus chromosome position (bp)")+
	theme_bw()
 
 plot
 
 ```