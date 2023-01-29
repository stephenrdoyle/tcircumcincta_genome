# Teladorsagia circumcincta genome analysis: BUSCOs

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk





```bash
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/BUSCO

conda activate busco_5.4.3

```

## WASHU
```bash
ln -s ../REFERENCES/WASHU/teladorsagia_circumcincta.PRJNA72569.WBPS17.genomic.fa


bsub.py --queue long --threads 20 60 busco_tc_washu_genome_nematoda_odb10 \
    "busco --in teladorsagia_circumcincta.PRJNA72569.WBPS17.genomic.fa --out washu_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_tc_washu_genome_eukaryota_odb10 \
    "busco --in teladorsagia_circumcincta.PRJNA72569.WBPS17.genomic.fa --out washu_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"
```


## DNAZOO
```bash
ln -s ../REFERENCES/DNAZOO/Tcircumcincta.DNAzoo.fa

bsub.py --queue long --threads 20 60 busco_tc_dnazoo_genome_nematoda_odb10 \
    "busco --in Tcircumcincta.DNAzoo.fa --out dnazoo_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_tc_dnazoo_genome_eukaryota_odb10 \
    "busco --in Tcircumcincta.DNAzoo.fa --out dnazoo_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"

	--------------------------------------------------
	|Results from dataset eukaryota_odb10             |
	--------------------------------------------------
	|C:62.3%[S:58.4%,D:3.9%],F:25.1%,M:12.6%,n:255    |
	|159	Complete BUSCOs (C)                       |
	|149	Complete and single-copy BUSCOs (S)       |
	|10	Complete and duplicated BUSCOs (D)        |
	|64	Fragmented BUSCOs (F)                     |
	|32	Missing BUSCOs (M)                        |
	|255	Total BUSCO groups searched               |
	--------------------------------------------------
```



## RAGTAG

```bash
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/PLAY/ragtag_output

bsub.py --queue long --threads 20 60 busco_tc_ragtag_genome_nematoda_odb10 \
    "busco --in ragtag.scaffold.fasta --out ragtag_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_tc_ragtag_genome_eukaryota_odb10 \
    "busco --in ragtag.scaffold.fasta --out ragtag_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"

```


## Haemonchus contortus
```bash
ln -s ../../../haemonchus_contortus/GENOME/REF/HAEM_V4_final.chr.fa

bsub.py --queue long --threads 20 60 busco_hcontortus_genome_nematoda_odb10 \
    "busco --in HAEM_V4_final.chr.fa --out hcontortus_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_hcontortus_genome_eukaryota_odb10 \
    "busco --in HAEM_V4_final.chr.fa --out hcontortus_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"


```





## Caenorhabditis elegans

```bash

ln -s ../../../REFERENCE_SEQUENCES/caenorhabditis_elegans/caenorhabditis_elegans.PRJNA13758.WBPS9.genomic_softmasked.fa

bsub.py --queue long --threads 20 60 busco_celegans_genome_nematoda_odb10 \
    "busco --in caenorhabditis_elegans.PRJNA13758.WBPS9.genomic_softmasked.fa --out celegans_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_celegans_genome_eukaryota_odb10 \
    "busco --in caenorhabditis_elegans.PRJNA13758.WBPS9.genomic_softmasked.fa --out celegans_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"


```



## Tcirc Canu_1.9

```bash 

ln -s /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/PACBIO_CANU1.9v2/tc_canu1.9.contigs.fasta

bsub.py --queue long --threads 20 30 busco_tc_canu1.9_genome_nematoda_odb10 "busco --in tc_canu1.9.contigs.fasta --out tc_canu1.9_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 30 busco_tc_canu1.9_genome_eukaryota_odb10 "busco --in tc_canu1.9.contigs.fasta --out tc_canu1.9_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"



```