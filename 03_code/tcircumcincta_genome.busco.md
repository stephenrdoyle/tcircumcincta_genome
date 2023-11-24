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

bsub.py --queue long --threads 20 60 busco_tc_washu_genome_metazoa_odb10 \
    "busco --in teladorsagia_circumcincta.PRJNA72569.WBPS17.genomic.fa --out washu_genome_metazoa_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/metazoa_odb10 --cpu 20 -f -r"
	--------------------------------------------------
	|Results from dataset metazoa_odb10               |
	--------------------------------------------------
	|C:49.7%[S:40.4%,D:9.3%],F:15.9%,M:34.4%,n:954    |
	|474	Complete BUSCOs (C)                       |
	|385	Complete and single-copy BUSCOs (S)       |
	|89	Complete and duplicated BUSCOs (D)        |
	|152	Fragmented BUSCOs (F)                     |
	|328	Missing BUSCOs (M)                        |
	|954	Total BUSCO groups searched               |
	--------------------------------------------------

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


bsub.py --queue long --threads 20 60 busco_tc_dnazoo_genome_metazoa_odb10 \
    "busco --in Tcircumcincta.DNAzoo.fa --out dnazoo_genome_metazoa_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/metazoa_odb10 --cpu 20 -f -r"

	--------------------------------------------------
	|Results from dataset metazoa_odb10               |
	--------------------------------------------------
	|C:47.8%[S:42.8%,D:5.0%],F:15.9%,M:36.3%,n:954    |
	|456	Complete BUSCOs (C)                       |
	|408	Complete and single-copy BUSCOs (S)       |
	|48	Complete and duplicated BUSCOs (D)        |
	|152	Fragmented BUSCOs (F)                     |
	|346	Missing BUSCOs (M)                        |
	|954	Total BUSCO groups searched               |
	--------------------------------------------------
```



## RAGTAG

```bash
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/PLAY/ragtag_output

bsub.py --queue long --threads 20 60 busco_tc_ragtag_genome_nematoda_odb10 \
    "busco --in ragtag.scaffold.fasta --out ragtag_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_tc_ragtag_genome_eukaryota_odb10 \
    "busco --in ragtag.scaffold.fasta --out ragtag_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_tc_ragtag_genome_metazoa_odb10 \
    "busco --in ragtag.scaffold.fasta --out ragtag_genome_metazoa_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/metazoa_odb10 --cpu 20 -f -r"
```


## Haemonchus contortus
```bash
ln -s ../../../haemonchus_contortus/GENOME/REF/HAEM_V4_final.chr.fa

bsub.py --queue long --threads 20 60 busco_hcontortus_genome_nematoda_odb10 \
    "busco --in HAEM_V4_final.chr.fa --out hcontortus_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_hcontortus_genome_eukaryota_odb10 \
    "busco --in HAEM_V4_final.chr.fa --out hcontortus_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"

bsub.py --queue long --threads 20 60 busco_hcontortus_genome_metazoa_odb10 \
    "busco --in HAEM_V4_final.chr.fa --out hcontortus_genome_metazoa_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/metazoa_odb10 --cpu 20 -f -r"


	--------------------------------------------------
	|Results from dataset metazoa_odb10               |
	--------------------------------------------------
	|C:61.4%[S:60.5%,D:0.9%],F:8.4%,M:30.2%,n:954     |
	|586	Complete BUSCOs (C)                       |
	|577	Complete and single-copy BUSCOs (S)       |
	|9	Complete and duplicated BUSCOs (D)        |
	|80	Fragmented BUSCOs (F)                     |
	|288	Missing BUSCOs (M)                        |
	|954	Total BUSCO groups searched               |
	--------------------------------------------------	
```





## Caenorhabditis elegans

```bash

ln -s ../../../REFERENCE_SEQUENCES/caenorhabditis_elegans/caenorhabditis_elegans.PRJNA13758.WBPS9.genomic_softmasked.fa

bsub.py --queue long --threads 20 60 busco_celegans_genome_nematoda_odb10 \
    "busco --in caenorhabditis_elegans.PRJNA13758.WBPS9.genomic_softmasked.fa --out celegans_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_celegans_genome_eukaryota_odb10 \
    "busco --in caenorhabditis_elegans.PRJNA13758.WBPS9.genomic_softmasked.fa --out celegans_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"

bsub.py --queue long --threads 20 60 busco_celegans_genome_metazoa_odb10 \
    "busco --in caenorhabditis_elegans.PRJNA13758.WBPS9.genomic_softmasked.fa --out celegans_genome_metazoa_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/metazoa_odb10 --cpu 20 -f -r"

	--------------------------------------------------
	|Results from dataset metazoa_odb10               |
	--------------------------------------------------
	|C:71.0%[S:70.0%,D:1.0%],F:3.2%,M:25.8%,n:954     |
	|678	Complete BUSCOs (C)                       |
	|668	Complete and single-copy BUSCOs (S)       |
	|10	Complete and duplicated BUSCOs (D)        |
	|31	Fragmented BUSCOs (F)                     |
	|245	Missing BUSCOs (M)                        |
	|954	Total BUSCO groups searched               |
	--------------------------------------------------
```



## Tcirc Canu_1.9

```bash 

ln -s /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/PACBIO_CANU1.9v2/tc_canu1.9.contigs.fasta

bsub.py --queue long --threads 20 30 busco_tc_canu1.9_genome_nematoda_odb10 "busco --in tc_canu1.9.contigs.fasta --out tc_canu1.9_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 30 busco_tc_canu1.9_genome_eukaryota_odb10 "busco --in tc_canu1.9.contigs.fasta --out tc_canu1.9_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"

bsub.py --queue long --threads 20 30 busco_tc_canu1.9_genome_metazoa_odb10 "busco --in tc_canu1.9.contigs.fasta --out tc_canu1.9_genome_metazoa_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/metazoa_odb10 --cpu 20 -f -r"

```




bsub.py --queue long --threads 20 60 busco_tc_genome_nematoda_odb10 \
    "busco --in Tc22_2_ragtag_mc.curated_primary.no_mt.unscrubbed.fa --out genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"