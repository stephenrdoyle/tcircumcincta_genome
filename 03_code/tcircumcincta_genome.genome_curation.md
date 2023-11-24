# Cleaning up the genome



## Recovering missing genes


- found a large discrepency between the number of BUSCOs recovered in the original tcirc pacbio assembly and the final assembly being curated
    - original: tc_canu1.9.contigs.fasta
    - curation: 
- want to try and recover these missing BUSCOs, and in turn, should be able to recover other nearby genes
- approach is to 
    - find out which BUSCOs are present in the original assembly, but missing in the curated assembly
    - extract the contigs from the original assembly
    - merge them with the new assembly as extra sequence (bag of bits)
    - check new BUSCO stats
- if this works, would be worth running HiC mapping again, to see if any of these addditonal contigs can be correctly integrated 



```bash 

cd ~/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/RECOVER_MISSING_BUSCOs

ln -s /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/MANUAL_CURATION/genome_nematoda_odb10/run_nematoda_odb10/missing_busco_list.tsv busco_ragtag_missing.txt

ln -s ~/lustre_link/teladorsagia_circumcincta/GENOME/BUSCO/tc_canu1.9_genome_nematoda_odb10/run_nematoda_odb10/missing_busco_list.tsv busco_tc_canu1.9_missing.txt

ln -s /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/BUSCO/tc_canu1.9_genome_nematoda_odb10/run_nematoda_odb10/full_table.tsv busco_tc_canu1.9_full_table.tsv

cp ../../../BUSCO/tc_canu1.9.contigs.fasta .

# get list of busco IDs present in tc_canu1.9 but missing in ragtag
diff -y  <(sort busco_ragtag_missing.txt) <(sort busco_tc_canu1.9_missing.txt) | grep "<" | cut -f1 | sort | uniq > tc_canu1.9.uniq.busco_ids.list


# extract contig data per gene, and when there is more than one hit, take the longest contig ID
>tmp2
while read -r BUSCO_ID; do 
    grep "$BUSCO_ID" busco_tc_canu1.9_full_table.tsv | cut -f3 > tmp 
    grep -f tmp tc_canu1.9.contigs.fasta.fai | sort -k2,2nr | head -n1 | cut -f1 >> tmp2
    done < tc_canu1.9.uniq.busco_ids.list 

# using the longest contig ID, extract the contig fasta    
> tc_canu1.9.uniq.busco_containing_contigs.all.fasta
cat tmp2 | sort | uniq | while read CONTIG_ID; do 
    samtools faidx tc_canu1.9.contigs.fasta ${CONTIG_ID} >> tc_canu1.9.uniq.busco_containing_contigs.all.fasta; 
    done
    

```
- 209 genes in 197 contigs spanning 21.8 Mb



```bash 
ln -s /nfs/users/nfs_j/jm62/lustre118_link/Tc22_ragtag/Tc22_ragtag_curated_out/Tc22_2_ragtag_mc.curated_primary.no_mt.unscrubbed.fa Tc_curated.fa

cat Tc_curated.fa tc_canu1.9.uniq.busco_containing_contigs.all.fasta > tmp.fasta

fastaq to_fasta -l0 tmp.fasta Tc_curated.busco_recovered.fasta

assembly-stats Tc_curated.busco_recovered.fasta
# stats for Tc_curated.busco_recovered.fasta
# sum = 661028114, n = 3566, ave = 185369.63, largest = 94624076
# N50 = 83760063, n = 4
# N60 = 77880069, n = 5
# N70 = 62249764, n = 6
# N80 = 747384, n = 23
# N90 = 59245, n = 437
# N100 = 1199, n = 3566
# N_count = 4063760
# Gaps = 12731

```


###§ busco 
```bash

bsub.py --queue long --threads 20 40 busco_Tc_curated.busco_recovered "busco --in Tc_curated.busco_recovered.fasta --out Tc_curated.busco_recovered_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"

```
	--------------------------------------------------
	|Results from dataset nematoda_odb10              |
	--------------------------------------------------
	|C:85.6%[S:73.9%,D:11.7%],F:7.8%,M:6.6%,n:3131    |
	|2680	Complete BUSCOs (C)                       |
	|2313	Complete and single-copy BUSCOs (S)       |
	|367	Complete and duplicated BUSCOs (D)        |
	|244	Fragmented BUSCOs (F)                     |
	|207	Missing BUSCOs (M)                        |
	|3131	Total BUSCO groups searched               |
	--------------------------------------------------

## Removing contigs and scaffolds based on significant overlap
- check whether the 
```bash

bsub.py 10 nucmer2 "nucmer Tc_curated.fa tc_canu1.9.uniq.busco_containing_contigs.all.fasta"

```

- remove - based on contained, significant overlap with chromosomes
scaffold_315
tig00019024
tig00054006
scaffold_558
scaffold_1960
tig00187062
scaffold_2251
scaffold_3317
scaffold_3733
tig00763292
tig00763963
tig00764941
tig00010307

## contigs and scaffolds based on presence of duplicated BUSCO data
- genes present in a chromosome and in extra sequences were compared
- busco data: /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/RECOVER_MISSING_BUSCOs/Tc_curated.busco_recovered.reduced.genome_nematoda_odb10/run_nematoda_odb10/full_table.tsv

- remove - based on duplicated buscos 
scaffold_454
scaffold_587
scaffold_3143
scaffold_677
tig00181848
scaffold_1353
scaffold_471
scaffold_1001
scaffold_3140
scaffold_1879
scaffold_2977
tig00763067
scaffold_237
scaffold_561
scaffold_2538
scaffold_598
scaffold_3337
scaffold_3732
tig00181713
scaffold_1603
scaffold_1197
scaffold_2412
scaffold_1914
scaffold_3565
scaffold_2102
scaffold_2255
scaffold_2339
scaffold_3719
scaffold_79
tig00013719
SUPER_4_unloc_7
scaffold_2585
scaffold_714
scaffold_202
scaffold_1332
tig00186286
scaffold_2079
scaffold_2754
scaffold_1671
scaffold_83
scaffold_3261
scaffold_2787
scaffold_928
scaffold_2054
scaffold_1675
SUPER_2_unloc_1
scaffold_598
scaffold_1172
SUPER_2_unloc_1
scaffold_3437
scaffold_388
scaffold_2667
SUPER_4_unloc_3
scaffold_2761
scaffold_2867
SUPER_4_unloc_5
scaffold_1982
scaffold_459
scaffold_1884
tig00760674
scaffold_1191
tig00182378
scaffold_1376
tig00028813
scaffold_1156
scaffold_3509
scaffold_909
scaffold_1700
scaffold_1176
scaffold_2410
scaffold_1583
scaffold_1260
scaffold_2851
scaffold_1429
scaffold_1209
scaffold_3274
scaffold_1463
scaffold_281
scaffold_3384
scaffold_3553
scaffold_147
scaffold_633
scaffold_2954
tig00001594
scaffold_1762
scaffold_2747
scaffold_2446
scaffold_3341
scaffold_554
SUPER_9
scaffold_900
scaffold_2873
scaffold_867
scaffold_1384
scaffold_2015
scaffold_2710
scaffold_2553
scaffold_726
scaffold_415
scaffold_2995
scaffold_3662
scaffold_2855
scaffold_3571
scaffold_3549
scaffold_1702
scaffold_77
scaffold_1624
scaffold_2041
scaffold_1885
scaffold_1760
scaffold_1314
tig00005828
SUPER_3_unloc_4
scaffold_941
scaffold_2372
scaffold_3750
scaffold_1172
scaffold_1376
scaffold_1156
tig00001010
tig00185468
tig00003764
tig00025056
tig00002106
tig00024900
tig00181765
tig00020405
tig00181744
tig00004772
tig00024403



- contigs highly duplicated BUSCO genes, but have some complete genes - can recover these genes and surrounding sequences, but ignore the duplicated genes, using the following coordinates 
tig00001010 15000   275000
tig00185468 40000   50000
tig00003764 1   11000
tig00025056 1   23500
tig00002106 80000   285066
tig00024900 0   20000
tig00181765 0   180000
tig00020405 20000   54150
tig00181744 0   100000
tig00004772 0   50000
tig00024403 50000   68799

- put these coordinates into a bed file "recovered.bed" for extracting sequences to keep

### extract nucleotide sequence on the recovered sections
```bash
>recovered.fa
while read -r NAME START END; do 
	samtools faidx Tc_curated.busco_recovered.fasta ${NAME}:${START}-${END} >> recovered.fa; 
	done < recovered.bed

```


## remove redundant scaffolds contigs based on nucmer overlap between chromosomes and extra sequences
- ran numcmer 
- workign dir: /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/FIX_MTDNA
- compared coords and percentage overlap 
	- clear why these didnt get picked up by other tools 
	- they seem divergent, & misassembled in general. Mostly rubbish really. 
	- some were clearly redundant. 
	- add these to the blocklist

```bash
# extract chromosomes
while read NAME; do 
    samtools faidx CURRENT_ASSEMBLY_mtDNA.fa ${NAME} >> chr.fa; 
    done < chr.list

# extract everything except the chromosomes
remove_ids=($(awk '{print $1}' CURRENT_ASSEMBLY_mtDNA.fa.fai | grep -x -v -f chr.list))

samtools faidx -o non-chromosomes.fa CURRENT_ASSEMBLY_mtDNA.fa "${remove_ids[@]}"


bsub.py --queue long 10 nucmer "nucmer chr.fa non-chromosomes.fa"
```


scaffold_1059
scaffold_1065
scaffold_1067
scaffold_1466
scaffold_150
scaffold_1789
scaffold_1810
scaffold_1963
scaffold_1994
scaffold_2084
scaffold_2085
scaffold_234
scaffold_24
scaffold_2485
scaffold_2554
scaffold_2680
scaffold_2695
scaffold_2699
scaffold_2825
scaffold_2839
scaffold_2884
scaffold_2903
scaffold_2978
scaffold_3035
scaffold_3154
scaffold_3155
scaffold_3201
scaffold_334
scaffold_3358
scaffold_3407
scaffold_3431
scaffold_3433
scaffold_3439
scaffold_3459
scaffold_3469
scaffold_3474
scaffold_3493
scaffold_3614
scaffold_3663
scaffold_37
scaffold_3701
scaffold_3763
scaffold_3768
scaffold_3791
scaffold_38
scaffold_470
scaffold_498
scaffold_546
scaffold_548
scaffold_1789
scaffold_3474
scaffold_3703
scaffold_56
scaffold_962
scaffold_766
scaffold_62
scaffold_604
scaffold_555
scaffold_506














## remove redundant sequences
```bash
# make a "blocklist.txt" using all of the listed scaffolds/contigs above



samtools faidx Tc_curated.busco_recovered.fasta

remove_ids=($(awk '{print $1}' Tc_curated.busco_recovered.fasta.fai | grep -x -v -f blocklist.txt))

samtools faidx -o output.fasta Tc_curated.busco_recovered.fasta "${remove_ids[@]}"

cat output.fasta recovered.fa > busco_recovered.fa

assembly-stats busco_recovered.fa

#stats for busco_recovered.fa
#sum = 644916545, n = 3382, ave = 190690.88, largest = 94624076
#N50 = 83760063, n = 4
#N60 = 77880069, n = 5
#N70 = 62249764, n = 6
#N80 = 1856228, n = 11
#N90 = 64366, n = 362
#N100 = 1199, n = 3382
#N_count = 3982291
#Gaps = 12486

```


## busco 
```bash

bsub.py --queue long --threads 20 40 busco_Tc_curated.busco_recovered.reduced "busco --in busco_recovered.fa --out Tc_curated.busco_recovered.reduced.genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"

```
	--------------------------------------------------
	|Results from dataset nematoda_odb10              |
	--------------------------------------------------
	|C:85.4%[S:79.1%,D:6.3%],F:7.6%,M:7.0%,n:3131     |
	|2675	Complete BUSCOs (C)                       |
	|2477	Complete and single-copy BUSCOs (S)       |
	|198	Complete and duplicated BUSCOs (D)        |
	|237	Fragmented BUSCOs (F)                     |
	|219	Missing BUSCOs (M)                        |
	|3131	Total BUSCO groups searched               |
	--------------------------------------------------




