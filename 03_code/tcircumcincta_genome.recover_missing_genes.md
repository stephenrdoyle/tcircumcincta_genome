# Recover missing genes


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

ln -s ~/lustre_link/teladorsagia_circumcincta/PLAY/ragtag_output/ragtag_genome_nematoda_odb10/run_nematoda_odb10/missing_busco_list.tsv busco_ragtag_missing.txt

ln -s ~/lustre_link/teladorsagia_circumcincta/GENOME/BUSCO/tc_canu1.9_genome_nematoda_odb10/run_nematoda_odb10/missing_busco_list.tsv busco_tc_canu1.9_missing.txt

ln -s /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/BUSCO/tc_canu1.9_genome_nematoda_odb10/run_nematoda_odb10/full_table.tsv busco_tc_canu1.9_full_table.tsv

cp ../../../BUSCO/tc_canu1.9.contigs.fasta .

# get list of busco IDs present in tc_canu1.9 but missing in ragtag
diff -y  <(sort busco_ragtag_missing.txt) <(sort busco_tc_canu1.9_missing.txt) | grep "<" | cut -f1 | sort | uniq > tc_canu1.9.uniq.busco_ids.list


# extract contig data per gene, and when there is more than one hit, take the longest contig ID
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
- 207 genes in 196 contigs spanning 21.7 Mb



```bash 
ln -s /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/PLAY/ragtag_output/ragtag.scaffold.fasta

cat ragtag.scaffold.fasta tc_canu1.9.uniq.busco_containing_contigs.all.fasta > tmp.fasta

fastaq to_fasta -l0 tmp.fasta ragtag.busco_recovered.fasta
```

sum = 661818020, n = 3490, ave = 189632.67, largest = 111853237
N50 = 86408422, n = 4
N60 = 79992886, n = 5
N70 = 79992886, n = 5
N80 = 65169149, n = 6
N90 = 62166, n = 365
N100 = 1199, n = 3490
N_count = 4119207
Gaps = 12722


## busco 
```bash

bsub.py --queue long --threads 20 40 busco_ragtag.busco_recovered "busco --in ragtag.busco_recovered.fasta --out ragtag.busco_recovered_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"

```


bsub.py 10 nucmer2 "nucmer ragtag.scaffold.fasta tc_canu1.9.uniq.busco_containing_contigs.all.fasta"

show-coords -ol out.delta | cat | grep "CONTAIN"

       1    13935  |   205844   219788  |    13935    13945  |    99.62  |    13935   510944  | scaffold_315	tig00000571	[CONTAINED]
72714653 72741640  |      683    27601  |    26988    26919  |    99.13  | 111853237    27601  | hcontortus_chr4_Celeg_TT_arrow_pilon_RagTag	tig00019024	[CONTAINS]
85771187 85795504  |    24383        1  |    24318    24383  |    98.87  | 86408422    24387  | hcontortus_chr3_Celeg_TT_arrow_pilon_RagTag	tig00054006	[CONTAINS]
       1    26861  |   974671  1001554  |    26861    26884  |    99.78  |    26861  1010514  | scaffold_558	tig00181744	[CONTAINED]
69209291 69231123  |        1    21872  |    21833    21872  |    98.37  | 86408422    21872  | hcontortus_chr3_Celeg_TT_arrow_pilon_RagTag	tig00187062	[CONTAINS]
       1    50864  |    52652     1705  |    50864    50948  |    99.37  |    50864   223259  | scaffold_2251	tig00760077	[CONTAINED]
     819    18261  |    17452        1  |    17443    17452  |    99.56  |    18746   174081  | scaffold_3317	tig00760163	[CONTAINED]
       1    13088  |   134492   121348  |    13088    13145  |    98.35  |    13088   231573  | scaffold_3733	tig00760233	[CONTAINED]
30064583 30096221  |    31640        1  |    31639    31640  |    98.77  | 86408422    31640  | hcontortus_chr3_Celeg_TT_arrow_pilon_RagTag	tig00763292	[CONTAINS]
111790110 111801042  |    10978        1  |    10933    10978  |    97.10  | 111853237    11119  | hcontortus_chr4_Celeg_TT_arrow_pilon_RagTag	tig00763963	[CONTAINS]
51079501 51104399  |    24896        1  |    24899    24896  |    99.35  | 86408422    24896  | hcontortus_chr3_Celeg_TT_arrow_pilon_RagTag	tig00764941	[CONTAINS]


to do
- remove based on contains
    - scaffold_315
    - tig00019024
    - tig00054006 - ?
    - scaffold_558
    - tig00187062
    - scaffold_2251
    - scaffold_3317
    - scaffold_3733
    - tig00763292
    - tig00763963
    - tig00764941
- to remove becasue sig overlap
    - tig00000547
    - tig00001010
    - tig00002238
    - tig00002291
    - tig00003764
    - tig00003930
    - tig00004772
    - tig00005297
    - tig00002106
    - tig00184120
    - tig00181817
