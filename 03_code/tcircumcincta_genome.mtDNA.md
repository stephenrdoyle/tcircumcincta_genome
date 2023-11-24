# Fix mitochondrial genome

### Stephen Doyle


- need to curate the mtDNA genome

```bash
# working directory
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/FIX_MTDNA

```

## Finding the mtDNA genome in the assembly
```bash
# get the final genome
cp ../RECOVER_MISSING_BUSCOs/CURRENT_ASSEMBLY.fa

# get a reference tcircumcincta mtDNA genome from NCBI
#- downloaded NC_013827.1


# run nucmer to identify the mitochondrial genome
nucmer CURRENT_ASSEMBLY.fa tc_mtDNA_NC_013827.1.fasta

# make coords
show-coords out.delta > mtDNA.coords
```

    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
=====================================================================================
       1    10483  |     3537    14066  |    10483    10530  |    97.05  | SUPER_33	NC_013827.1
   10485    24522  |        2    14066  |    14038    14065  |    97.76  | SUPER_33	NC_013827.1
49861218 49861286  |     7765     7697  |       69       69  |    98.55  | SUPER_5	NC_013827.1

- mtDNA is on SUPER_33,
- SUPER_33 is bigger than expected, likely caused by more than one mtDNA genome in a tandem arrau
- a single copy lies between 10485 and 24522


```bash
# samtools to extract the mtDNA 
samtools faidx CURRENT_ASSEMBLY.fa SUPER_33:10484-24522 > tc_assembly_mtDNA_genome.fa

assembly-stats tc_assembly_mtDNA_genome.fa
#stats for tc_assembly_mtDNA_genome.fa
#sum = 14039, n = 1, ave = 14039.00, largest = 14039
#N50 = 14039, n = 1
#N60 = 14039, n = 1
#N70 = 14039, n = 1
#N80 = 14039, n = 1
#N90 = 14039, n = 1
#N100 = 14039, n = 1
#N_count = 0
#Gaps = 0
```

- genome is slightly smaller than the reference of 14066
    - might contain indels, and so size might change after correction

- changed the sequence name to ">tcircumcincta_chr_mtDNA"
- checked the starting coords - starts at CO1 as per the original draft. Had to fix the start codon from TTT to ATT (Met)



## Replace mtDNA in genome

```bash
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/PRE_POLISH_COMPLETE 

# get mtDNA
cp ../FIX_MTDNA/tc_assembly_mtDNA_genome.fa .

# get genome 
cp /lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/RECOVER_MISSING_BUSCOs/busco_recovered.fa CURRENT_ASSEMBLY.fa

samtools faidx CURRENT_ASSEMBLY.fa

echo "SUPER_33" > blocklist.txt

remove_ids=($(awk '{print $1}' CURRENT_ASSEMBLY.fa.fai | grep -x -v -f blocklist.txt))

samtools faidx -o output.fasta CURRENT_ASSEMBLY.fa "${remove_ids[@]}"

cat output.fasta tc_assembly_mtDNA_genome.fa > CURRENT_ASSEMBLY_mtDNA.fa
```



stats for CURRENT_ASSEMBLY_mtDNA.fa
sum = 644904321, n = 3382, ave = 190687.26, largest = 94624076
N50 = 83760063, n = 4
N60 = 77880069, n = 5
N70 = 62249764, n = 6
N80 = 1856228, n = 11
N90 = 64366, n = 362
N100 = 1199, n = 3382
N_count = 3982291
Gaps = 12486


