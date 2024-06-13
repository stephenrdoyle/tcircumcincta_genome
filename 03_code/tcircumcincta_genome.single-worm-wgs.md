# Teladorsagia circumcincta genome: single-worm whole genome sequencing

### Stephen Doyle



## Notes on Tcirc single worm mapping
- 48 samples sequenced
- Individual worms from the following strains:
- MTci1
    - single adults, from a post-BZ population (MTci1_single-adult_post-BZ) x10
- MTci2
    - single L4s, that are drug naive and drug susceptible (MTci2_single-L4) x10
- MTci5
    - single adults, from a post-BZ population (MTci5_single-adult_post-BZ) x9 
- MTci5
    single adults, from a post-BZ population (MTci5_single-adultM_post-IVM) x9
- MTci7_single-adult_MOX-R
    - single adults, from a MOX-R population (MTci7_single-adult_MOX-R) x10


## Get the data from iRODs
```bash
cd ~/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW


# need to log in to iRODs to refresh access so bato will work

iinit

conda deactivate

module load ISG/singularity nextflow

nextflow run /nfs/users/nfs_s/sd21/lustre_link/software/nextflow-tools/baton-extractor/main.nf --study 7548 --runid 48232

# fastq.gz files end up here: /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq
```


## Mapping to the reference genome

```bash
module load mapping-helminth/v1.0.8

bsub.py --queue long 10 mapping "mapping-helminth --input tc_singleworm_mapping_manifest.txt --reference teladorsagia_circumcincta_tci2_wsi2.4.fa"
```

``` 
- where "tc_singleworm_mapping_manifest.txt" is: 
```bash
ID,R1,R2
MTci1_single-adult_post-BZ_S01,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#1_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#1_2.fastq.gz
MTci1_single-adult_post-BZ_S02,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#2_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#2_2.fastq.gz
MTci1_single-adult_post-BZ_S03,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#3_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#3_2.fastq.gz
MTci1_single-adult_post-BZ_S04,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#4_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#4_2.fastq.gz
MTci1_single-adult_post-BZ_S05,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#5_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#5_2.fastq.gz
MTci1_single-adult_post-BZ_S06,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#6_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#6_2.fastq.gz
MTci1_single-adult_post-BZ_S07,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#7_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#7_2.fastq.gz
MTci1_single-adult_post-BZ_S08,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#8_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#8_2.fastq.gz
MTci1_single-adult_post-BZ_S09,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#9_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#9_2.fastq.gz
MTci1_single-adult_post-BZ_S10,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#10_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#10_2.fastq.gz
MTci2_single-L4_S01,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#11_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#11_2.fastq.gz
MTci2_single-L4_S02,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#12_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#12_2.fastq.gz
MTci2_single-L4_S03,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#13_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#13_2.fastq.gz
MTci2_single-L4_S04,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#14_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#14_2.fastq.gz
MTci2_single-L4_S05,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#15_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#15_2.fastq.gz
MTci2_single-L4_S06,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#16_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#16_2.fastq.gz
MTci2_single-L4_S07,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#17_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#17_2.fastq.gz
MTci2_single-L4_S08,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#18_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#18_2.fastq.gz
MTci2_single-L4_S09,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#19_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#19_2.fastq.gz
MTci2_single-L4_S10,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#20_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#20_2.fastq.gz
MTci5_single-adult_post-BZ_S02,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#21_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#21_2.fastq.gz
MTci5_single-adult_post-BZ_S03,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#22_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#22_2.fastq.gz
MTci5_single-adult_post-BZ_S04,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#23_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#23_2.fastq.gz
MTci5_single-adult_post-BZ_S05,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#24_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#24_2.fastq.gz
MTci5_single-adult_post-BZ_S06,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#25_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#25_2.fastq.gz
MTci5_single-adult_post-BZ_S07,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#26_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#26_2.fastq.gz
MTci5_single-adult_post-BZ_S08,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#27_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#27_2.fastq.gz
MTci5_single-adult_post-BZ_S09,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#28_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#28_2.fastq.gz
MTci5_single-adult_post-BZ_S10,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#29_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#29_2.fastq.gz
MTci5_single-adultM_post-IVM_S01,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#30_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#30_2.fastq.gz
MTci5_single-adultM_post-IVM_S02,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#31_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#31_2.fastq.gz
MTci5_single-adultM_post-IVM_S04,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#32_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#32_2.fastq.gz
MTci5_single-adultM_post-IVM_S05,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#33_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#33_2.fastq.gz
MTci5_single-adultM_post-IVM_S06,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#34_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#34_2.fastq.gz
MTci5_single-adultM_post-IVM_S07,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#35_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#35_2.fastq.gz
MTci5_single-adultM_post-IVM_S08,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#36_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#36_2.fastq.gz
MTci5_single-adultM_post-IVM_S09,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#37_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#37_2.fastq.gz
MTci5_single-adultM_post-IVM_S10,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#38_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#38_2.fastq.gz
MTci7_single-adult_MOX-R_S01,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#39_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#39_2.fastq.gz
MTci7_single-adult_MOX-R_S02,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#40_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#40_2.fastq.gz
MTci7_single-adult_MOX-R_S03,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#41_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#41_2.fastq.gz
MTci7_single-adult_MOX-R_S04,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#42_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#42_2.fastq.gz
MTci7_single-adult_MOX-R_S05,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#43_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#43_2.fastq.gz
MTci7_single-adult_MOX-R_S06,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#44_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#44_2.fastq.gz
MTci7_single-adult_MOX-R_S07,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#45_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#45_2.fastq.gz
MTci7_single-adult_MOX-R_S08,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#46_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#46_2.fastq.gz
MTci7_single-adult_MOX-R_S09,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#47_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#47_2.fastq.gz
MTci7_single-adult_MOX-R_S10,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#48_1.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/RAW/results/raw_fastq/raw_48232_4#48_2.fastq.gz
```


## Sequnecing QC
```bash
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM/results

multiqc .

```
[Multiqc report](../04_analysis/tc_single-worm_multiqc_report.html) 

- good sequencing results for most samples
- very poor for the Tci2 L4 
    - these were significantly smaller than the other adult samples, and so perhaps makes sense.



## Variant calling

#run_genotyping.sh 
```bash
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/DRUG_POPGEN/SINGLE_WORM

ls -1 $PWD/results/MTci*/*bam > bam.list

# edit the header to add the prefix, reference, and path to bam list
./run_genotyping.sh