# Teladorsagia circumcincta genome annotation

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk


- Approaches
    - Braker2 with Tcirc RNAseq data from Moredun
    - Braker2 with metazoan proteins
    - Braker2 with Haemonchus proteins



```bash

# working directory:
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION 

```

## Raw data

- data
    - RNAseq from Moredun

| sequence_centre_ID | Species                   | Lifestage    | dpi  |
|--------------------|---------------------------|--------------|------|
| DP_1               | Teladorsagia circumcincta | xL3          | na   |
| DP_2               | Teladorsagia circumcincta | xL3          | na   |
| DP_3               | Teladorsagia circumcincta | xL3          | na   |
| DP_4               | Teladorsagia circumcincta | L4           | 7    |
| DP_5               | Teladorsagia circumcincta | L4           | 7    |
| DP_6               | Teladorsagia circumcincta | L4           | 7    |
| DP_7               | Teladorsagia circumcincta | L4           | 7    |
| DP_8               | Teladorsagia circumcincta | adult male   | 28   |
| DP_9               | Teladorsagia circumcincta | adult female | 28   |


```bash
mkdir RAW_DATA

# copied Moredun data from SFTP, into a directory called "Tcirc_DNBseq"

cd Tcirc_DNBseq 

# collate RNAseq read files
find . -name "*.gz" -print0 | xargs -0 -I {} mv {} .

```




## Mask the genome
```bash
module load tetools/1.1-c3

cp ../ASSEMBLY/FINAL_GENOME/teladorsagia_circumcincta_tci2_wsi1.0.fa .

bsub.py --queue long --threads 16 20 rm_build "BuildDatabase -name tcirc -engine ncbi teladorsagia_circumcincta_tci2_wsi1.0.fa"

#  RepeatModeler v2.0.1
bsub.py --queue long --threads 16 20 rm_modeller "RepeatModeler -database tcirc -engine ncbi -pa 16"

#Run Repeatmasker
bsub.py --queue long --threads 16 20 rm_modeller "RepeatMasker -pa 16 -gff -nolow -xsmall -lib RM_47333.WedMay312009502023/consensi.fa.classified teladorsagia_circumcincta_tci2_wsi1.0.fa"

```


## Map RNAseq reads
```

cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING

module load hisat2/2.1.0--py36pl5.22.0_0


hisat2-build teladorsagia_circumcincta_tci2_wsi2.0.fa tcirc 

for i in `cd ../RAW_DATA/Tcirc_DNBseq ; ls -1 *_1.fq.gz | sed 's/_1.fq.gz//g'`; do

    hisat2 -p 7 -x tcirc  --dta -1 ../RAW_DATA/Tcirc_DNBseq/${i}_1.fq.gz -2 ../RAW_DATA/Tcirc_DNBseq/${i}_2.fq.gz -S ${i}.tcirc.sam

    samtools view --threads 7 -b -o ${i}.tcirc.bam ${i}.tcirc.sam

    samtools sort -m 7G -o ${i}.tcirc_sorted.bam -T ${i}.tcirc_temp --threads 7 ${i}.tcirc.bam;
    done

rm *tcirc.bam *tcirc.sam

```



## Braker3


```bash
module load ISG/singularity/3.10.0

ANNOTATION_DIR=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION

PREFIX=tcirc_braker
WORKDIR=${ANNOTATION_DIR}/${PREFIX}/data
DESTDIR=${ANNOTATION_DIR}/${PREFIX}/output

mkdir -p ${WORKDIR} ${DESTDIR}

export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output,"
export BRAKER_SIF=/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/braker3_singularity/braker3.sif 

#ln -s ../ASSEMBLY/FINAL_GENOME/teladorsagia_circumcincta_tci2_wsi1.0.fa ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa
cp teladorsagia_circumcincta_tci2_wsi1.0.fa.masked ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked

# w RNAseq test - mapping
echo 'singularity exec \
${BRAKER_SIF} braker.pl \
--species=teladorsagia_circumcincta \
--genome=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked \
--workingdir=${WORKDIR} \
--useexisting \
--rnaseq_sets_ids=DP_1 \
--rnaseq_sets_dirs=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RAW_DATA/Tcirc_DNBseq' 
```
- put all of above in the script "run_braker.sh"

```bash
# run run_braker.sh
bsub.py --queue long --threads 7 20 braker "bash ./run_braker.sh"
```


# path to bam files
braker.pl --species=yourSpecies --genome=genome.fasta \
       --rnaseq_sets_ids=BAM_ID1,BAM_ID2 \
       --rnaseq_sets_dirs=/path/to/local/bam/files/

# actual bam files
braker.pl --species=yourSpecies --genome=genome.fasta \
       --bam=file1.bam,file2.bam


# RNAseq data from Dan Price, mapped
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_10.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_11.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_12.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_13.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_14.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_15.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_1.tcirc.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_2.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_3.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_4.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_5.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_6.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_7.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_8.tcirc_sorted.bam
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_9.tcirc_sorted.bam


```bash
module load ISG/singularity/3.10.0

ANNOTATION_DIR=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION

PREFIX=tcirc_braker_bams
WORKDIR=${ANNOTATION_DIR}/${PREFIX}/data
DESTDIR=${ANNOTATION_DIR}/${PREFIX}/output

mkdir -p ${WORKDIR} ${DESTDIR}

export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output,"
export BRAKER_SIF=/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/braker3_singularity/braker3.sif 

#ln -s ../ASSEMBLY/FINAL_GENOME/teladorsagia_circumcincta_tci2_wsi1.0.fa ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa
cp teladorsagia_circumcincta_tci2_wsi1.0.fa.masked ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked

# w RNAseq test - mapping
singularity exec \
${BRAKER_SIF} braker.pl \
--species=teladorsagia_circumcincta2 \
--genome=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked \
--workingdir=${WORKDIR} \
--useexisting \
--threads=8 \
--bam=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_10.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_11.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_12.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_13.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_14.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_15.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_1.tcirc.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_2.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_3.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_4.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_5.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_6.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_7.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_8.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_9.tcirc_sorted.bam 
```
- put all of above in the script "run_braker_bams.sh"

```bash
# run run_braker.sh
bsub.py --queue long --threads 8 20 braker_bams "bash ./run_braker_bams.sh"
```






# proteins and bams
braker.pl --genome=genome.fa --prot_seq=orthodb.fa \
    --bam=/path/to/SRA_ID1.bam,/path/to/SRA_ID2.bam


# Haemonchus proteins, from WBPS
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.protein.fa.gz

gunzip haemonchus_contortus.PRJEB506.WBPS18.protein.fa.gz




```bash
module load ISG/singularity/3.10.0

ANNOTATION_DIR=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION

PREFIX=tcirc_braker_bams_prot
WORKDIR=${ANNOTATION_DIR}/${PREFIX}/data
DESTDIR=${ANNOTATION_DIR}/${PREFIX}/output

mkdir -p ${WORKDIR} ${DESTDIR}

export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output,"
export BRAKER_SIF=/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/braker3_singularity/braker3.sif 

#ln -s ../ASSEMBLY/FINAL_GENOME/teladorsagia_circumcincta_tci2_wsi1.0.fa ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa
cp teladorsagia_circumcincta_tci2_wsi1.0.fa.masked ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked

# w RNAseq , proteins
singularity exec \
${BRAKER_SIF} braker.pl \
--species=teladorsagia_circumcincta3 \
--genome=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked \
--workingdir=${WORKDIR} \
--useexisting \
--threads=8 \
--prot_seq=/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/haemonchus_contortus.PRJEB506.WBPS18.protein.fa \
--bam=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_10.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_11.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_12.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_13.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_14.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_15.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_1.tcirc.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_2.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_3.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_4.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_5.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_6.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_7.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_8.tcirc_sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_9.tcirc_sorted.bam 
```
- put all of above in the script "run_braker_bams_prot.sh"

```bash
# run run_braker.sh
bsub.py --queue long --threads 8 20 braker_bams_prot "bash ./run_braker_bams_prot.sh"
```





module load ISG/singularity/3.10.0

ANNOTATION_DIR=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION

PREFIX=test
WORKDIR=${ANNOTATION_DIR}/${PREFIX}/data
DESTDIR=${ANNOTATION_DIR}/${PREFIX}/output

mkdir -p ${WORKDIR} ${DESTDIR}

export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output,"
export BRAKER_SIF=/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/braker3_singularity/braker3.sif 

#ln -s ../ASSEMBLY/FINAL_GENOME/teladorsagia_circumcincta_tci2_wsi1.0.fa ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa
cp teladorsagia_circumcincta_tci2_wsi1.0.fa.masked ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked

# w RNAseq test - mapping
singularity exec \
${BRAKER_SIF} braker.pl \
--species=teladorsagia_circumcincta4 \
--genome=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked \
--workingdir=${WORKDIR} \
--useexisting \
--threads=8 \
--prot_seq=/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/haemonchus_contortus.PRJEB506.WBPS18.protein.fa2 \
--bam=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_1.tcirc_n-sorted.bam






```bash
module load ISG/singularity/3.10.0

ANNOTATION_DIR=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION

PREFIX=tcirc_braker_prot
WORKDIR=${ANNOTATION_DIR}/${PREFIX}/data
DESTDIR=${ANNOTATION_DIR}/${PREFIX}/output

mkdir -p ${WORKDIR} ${DESTDIR}

export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output,"
export BRAKER_SIF=/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/braker3_singularity/braker3.sif 

#ln -s ../ASSEMBLY/FINAL_GENOME/teladorsagia_circumcincta_tci2_wsi1.0.fa ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa
cp teladorsagia_circumcincta_tci2_wsi1.0.fa.masked ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked

# w  proteins
singularity exec \
${BRAKER_SIF} braker.pl \
--species=teladorsagia_circumcincta3 \
--genome=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked \
--workingdir=${WORKDIR} \
--useexisting \
--gff3 \
--threads=8 \
--prot_seq=/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/haemonchus_contortus.PRJEB506.WBPS18.protein.fa
```
- put all of above in the script "run_braker_bams_prot.sh"

```bash
# run run_braker.sh
bsub.py --queue long --threads 8 20 braker_prot "bash ./run_braker_prot.sh"
```




```bash
module load ISG/singularity/3.10.0

ANNOTATION_DIR=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION

PREFIX=tcirc_braker_prot_rna-all
WORKDIR=${ANNOTATION_DIR}/${PREFIX}/data
DESTDIR=${ANNOTATION_DIR}/${PREFIX}/output

mkdir -p ${WORKDIR} ${DESTDIR}

export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output,"
export BRAKER_SIF=/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/braker3_singularity/braker3.sif 

#ln -s ../ASSEMBLY/FINAL_GENOME/teladorsagia_circumcincta_tci2_wsi1.0.fa ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa
cp teladorsagia_circumcincta_tci2_wsi1.0.fa.masked ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked

# w  proteins
singularity exec \
${BRAKER_SIF} braker.pl \
--species=teladorsagia_circumcincta3 \
--genome=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked \
--workingdir=${WORKDIR} \
--useexisting \
--gff3 \
--threads=8 \
--prot_seq=/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/haemonchus_contortus.PRJEB506.WBPS18.protein.fa2 \
--bam=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_10.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_11.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_12.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_13.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_14.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_15.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_1.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_2.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_3.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_4.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_5.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_6.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_7.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_8.tcirc_n-sorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_9.tcirc_n-sorted.bam
```
- put all of above in the script "run_braker_complete.sh"

```bash
# run run_braker.sh
bsub.py --queue hugemem --threads 8 200 braker_complete "bash ./run_braker_complete.sh"
```




/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/TSEBRA/bin/tsebra.py \
--gtf /lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/tcirc_braker_prot_rna-all/data/augustus.hints.gtf,/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/tcirc_braker_prot_rna-all/data/GeneMark-ETP/genemark.gtf --keep_gtf /lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/tcirc_braker_prot_rna-all/data/GeneMark-ETP/training.gtf --hintfiles /lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/tcirc_braker_prot_rna-all/data/hintsfile.gff --filter_single_exon_genes --cfg /nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/TSEBRA/bin/../config/braker3.cfg --out /lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/tcirc_braker_prot_rna-all/data/braker.gtf -q 2>/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/tcirc_braker_prot_rna-all/data/errors/tsebra.stderr




```bash
conda activate busco_5.4.3

bsub.py --queue yesterday --threads 1 10 busco_nematoda_odb10 \
    "busco --in braker.aa --out BUSCO_nematoda_odb10 --mode protein --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 -f -r"



module load hisat2/2.1.0--py36pl5.22.0_0
for i in ` ls -1 *_1.fastq.gz | sed 's/_1.fastqq.gz//g'`; do

    hisat2 -p 7 -x tcirc  --dta -1 ${i}_1.fastq.gz -2 ${i}_2.fastq.gz -S ${i}.tcirc.sam

    samtools view --threads 7 -b -o ${i}.tcirc.bam ${i}.tcirc.sam

    samtools sort -m 7G -o ${i}.tcirc_sorted.bam -T ${i}.tcirc_temp --threads 7 ${i}.tcirc.bam;
    done

rm *tcirc.bam *tcirc.sam
```





# merged RNAseq data

```bash
module load ISG/singularity/3.10.0

ANNOTATION_DIR=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION

PREFIX=tcirc_braker_rna_new
WORKDIR=${ANNOTATION_DIR}/${PREFIX}/data
DESTDIR=${ANNOTATION_DIR}/${PREFIX}/output

mkdir -p ${WORKDIR} ${DESTDIR}

export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output,"
export BRAKER_SIF=/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/braker3_singularity/braker3.sif 

#ln -s ../ASSEMBLY/FINAL_GENOME/teladorsagia_circumcincta_tci2_wsi1.0.fa ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa
cp teladorsagia_circumcincta_tci2_wsi1.0.fa.masked ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked

# w 
singularity exec \
${BRAKER_SIF} braker.pl \
--genome=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked \
--workingdir=${WORKDIR} \
--gff3 \
--threads=8 \
--bam=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/merged.n-storted.bam
```

```bash
module load ISG/singularity/3.10.0

ANNOTATION_DIR=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION

PREFIX=tcirc_braker_rna_old
WORKDIR=${ANNOTATION_DIR}/${PREFIX}/data
DESTDIR=${ANNOTATION_DIR}/${PREFIX}/output

mkdir -p ${WORKDIR} ${DESTDIR}

export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output,"
export BRAKER_SIF=/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/braker3_singularity/braker3.sif 

#ln -s ../ASSEMBLY/FINAL_GENOME/teladorsagia_circumcincta_tci2_wsi1.0.fa ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa
cp teladorsagia_circumcincta_tci2_wsi1.0.fa.masked ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked

# w 
singularity exec \
${BRAKER_SIF} braker.pl \
--genome=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked \
--workingdir=${WORKDIR} \
--gff3 \
--threads=8 \
--bam=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/merged.old-storted.bam
```









```bash
module load ISG/singularity/3.10.0

ANNOTATION_DIR=/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ASSEMBLY/FINAL_GENOME/RE_POLCA/POLCA_R2

PREFIX=tcirc_braker_prot
WORKDIR=${ANNOTATION_DIR}/${PREFIX}/data
DESTDIR=${ANNOTATION_DIR}/${PREFIX}/output

mkdir -p ${WORKDIR} ${DESTDIR}

export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output,"
export BRAKER_SIF=/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/braker3_singularity/braker3.sif 


# w  proteins
singularity exec \
${BRAKER_SIF} braker.pl \
--species=teladorsagia_circumcincta3 \
--genome=/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ASSEMBLY/FINAL_GENOME/RE_POLCA/POLCA_R2/teladorsagia_circumcincta_tci2_wsi1.0.fa.PolcaCorrected.fa.PolcaCorrected.fa \
--workingdir=${WORKDIR} \
--useexisting \
--gff3 \
--threads=8 \
--prot_seq=/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/haemonchus_contortus.PRJEB506.WBPS18.protein.fa
```
- put all of above in the script "run_braker_bams_prot.sh"

```bash
# run run_braker.sh
bsub.py --queue long --threads 8 20 braker_prot "bash ./run_braker_proteins-only.sh"
```


```bash
module load ISG/singularity/3.10.0

ANNOTATION_DIR=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION

PREFIX=tcirc_braker_prot_rna-all
WORKDIR=${ANNOTATION_DIR}/${PREFIX}/data
DESTDIR=${ANNOTATION_DIR}/${PREFIX}/output

mkdir -p ${WORKDIR} ${DESTDIR}

export SINGULARITY_BIND="/nfs:/nfs,/lustre:/lustre,$WORKDIR:/data,$DESTDIR:/output,"
export BRAKER_SIF=/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/braker3_singularity/braker3.sif

#ln -s ../ASSEMBLY/FINAL_GENOME/teladorsagia_circumcincta_tci2_wsi1.0.fa ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa
cp teladorsagia_circumcincta_tci2_wsi1.0.fa.masked ${WORKDIR}/teladorsagia_circumcincta_tci2_wsi1.0.fa.masked

# w  proteins
singularity exec \
${BRAKER_SIF} braker.pl \
--genome=/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ASSEMBLY/FINAL_GENOME/RE_POLCA/POLCA_R2/teladorsagia_circumcincta_tci2_wsi2.0.fa \
--workingdir=${WORKDIR} \
--useexisting \
--gff3 \
--threads=20 \
--prot_seq=/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/haemonchus_contortus.PRJEB506.WBPS18.protein.fa2 \
--bam=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_10.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_11.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_12.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_13.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_14.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_15.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_1.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_2.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_3.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_4.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_5.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_6.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_7.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_8.tcirc_nsorted.bam,/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/RNASEQ_MAPPING/DP_9.tcirc_nsorted.bam
```






```bash
# combining braker and augustus intermediate GFFs

braker.gff3
augustus.hints.gff3


# identify genes present in the augustus data that is missing from the braker data
bedtools subtract -s -a augustus.hints.gff3 -b braker.gff3 > augustus_genes_missing_from_braker.gff3

# remove single exon genes
gffread -U augustus_genes_missing_from_braker.gff3 -o augustus_genes_missing_from_braker.filtered.gff3



# tcirc 2.4 genome
teladorsagia_circumcincta_tci2_wsi2.4.fa




# filter the gffs to include only contigs/scaffolds present in the new tc2.4 assembly
grep ">" teladorsagia_circumcincta_tci2_wsi2.4.fa | sed 's/>//g' > teladorsagia_circumcincta_tci2_wsi2.4.names

grep -w -f teladorsagia_circumcincta_tci2_wsi2.4.names augustus_genes_missing_from_braker.filtered.gff3 > augustus_genes_missing_from_braker.filtered.tc2.4.gff3

grep -w -f teladorsagia_circumcincta_tci2_wsi2.4.names braker.gff3 > braker.tc2.4.gff3


# combine braker and augustus gffs
cat augustus_genes_missing_from_braker.filtered.tc2.4.gff3 braker.tc2.4.gff3 | sort -k1,1 -k4,4n > braker_augustus.tc2.4.gff



# make proteins from each gff3 iteration
gffread braker.gff3 -g ../teladorsagia_circumcincta_tci2_wsi2.0.fa -y braker.proteins.fa
gffread braker.tc2.4.gff3 -g teladorsagia_circumcincta_tci2_wsi2.4.fa -y braker.tc2.4.proteins.fa
gffread braker_augustus.tc2.4.gff -g teladorsagia_circumcincta_tci2_wsi2.4.fa -y braker_augustus.tc2.4.proteins.fa





bsub.py --queue long --threads 20 40 busco_braker "busco --in braker.proteins.fa --out braker.proteins --mode protein --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"

bsub.py --queue long --threads 20 40 busco_braker_augustus.tc2.4.proteins "busco --in braker_augustus.tc2.4.proteins.fa --out braker_augustus.tc2.4.proteins --mode protein --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"

bsub.py --queue long --threads 20 40 busco_braker.tc2.4.proteins "busco --in braker.tc2.4.proteins.fa --out braker.tc2.4.proteins --mode protein --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"

bsub.py --queue long --threads 20 40 busco_hon.proteins "busco --in haemonchus_contortus.PRJEB506.WBPS18.protein.fa --out hcon.proteins --mode protein --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"
```


- braker 
    - C:89.1%[S:42.7%,D:46.4%],F:0.7%,M:10.2%,n:3131

- braker 2.4
    - C:88.7%[S:43.6%,D:45.1%],F:0.6%,M:10.7%,n:3131

- braker + augustus 2.4
    - C:95.3%[S:46.9%,D:48.4%],F:1.0%,M:3.7%,n:3131

- Haemonchus
    - C:96.2%[S:85.1%,D:11.1%],F:0.5%,M:3.3%,n:3131




```bash
# clean up the braker gff
gt gff3 -force -tidy -setsource braker -o braker.tc2.4.gt.gff3 braker.tc2.4.gff3
#> gene = 16694
#> mRNA = 30031








# extract mRNAs and convert them to genes for augustus gff - they are missing.


awk -F '[\t;]' '{if($3=="mRNA") print $10}' augustus_genes_missing_from_braker.filtered.tc2.4.gff3 | sort | uniq |  while read -r ID; do 
    grep -w $ID augustus_genes_missing_from_braker.filtered.tc2.4.gff3 |\
    sort -k4,4 -k5,5 > tmp ;\
    start=$(head -n1 tmp | awk '{print $4}');
    end=$(tail -n1 tmp | awk '{print $5}');
    grep -w $ID augustus_genes_missing_from_braker.filtered.tc2.4.gff3 | awk -v name=$ID -v start=$start -v end=$end '{print $1,$2,"gene",start,end,1,$7,".",name}' OFS="\t";
    done | sort | uniq | sed 's/geneID/ID/g' > augustus_genes.genes.gff

# convert "geneID" to "Parent" in mRNA rows
sed -i 's/geneID/Parent/g'    augustus_genes_missing_from_braker.filtered.tc2.4.gff3

# merge mRNA and progeny with gene IDs
cat augustus_genes_missing_from_braker.filtered.tc2.4.gff3 augustus_genes.genes.gff > tmp; mv tmp augustus_genes_missing_from_braker.filtered.tc2.4.gff3 

# clean up
gt gff3 -force -tidy -setsource augustus -o augustus_genes_missing_from_braker.filtered.tc2.4.gt.gff3 augustus_genes_missing_from_braker.filtered.tc2.4.gff3


#> gene = 7024
#> mRNA = 8271

# making the genesets unique within braker and augustus gffs, before merging

cat braker.tc2.4.gt.gff3 |\
    sed -e 's/=mRNA/=b.mRNA/g' -e 's/=gene/=b.gene/g' > braker.tc2.4.gt.renamed.gff3 

cat augustus_genes_missing_from_braker.filtered.tc2.4.gt.gff3 |\
    sed -e 's/=mRNA/=a.mRNA/g' -e 's/=gene/=a.gene/g' > augustus_genes_missing_from_braker.filtered.tc2.4.gt.renamed.gff3

gt gff3 -tidy -sort -force -o braker_augustus.2.4.gt.gff3 braker.tc2.4.gt.renamed.gff3 augustus_genes_missing_from_braker.filtered.tc2.4.gt.renamed.gff3


 437385 CDS
 437386 exon
  23718 gene
 327667 intron
  38302 mRNA
  30022 start_codon
  30023 stop_codon


  # integrating apollo manual curation
  #- identify overlapping positions of curated genes in current reference, and remove them
  #- format new genes
  #- merge

# download and transfer apollo annotations
  scp Annotations.gff3.gz sd21@farm5-head2:/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/BRAKER_V_AUGUSTUS/
  gunzip Annotations.gff3.gz


cat Annotations.gff3 | cut -f3 | sort | uniq -c | grep -v "#"
  15631 CDS
  15847 exon
   1209 gene
   1342 mRNA
    201 non_canonical_five_prime_splice_site
     64 non_canonical_three_prime_splice_site
      1 pseudogene
      2 rRNA
      1 transcript


grep -w -f teladorsagia_circumcincta_tci2_wsi2.4.names Annotations.gff3  > Annotations.tc2.4.gff3 

gt gff3 -force -tidy -sort -o Annotations.tc2.4.gt.gff3 Annotations.tc2.4.gff3 


# identify genes present in the apollo data that is are present in the braker data

bedtools subtract -A -s -a braker_augustus.2.4.gt.gff3 -b Annotations.gt.gff3 | awk '{if($3=="gene") print $9}' | sed 's/ID=//g' > gene_ids_to_keep.list

#--- TEST FIX
bedtools subtract -A -s -a braker_augustus.2.4.gt.gff3 -b Annotations.tc2.4.gt.gff3 | awk '{if($3=="gene") print $9}' | sed 's/ID=//g' > test.gene_ids_to_keep.list


# filter the data
conda activate agat

/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/AGAT/bin/agat_sp_filter_feature_from_keep_list.pl --gff braker_augustus.2.4.gt.gff3 --keep_list gene_ids_to_keep.list --output braker_augustus.2.4.minus.curated.gff3

#--- TEST FIX
/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/AGAT/bin/agat_sp_filter_feature_from_keep_list.pl --gff braker_augustus.2.4.gt.gff3 --keep_list test.gene_ids_to_keep.list --output test.braker_augustus.2.4.minus.curated.gff3



# add start and stops
/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/AGAT/bin/agat_sp_add_start_and_stop.pl --gff braker_augustus.2.4.minus.curated.gff3 --fasta teladorsagia_circumcincta_tci2_wsi2.4.fa --out braker_augustus.2.4.minus.curated_start-stop-fix.gff3





# merge braker+augustus with apollo
gt gff3 -tidy -sort -force -o braker_augustus.apollo.2.4.gff3 braker_augustus.2.4.minus.curated_start-stop-fix.gff3 Annotations.tc2.4.gt.gff3

#--- TEST FIX
gt gff3 -tidy -sort -force -o test.braker_augustus.apollo.2.4.gff3 test.braker_augustus.2.4.minus.curated.gff3 Annotations.tc2.4.gt.gff3







# sort
cat braker_augustus.apollo.2.4.gff3 | sort -k1,1 -k4,4n > braker_augustus.apollo.2.4.sorted.gff3

cat braker_augustus.apollo.2.4.sorted.gff3 | cut -f3 | sort | uniq -c | grep -v "#"

 436351 CDS
 436555 exon
  24010 gene
 318556 intron
  38163 mRNA
    191 non_canonical_five_prime_splice_site
     62 non_canonical_three_prime_splice_site
      1 pseudogene
      2 rRNA
  36227 start_codon
  35886 stop_codon
      1 transcript

/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/AGAT/bin/agat_sp_add_start_and_stop.pl --gff braker_augustus.apollo.2.4.sorted.gff3 --fasta teladorsagia_circumcincta_tci2_wsi2.4.fa --out tmp.braker_augustus.apollo.2.4.sorted.gff3

gt gff3 -tidy -sort tmp.braker_augustus.apollo.2.4.sorted.gff3 | grep -v "###" | sort -k1,1 -k4,4n > braker_augustus.apollo.2.4.sorted.gff3

cat braker_augustus.apollo.2.4.sorted.gff3 | cut -f3 | sort | uniq -c | grep -v "#"

 436351 CDS
 436555 exon
    715 five_prime_UTR
  24010 gene
 318556 intron
  38163 mRNA
    191 non_canonical_five_prime_splice_site
     62 non_canonical_three_prime_splice_site
      1 pseudogene
      2 rRNA
  37440 start_codon
  37162 stop_codon
   1117 three_prime_UTR
      1 transcript

# compress for viewing in apollo
bgzip braker_augustus.apollo.2.4.sorted.gff3
tabix braker_augustus.apollo.2.4.sorted.gff3.gz
```



```bash
# check protein completeness of original, and filtered for longest isoform
#- generate protein sequences from gff
grep -v "rRNA" braker_augustus.apollo.2.4.sorted.gff3 | gffread - -g teladorsagia_circumcincta_tci2_wsi2.4.fa -y braker_augustus.apollo.2.4.proteins.fa

#-extract longest isoforms from GFF, and then create proteins
conda activate agat
/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/AGAT/bin/agat_sp_keep_longest_isoform.pl -gff braker_augustus.apollo.2.4.sorted.gff3  -out braker_augustus.apollo.2.4.longest_isoform.gff3

gffread braker_augustus.apollo.2.4.longest_isoform.gff3 -g teladorsagia_circumcincta_tci2_wsi2.4.fa -y braker_augustus.apollo.2.4.longest_isoform.proteins.fa

#- run busco
conda activate busco_5.4.3
bsub.py --queue long --threads 20 40 busco_combined "busco --in braker_augustus.apollo.2.4.proteins.fa --out braker_augustus.apollo.2.4.proteins --mode protein --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"

bsub.py --queue long --threads 20 40 busco_combined_longest "busco --in braker_augustus.apollo.2.4.longest_isoform.proteins.fa --out braker_augustus.apollo.2.4.longest_isoform.proteins --mode protein --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


# braker + augustus + apollo 
C:95.7%[S:47.2%,D:48.5%],F:0.9%,M:3.4%,n:3131

# braker + augustus + apollo - longest isoform
C:95.7%[S:81.2%,D:14.5%],F:0.9%,M:3.4%,n:3131

```




```bash
# renaming GFF

# STEP 1 - generate using IDs for each gene
# - extract gene name from lines in gff matching gene|GENE
# - remove surrounding characters
# - sort by gene ID and remove duplicates
# - print new gene id (species prefix with 8 digit ID that increases by 10 for each gene), and old gene id

gff=braker_augustus.apollo.2.4.sorted.gff3
species_prefix=TCIR_

awk -F'[\t;]' '{if($3=="gene" || $3=="GENE") print $9}' ${gff} | sed -e 's/ID=//g' -e 's/\;//g' | sort -V | uniq | awk -v species_prefix="$species_prefix" '{fmt=species_prefix"%08d\t%s\n"; printf fmt,NR*10,$0}' > genes_renames.list


# correct gene IDs
cp ${gff} ${gff}.tmp

while read new_gene old_gene; do 
    sed -i \
    -e "s/ID=${old_gene}$/ID=${new_gene}/g" \
    -e "s/ID=${old_gene}\;/ID=${new_gene}\;/g" \
    -e "s/Parent=${old_gene}$/Parent=${new_gene}/g" \
    -e "s/Parent=${old_gene}\;/Parent=${new_gene}\;/g" \
    ${gff}.tmp; 
done < genes_renames.list


# fix mRNA IDs
cat ${gff}.tmp | grep "mRNA" | awk -F '[\t;]' '{if($3=="mRNA") print $9,$10}' OFS="\t" | sed -e 's/ID=//g' -e 's/Parent=//g' > mRNA_gene_IDs.list

#- add transcript identifiers
>mRNA_IDs_NAMEs_transcriptIDs.txt
cut -f2 mRNA_gene_IDs.list | sort | uniq | while read -r GENE; do 
    grep -w ${GENE} mRNA_gene_IDs.list | cat -n | awk '{print $2,$3,$3"-0000"$1}' OFS="\t" >> mRNA_IDs_NAMEs_transcriptIDs.txt; 
    done

# substitute mRNA ids
while read mRNA_id gene_id transcript_id; do 
    sed -i \
    -e "s/ID=${mRNA_id}$/ID=${transcript_id}/g" \
    -e "s/ID=${mRNA_id}\;/ID=${transcript_id}\;/g" \
    -e "s/Parent=${mRNA_id}$/Parent=${transcript_id}/g" \
    -e "s/Parent=${mRNA_id}\;/Parent=${transcript_id}\;/g" \
    ${gff}.tmp; 
done < mRNA_IDs_NAMEs_transcriptIDs.txt



mv braker_augustus.apollo.2.4.sorted.gff3.tmp braker_augustus.apollo.2.4.renamed.gff3

bgzip braker_augustus.apollo.2.4.renamed.gff3
tabix braker_augustus.apollo.2.4.renamed.gff3.gz



conda activate agat
/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/AGAT/bin/agat_sp_manage_attributes.pl -gff braker_augustus.apollo.2.4.renamed.gff3 --tag all_attributes -p level1,level2,level2 --out braker_augustus.apollo.2.4.clean.gff3

cat braker_augustus.apollo.2.4.clean.gff3 | sort -k1,1 -k4,4n | bgzip > braker_augustus.apollo.2.4.clean.gff3.gz
tabix braker_augustus.apollo.2.4.clean.gff3.gz


/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/AGAT/bin/agat_sp_webApollo_compliant.pl --gff braker_augustus.apollo.2.4.clean.gff3 --output braker_augustus.apollo.2.4.apollo.gff3
```





## interproscan
```bash
module load interproscan/5.57-90.0

### interproscan

```bash
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ANNOTATION/INTERPROSCAN

ln -s ../BRAKER_V_AUGUSTUS/braker_augustus.apollo.2.4.clean.gff3
ln -s ../BRAKER_V_AUGUSTUS/teladorsagia_circumcincta_tci2_wsi2.4.fa

# generate a protein fasta from annotation and reference genome
gffread -y PROTEINS.fa -g teladorsagia_circumcincta_tci2_wsi2.4.fa braker_augustus.apollo.2.4.clean.gff3
sed -e 's/\.//g' PROTEINS.fa > tmp; mv tmp PROTEINS.fa

# run interproscan
farm_interproscan -a PROTEINS.fa -o IPS.output.gff

# lift over GO terms from interproscan to GFF
extract_interproscan_go_terms -i IPS.output.gff -e braker_augustus.apollo.2.4.clean.gff3

# filter and rename
grep ^'\#\#\|hc' tmp.gff.go.gff | grep -v "mtDNA" | grep -v "FASTA" | grep -v ">" > HCON_V4_WBP11plus_190125.ips.gff3
```


```bash
fastaq split_by_base_count PROTEINS.fa tc_proteins_ --max_seqs 100 10000000


mkdir IPS_OUT

for i in tc_proteins*; do
interproscan.sh  --input ${i} --goterms --iprlookup --output-dir IPS_OUT; 
done

--applications

bsub.py --queue basement --threads 30 10 ips "interproscan.sh  --input PROTEINS.fa --cpu 30 --goterms --iprlookup --output-dir IPS_OUT"

conda activate agat

cat IPS_OUT/*tsv > ips_collated.tsv

agat_sp_manage_functional_annotation.pl -f braker_augustus.apollo.2.4.clean.gff3 -i ips_collated.tsv --output braker_augustus.apollo.2.4.ips
```







```bash
####----- TESTING


cat test.gene_ids_to_keep.list |\
    while read geneID; do 
        grep "${geneID}$" braker_augustus.2.4.gt.gff3 |\
        grep "mRNA" | cut -f9 | sed -e 's/ID=//g' -e 's/;Parent=/\t/g' ; 
    done > test.orignal.mrna.gene.list


mkdir kept_gene_models

while read mrna gene; do
    > kept_gene_models/${gene}_${mrna}.model.gff
    #gene
    grep "ID=${gene}$" test.braker_augustus.2.4.minus.curated.gff3 >> kept_gene_models/${gene}_${mrna}.model.gff
    #mrna
    grep "ID=${mrna};" test.braker_augustus.2.4.minus.curated.gff3 >> kept_gene_models/${gene}_${mrna}.model.gff
    #body
    grep "${mrna}$" test.braker_augustus.2.4.minus.curated.gff3 >> kept_gene_models/${gene}_${mrna}.model.gff

done < test.orignal.mrna.gene.list




cat kept_gene_models/* |\
    sed -e 's/=mRNA/=b.mRNA/g' -e 's/=gene/=b.gene/g' > kept_models.gff3 

cat Annotations.tc2.4.gt.gff3 |\
    sed -e 's/=mRNA/=a.mRNA/g' -e 's/=gene/=a.gene/g' > updated.models.gff3

cat kept_models.gff3 updated.models.gff3 | sort | uniq | sort -k1,1 -k4,4n > test.new_annotation.gff

gt gff3 -tidy -sort -force -o test.new_annotation.gt.gff test.new_annotation.gff

cat test.new_annotation.gt.gff | sort -k1,1 -k4,4n | bgzip > test.new_annotation.gt.gff.gz


#- generate protein sequences from gff
agat_sp_extract_sequences.pl -g test.new_annotation.gt.gff -f teladorsagia_circumcincta_tci2_wsi2.4.fa -p --output test.new_annotation.gt.proteins.fa




# finding genes with internal stop codons
fastaq to_fasta -l 0 test.new_annotation.gt.proteins.fa test.new_annotation.gt.proteins.l0.fa

cat test.new_annotation.gt.proteins.l0.fa | sed 's/\*$//g' | grep -A1 "\*" | grep ">" > genes_w_stops.list

cat genes_w_stops.list | sed 's/>//g' | while read MRNA GENE SEQ TYPE; do 
    grep "ID=${MRNA};" test.new_annotation.gt.gff; 
    done



# download and transfer apollo annotations
scp Annotations_R2.gff3.gz sd21@farm5-head2:/lustre/scratch125/pam/teams/team333/sd21/teladorsagia_circumcincta/GENOME/ANNOTATION/BRAKER_V_AUGUSTUS/
gunzip Annotations_R2.gff3.gz


cat Annotations_R2.gff3 | cut -f3 | sort | uniq -c | grep -v "#"
   6276 CDS
   6446 exon
    580 gene
    631 mRNA
     88 non_canonical_five_prime_splice_site
     35 non_canonical_three_prime_splice_site


grep -w -f teladorsagia_circumcincta_tci2_wsi2.4.names Annotations_R2.gff3  > Annotations_R2.tc2.4.gff3 

gt gff3 -force -tidy -sort -o Annotations_R2.tc2.4.gt.gff3 Annotations_R2.tc2.4.gff3

cat Annotations_R2.tc2.4.gt.gff3 | cut -f3 | sort | uniq -c | grep -v "#"
   5857 CDS
   6023 exon
    541 gene
    590 mRNA
     82 non_canonical_five_prime_splice_site
     34 non_canonical_three_prime_splice_site


bedtools subtract -A -s -a test.new_annotation.gt.gff -b Annotations_R2.tc2.4.gt.gff3 | awk -F '[\t;]' '{if($3=="gene") print $9}' | sed 's/ID=//g' > R2.gene_ids_to_keep.list

# extract the gene coords from both annotations, and then compare with bedtools subtract
# cat test.new_annotation.gt.gff | awk '{if($3=="gene") print}' > test.new_annotation.gt.gff.genes
# cat Annotations_R2.tc2.4.gt.gff3  | awk '{if($3=="gene") print}' > Annotations_R2.tc2.4.gt.gff3.genes
```


```bash
wc -l R2.gene_ids_to_keep.list
22963 R2.gene_ids_to_keep.list


cat R2.gene_ids_to_keep.list |\
    while read geneID; do 
        grep "${geneID}$" test.new_annotation.gt.gff |\
        grep "mRNA" | cut -f9 | sed -e 's/ID=//g' -e 's/;Parent=/\t/g' ; 
    done > R2.orignal.mrna.gene.list


mkdir R2.kept_gene_models

while read mrna gene; do
    > R2.kept_gene_models/${gene}_${mrna}.model.gff
    #gene
    grep "ID=${gene}$" test.new_annotation.gt.gff >> R2.kept_gene_models/${gene}_${mrna}.model.gff
    #mrna
    grep "ID=${mrna};" test.new_annotation.gt.gff >> R2.kept_gene_models/${gene}_${mrna}.model.gff
    #body
    grep "${mrna}$" test.new_annotation.gt.gff >> R2.kept_gene_models/${gene}_${mrna}.model.gff

done < R2.orignal.mrna.gene.list


cat R2.kept_gene_models/* |\
    sed -e 's/=mRNA/=b.mRNA/g' -e 's/=gene/=b.gene/g' > R2.kept_models.gff3 

cat Annotations_R2.tc2.4.gt.gff3 |\
    sed -e 's/=mRNA/=a.mRNA/g' -e 's/=gene/=a.gene/g' > R2.updated.models.gff3

cat R2.kept_models.gff3 R2.updated.models.gff3 | sort | uniq | sort -k1,1 -k4,4n > R2.new_annotation.gff

gt gff3 -tidy -sort -force -o R2.new_annotation.gt.gff R2.new_annotation.gff

cat R2.new_annotation.gt.gff | sort -k1,1 -k4,4n | bgzip > R2.new_annotation.gt.gff.gz

tabix R2.new_annotation.gt.gff.gz

cut -f3 R2.new_annotation.gt.gff | grep -v "#" | sort | uniq -c
 413747 CDS
 413913 exon
  22399 gene
 308653 intron
  36182 mRNA
     82 non_canonical_five_prime_splice_site
     34 non_canonical_three_prime_splice_site
  28212 start_codon
  28212 stop_codon



agat_sp_add_start_and_stop.pl --gff R2.new_annotation.gt.gff --fasta teladorsagia_circumcincta_tci2_wsi2.4.fa --out R2.new_annotation.gt.start.stop.gff
```bash




```bash
# renaming GFF

# STEP 1 - generate using IDs for each gene
# - extract gene name from lines in gff matching gene|GENE
# - remove surrounding characters
# - sort by gene ID and remove duplicates
# - print new gene id (species prefix with 8 digit ID that increases by 10 for each gene), and old gene id

gff=R2.new_annotation.gt.start.stop.gff
species_prefix=TCIR_

awk -F'[\t;]' '{if($3=="gene" || $3=="GENE") print $9}' ${gff} | sed -e 's/ID=//g' -e 's/\;//g' | sort -V | uniq | awk -v species_prefix="$species_prefix" '{fmt=species_prefix"%08d\t%s\n"; printf fmt,NR*10,$0}' > genes_renames.list


# correct gene IDs
cp ${gff} ${gff}.tmp

while read new_gene old_gene; do 
    sed -i \
    -e "s/ID=${old_gene}$/ID=${new_gene}/g" \
    -e "s/ID=${old_gene}\;/ID=${new_gene}\;/g" \
    -e "s/Parent=${old_gene}$/Parent=${new_gene}/g" \
    -e "s/Parent=${old_gene}\;/Parent=${new_gene}\;/g" \
    ${gff}.tmp; 
done < genes_renames.list


# fix mRNA IDs
cat ${gff}.tmp | grep "mRNA" | awk -F '[\t;]' '{if($3=="mRNA") print $9,$10}' OFS="\t" | sed -e 's/ID=//g' -e 's/Parent=//g' > mRNA_gene_IDs.list

#- add transcript identifiers
>mRNA_IDs_NAMEs_transcriptIDs.txt
cut -f2 mRNA_gene_IDs.list | sort | uniq | while read -r GENE; do 
    grep -w ${GENE} mRNA_gene_IDs.list | cat -n | awk '{print $2,$3,$3"-0000"$1}' OFS="\t" >> mRNA_IDs_NAMEs_transcriptIDs.txt; 
    done

# substitute mRNA ids
while read mRNA_id gene_id transcript_id; do 
    sed -i \
    -e "s/ID=${mRNA_id}$/ID=${transcript_id}/g" \
    -e "s/ID=${mRNA_id}\;/ID=${transcript_id}\;/g" \
    -e "s/Parent=${mRNA_id}$/Parent=${transcript_id}/g" \
    -e "s/Parent=${mRNA_id}\;/Parent=${transcript_id}\;/g" \
    ${gff}.tmp; 
done < mRNA_IDs_NAMEs_transcriptIDs.txt


agat_sp_manage_attributes.pl -gff R2.new_annotation.gt.start.stop.gff.tmp -p level1,level2 --att ID/Name,date_creation,date_last_modified,orig_id,owner --cp --overwrite --out R2.new_annotation.gt.updatenames.gff.tmp

mv R2.new_annotation.gt.updatenames.gff.tmp braker.augustus.apollo.2.4.2.renamed.gff3








#- generate protein sequences from gff
agat_sp_extract_sequences.pl -g braker.augustus.apollo.2.4.2.renamed.gff3  -f teladorsagia_circumcincta_tci2_wsi2.4.fa -p --output braker.augustus.apollo.2.4.2.renamed.proteins.fa




# finding genes with internal stop codons
fastaq to_fasta -l 0 braker.augustus.apollo.2.4.2.renamed.proteins.fa braker.augustus.apollo.2.4.2.renamed.proteins.l0.fa

cat braker.augustus.apollo.2.4.2.renamed.proteins.l0.fa | sed 's/\*$//g' | grep -A1 "\*" | grep ">" > genes_w_stops.list

cat genes_w_stops.list | sed 's/>//g' | while read MRNA GENE SEQ TYPE; do 
    grep "ID=${MRNA};" braker.augustus.apollo.2.4.2.renamed.gff3; 
    done



cat braker.augustus.apollo.2.4.2.renamed.gff3 |\
    grep -v "TCIR_00076170" |\
    grep -v "TCIR_00192170" |\
    grep -v "TCIR_00192180" > braker.augustus.apollo.2.4.2.renamed.gff3.tmp



cat braker.augustus.apollo.2.4.2.renamed.gff3.tmp | sort -k1,1 -k4,4n | bgzip > braker.augustus.apollo.2.4.2.renamed.gff3.gz
tabix braker.augustus.apollo.2.4.2.renamed.gff3.gz





# generate protein sequences from gff
agat_sp_extract_sequences.pl -g braker.augustus.apollo.2.4.2.renamed.gff3.tmp  -f teladorsagia_circumcincta_tci2_wsi2.4.fa -p --output braker.augustus.apollo.2.4.2.renamed.proteins.fa


# check protein completeness of original, and filtered for longest isoform

#-extract longest isoforms from GFF, and then create proteins

agat_sp_keep_longest_isoform.pl -gff braker.augustus.apollo.2.4.2.renamed.gff3.tmp  -out braker.augustus.apollo.2.4.2.renamed.longest-isoform.gff3

agat_sp_extract_sequences.pl -g braker.augustus.apollo.2.4.2.renamed.longest-isoform.gff3  -f teladorsagia_circumcincta_tci2_wsi2.4.fa -p --output braker.augustus.apollo.2.4.2.renamed.proteins.longest-isoform.fa


#- run busco
conda activate busco_5.4.3
bsub.py --queue long --threads 20 40 busco_combined "busco --in braker.augustus.apollo.2.4.2.renamed.proteins.fa --out braker.augustus.apollo.2.4.2.renamed.proteins --mode protein --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"

bsub.py --queue long --threads 20 40 busco_combined_longest "busco --in braker.augustus.apollo.2.4.2.renamed.proteins.longest-isoform.fa --out braker.augustus.apollo.2.4.2.renamed.proteins.longest-isoform --mode protein --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"

```bash



### compleasm
```bash
conda activate compleasm

braker.augustus.apollo.2.4.2.renamed.proteins.fa

# miniprot
miniprot --trans -u -I --outs=0.95 --gff -t 8 teladorsagia_circumcincta_tci2_wsi2.4.fa braker.augustus.apollo.2.4.2.renamed.proteins.fa > tc_2.4.2.gff

# analysis with miniprot output gff file
bsub.py 10 compleasm_tc_2.4.2 "compleasm analyze -g tc_2.4.2.gff -o compleasm_tc_2.4.2 -l nematoda -t 8 --library_path /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/BUSCO/mb_downloads"

```



```bash
#### ---- TESTING 

mkdir FIX_MISSING_GENES

bedtools intersect -f 0.75 -s -v -a braker_augustus.apollo.2.4.clean.gff3 -b braker.augustus.apollo.2.4.2.renamed.gff3 > FIX_MISSING_GENES/missing_genes.gff

cd FIX_MISSING_GENES

cat missing_genes.gff | awk '{if($3=="gene") print}' | cut -f9 | sed 's/ID=//g' > missing.genes.list

wc -l missing.genes.list
#> 1489 missing.genes.list

while read GENEID; do
    grep "${GENEID}" missing_genes.gff >> missing_genes.filtered.gff;
    done  < missing.genes.list


conda activate agat

agat_sp_fix_overlaping_genes.pl -f missing_genes.filtered.gff -o missing_genes.agat.remove_overlaps.gff

# fix a problematic gene
grep -v "nbis-gene" missing_genes.agat.remove_overlaps.gff | grep -v "TCIR_00153840" > missing_genes.agat.remove_overlaps.no-nbis.gff

# after filter
awk '{if($3=="gene") print}' missing_genes.agat.remove_overlaps.no-nbis.gff | wc -l
#> 1120

# check in apollo - bgzip/tabix
cat missing_genes.agat.remove_overlaps.no-nbis.gff | sort -k1,1 -k4,4n | bgzip > missing_genes.agat.remove_overlaps.gff.gz
tabix missing_genes.agat.remove_overlaps.gff.gz

#> all looked ok


# fix names so no overlaps with exisiting gene IDs
cat missing_genes.agat.remove_overlaps.no-nbis.gff | sed 's/TCIR_0/TCIR_1/g' > missing_genes.agat.remove_overlaps.no-nbis.fix-names.gff


cp ../braker.augustus.apollo.2.4.2.renamed.gff3 .

agat_sp_merge_annotations.pl --gff braker.augustus.apollo.2.4.2.renamed.gff3 --gff missing_genes.agat.remove_overlaps.no-nbis.fix-names.gff --out tcircumcincta.v2.5.gff3.tmp




# finding genes with internal stop codons

# generate protein sequences from gff
agat_sp_extract_sequences.pl -g tcircumcincta.v2.5.gff3.tmp  -f ../teladorsagia_circumcincta_tci2_wsi2.4.fa -p --output tcircumcincta.v2.5.proteins.fa

fastaq to_fasta -l 0 tcircumcincta.v2.5.proteins.fa tcircumcincta.v2.5.proteins.l0.fa

cat tcircumcincta.v2.5.proteins.l0.fa | sed 's/\*$//g' | grep -A1 "\*" | grep ">" > genes_w_stops.list

cat genes_w_stops.list | sed 's/>//g' | while read MRNA GENE SEQ TYPE; do 
    grep "ID=${MRNA};" tcircumcincta.v2.5.gff3.tmp; 
    done

more genes_w_stops.list
>TCIR_00192180-00001 gene=TCIR_00192180 seq_id=tci2_wsi2.0_lg_2.1 type=cds
>TCIR_10038850-00002 gene=TCIR_10038850 seq_id=tci2_wsi2.0_chr_2 type=cds
>TCIR_00035130-00001 gene=TCIR_00035130 seq_id=tci2_wsi2.0_chr_2 type=cds
>TCIR_00076180-00001 gene=TCIR_00076180 seq_id=tci2_wsi2.0_chr_3 type=cds
>TCIR_00211160-00001 gene=TCIR_00211160 seq_id=tci2_wsi2.0_scaf_2039 type=cds

# checked each protein sequence, and there is only one problematic one (TCIR_10038850), so will remove it. Also removed scaffs with genes with a single duplicated gene present elsewhere in the genome
cat tcircumcincta.v2.5.gff3.tmp | 
    grep -v "TCIR_10038850" |\
    grep -v -w "tci2_wsi2.0_scaf_683" |\
    grep -v -w "tci2_wsi2.0_scaf_977" |\
    grep -v -w "tci2_wsi2.0_scaf_1891" |\
    grep -v -w "tci2_wsi2.0_scaf_1688" |\
    grep -v -w "tci2_wsi2.0_scaf_1612" |\
    grep -v -w "tci2_wsi2.0_scaf_1137" > tcircumcincta.v2.5.gff3






# generate protein sequences from gff
conda activate agat

agat_sp_extract_sequences.pl -g tcircumcincta.v2.5.gff3  -f ../teladorsagia_circumcincta_tci2_wsi2.4.fa -p --output tcircumcincta.v2.5.proteins.fa


# check protein completeness of original, and filtered for longest isoform

# -- extract longest isoforms from GFF, and then create proteins

agat_sp_keep_longest_isoform.pl -gff tcircumcincta.v2.5.gff3  -out tcircumcincta.v2.5.longest-isoform.gff3

agat_sp_extract_sequences.pl -g tcircumcincta.v2.5.longest-isoform.gff3  -f ../teladorsagia_circumcincta_tci2_wsi2.4.fa -p --output tcircumcincta.v2.5.longest-isoform.fa

conda deactivate

#- run busco
conda activate busco_5.6.1
busco --in tcircumcincta.v2.5.proteins.fa --out busco_tcircumcincta.v2.5.proteins --mode protein --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 -f -r

#> C:96.0%[S:48.1%,D:47.9%],F:0.9%,M:3.1%,n:3131

busco --in tcircumcincta.v2.5.longest-isoform.fa --out busco_tcircumcincta.v2.5.longest-isoform --mode protein --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 -f -r

#> C:95.9%[S:82.5%,D:13.4%],F:0.8%,M:3.3%,n:3131


cat tcircumcincta.v2.5.gff3 | sort -k1,1 -k4,4n | bgzip > tcircumcincta.v2.5.gff3.gz
tabix tcircumcincta.v2.5.gff3.gz


cat tcircumcincta.v2.5.gff3 | awk -F'[\t;]' '{if($3=="gene") print $9,$1,$4,$5,$7}' OFS="\t" | sed 's/ID=//g' | sort -k2,2V -k3,3n > tcircumcincta.v2.5.genes.list

```