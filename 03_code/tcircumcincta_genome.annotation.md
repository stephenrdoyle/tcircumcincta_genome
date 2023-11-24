# Teladorsagia circumcincta genome annotation




Approaches
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






# merged RNAseq data

```
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

```
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


# braker 
C:89.1%[S:42.7%,D:46.4%],F:0.7%,M:10.2%,n:3131

# braker 2.4
C:88.7%[S:43.6%,D:45.1%],F:0.6%,M:10.7%,n:3131

# braker + augustus 2.4
C:95.3%[S:46.9%,D:48.4%],F:1.0%,M:3.7%,n:3131

# Haemonchus
C:96.2%[S:85.1%,D:11.1%],F:0.5%,M:3.3%,n:3131





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


# filter the data
conda activate agat

/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/AGAT/bin/agat_sp_filter_feature_from_keep_list.pl --gff braker_augustus.2.4.gt.gff3 --keep_list gene_ids_to_keep.list --output braker_augustus.2.4.minus.curated.gff3

# add start and stops
/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/AGAT/bin/agat_sp_add_start_and_stop.pl --gff braker_augustus.2.4.minus.curated.gff3 --fasta teladorsagia_circumcincta_tci2_wsi2.4.fa --out braker_augustus.2.4.minus.curated_start-stop-fix.gff3

# merge braker+augustus with apollo
gt gff3 -tidy -sort -force -o braker_augustus.apollo.2.4.gff3 braker_augustus.2.4.minus.curated_start-stop-fix.gff3 Annotations.tc2.4.gt.gff3

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






## interproscan

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


