# Teladorsagia circumcincta genome analysis: polishing the genome

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk



- polishing the genome in two ways, 
    - first with pacbio reads
    - second with pilon



```bash 
# working dir
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/POLISH

ln -s ~/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/PRE_POLISH_COMPLETE/ASSEMBLY_230308.fa .

./run_arrow_pacbio_genome_polish.sh

```

- where "run_arrow_pacbio_genome_polish.sh" is a wrapper script to speed up the arrow steps by mapping and them splitting to polish each sequence in an array:
```bash
#!/bin/bash

#-------------------------------------------------------------------------------
# run_arrow_pacbio_genome_polish.sh
#-------------------------------------------------------------------------------

# stephen doyle
# Jan 2023


export PREFIX=tc_ASSEMBLY_230308
export REFERENCE=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/PRE_POLISH_COMPLETE/ASSEMBLY_230308.fa
export PB_READ_DATA_FOFN=/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/POLISH/subread_bams.fofn

# load required modules
module load pacbio-smrttools/7.0.1.66768
module load samtools/1.14--hb421002_0

## file locations
export LOG_FILES="$PWD/pb_arrow_polish_${PREFIX}_out/LOG_FILES"
export REFERENCE_FILES="$PWD/pb_arrow_polish_${PREFIX}_out/REFERENCE_FILES"
export MAPPING="$PWD/pb_arrow_polish_${PREFIX}_out/MAPPING"
export SPLIT_FASTAS="$PWD/pb_arrow_polish_${PREFIX}_out/SPLIT_FASTAS"
export SPLIT_BAMS="$PWD/pb_arrow_polish_${PREFIX}_out/SPLIT_BAMS"
export SPLIT_POLISHED_FASTAS="$PWD/pb_arrow_polish_${PREFIX}_out/SPLIT_POLISHED_FASTAS"


### Check output directories exist & create them as needed
[ -d ${LOG_FILES} ] || mkdir -p ${LOG_FILES}
[ -d ${REFERENCE_FILES} ] || mkdir -p ${REFERENCE_FILES}
[ -d ${MAPPING} ] || mkdir -p ${MAPPING}
[ -d ${SPLIT_FASTAS} ] || mkdir -p ${SPLIT_FASTAS}
[ -d ${SPLIT_BAMS} ] || mkdir -p ${SPLIT_BAMS}
[ -d ${SPLIT_POLISHED_FASTAS} ] || mkdir -p ${SPLIT_POLISHED_FASTAS}



# save current script in run folder to reproduce the exact output
cp ${PWD}/run_arrow_pacbio_genome_polish.sh ${PWD}/pb_arrow_polish_${PREFIX}_out/commands.$(date -Iminutes).txt

# make a progress file
#> ${PWD}/pb_arrow_polish_${PREFIX}.progress.log





#-------------------------------------------------------------------------------
### 01. Prepare reference files
#-------------------------------------------------------------------------------
# pbmm2: https://github.com/PacificBiosciences/pbmm2/

func_build_reference() {

if [ -f "${REFERENCE_FILES}/REF.fa" ]; then
        echo "Reference is already setup. Moving on."
        exit 0
    else

    cp ${REFERENCE} ${REFERENCE_FILES}/REF.fa

    pbmm2 index ${REFERENCE_FILES}/REF.fa ${REFERENCE_FILES}/REF.mmi

    samtools faidx ${REFERENCE_FILES}/REF.fa > ${REFERENCE_FILES}/REF.fa.fai

    cat ${REFERENCE_FILES}/REF.fa.fai | cut -f1 > ${REFERENCE_FILES}/sequences.list

    while read SEQUENCE; do
       echo "samtools faidx ${REFERENCE_FILES}/REF.fa ${SEQUENCE} > ${SPLIT_FASTAS}/${SEQUENCE}.fa";
    done < ${REFERENCE_FILES}/sequences.list | parallel -j20
fi
}

export -f func_build_reference



#-------------------------------------------------------------------------------
### 02. Map pacbio reads to reference
#-------------------------------------------------------------------------------
func_map_reads () {
    # check if bam file exits - if yes, then exit. Else, run mapping
    if [ -s "${MAPPING}/REF.bam" ]; then
            echo "Bam file is already setup. Moving on."
            exit 0
        else

        pbmm2 align ${REFERENCE_FILES}/REF.fa ${PB_READ_DATA_FOFN} ${MAPPING}/REF.bam --preset SUBREAD --sort -j 20 -J 20

        pbindex ${MAPPING}/REF.bam

        samtools index ${MAPPING}/REF.bam
    fi
}

export -f func_map_reads



#-------------------------------------------------------------------------------
### 03. Split bams per reference
#-------------------------------------------------------------------------------
func_split_bams () {

while read SEQUENCE; do
    echo "samtools view -b --reference ${SPLIT_FASTAS}/${SEQUENCE}.fa -o ${SPLIT_BAMS}/${SEQUENCE}.bam ${MAPPING}/REF.bam ${SEQUENCE}";
    done < ${REFERENCE_FILES}/sequences.list  | parallel -j20

while read SEQUENCE; do
    echo "pbindex ${SPLIT_BAMS}/${SEQUENCE}.bam";
    done < ${REFERENCE_FILES}/sequences.list | parallel -j20

while read SEQUENCE; do
    echo "samtools index ${SPLIT_BAMS}/${SEQUENCE}.bam";
    done < ${REFERENCE_FILES}/sequences.list | parallel -j20

}

export -f func_split_bams


#-------------------------------------------------------------------------------
### 04. Polish split bams and fasta sequences using arrow
#-------------------------------------------------------------------------------
func_arrow_polish () {

while read SEQUENCE; do
    echo "arrow --algorithm arrow --threaded -j 20 --referenceFilename ${SPLIT_FASTAS}/${SEQUENCE}.fa -o ${SPLIT_POLISHED_FASTAS}/${SEQUENCE}.arrow.fa ${SPLIT_BAMS}/${SEQUENCE}.bam";
    done < ${REFERENCE_FILES}/sequences.list | parallel -j20

ls -v ${SPLIT_POLISHED_FASTAS}/*.arrow.fa | xargs cat > $PWD/pb_arrow_polish_${PREFIX}_out/${PREFIX}.arrow.fa

}

export -f func_arrow_polish



#-------------------------------------------------------------------------------
# running the pipeline
#-------------------------------------------------------------------------------

# func_build_reference
bsub -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>7000] rusage[mem=7000]" -M7000 -o ${LOG_FILES}/01_build_reference.o -e ${LOG_FILES}/01_build_reference.e -J 01_build_reference_${PREFIX} func_build_reference

# func_map_reads
bsub -w "01_build_reference_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/02_map_reads.o -e ${LOG_FILES}/02_map_reads.e -J 02_map_reads_${PREFIX} func_map_reads

# func_split_bams
bsub -w "02_map_reads_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>10000] rusage[mem=10000]" -q long -M10000 -o ${LOG_FILES}/03_split_bams.o -e ${LOG_FILES}/03_split_bams.e -J 03_split_bams_${PREFIX} func_split_bams

# func_arrow_polish
bsub -w "03_split_bams_${PREFIX}" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>200000] rusage[mem=200000]" -q hugemem -n 20 -M200000 -o ${LOG_FILES}/04_arrow_polish.o -e ${LOG_FILES}/04_arrow_polish.e -J 04_arrow_polish_${PREFIX} func_arrow_polish
```













## Polishing with RNAseq reads
- want to determine if mapped RNAseq reads give an indication if there are unpolished indels in the genome
- will map RNAseq reads to 
    - unpolished
    - arrow polished
    - arrow + polca polished
    - arrow + polca + clean polished
- will determine 
    - mapping rates
    - presence and count of indels in CIGAR strings of mapped reads


```bash
# working dir
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/POLISH/RNAseq_test


# get references
ln -s /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/PRE_POLISH_COMPLETE/ASSEMBLY_230308.fa pre_polish.fa
ln -s ../pb_arrow_polish_tc_ASSEMBLY_230308_out/tc_ASSEMBLY_230308.arrow.fa arrow_polished.fa
ln -s ../pb_arrow_polish_tc_ASSEMBLY_230308_out/POLCA/tc_ASSEMBLY_230308.arrow.fa.PolcaCorrected.fa arrow_polca_polished.fa
ln -s ../pb_arrow_polish_tc_ASSEMBLY_230308_out/POLCA/ASSEMBLY_230317.fa arrow_polca_polished_plus_clean.fa

# get data - old RNAseq data from adult male, female and L3
zcat ../../../../TRANSCRIPTOME/DATA_ILLUMINA/*_1.fastq.gz | gzip > rnaseq_merged_1.fastq.gz
zcat ../../../../TRANSCRIPTOME/DATA_ILLUMINA/*_2.fastq.gz | gzip > rnaseq_merged_2.fastq.gz

```

# make index for each reference and map with STAR
```bash
for i in *fa; do
mkdir ${i%.fa}.index
bsub.py --threads 8 20 01_star_index \
"~sd21/lustre_link/software/TRANSCRIPTOME/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
--runMode genomeGenerate \
--genomeSAindexNbases 13 \
--runThreadN 8 \
--genomeDir ${i%.fa}.index \
--genomeFastaFiles ${i}";
done


for i in *fa; do \
bsub.py 10 --queue yesterday --threads 7 starmap_${i} ~sd21/lustre_link/software/TRANSCRIPTOME/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
--runThreadN 7 \
--genomeDir ${i%.fa}.index \
--readFilesIn rnaseq_merged_10m_1.fastq.gz rnaseq_merged_10m_2.fastq.gz \
--readFilesCommand zcat \
--alignIntronMin 10 \
--outTmpDir ${i%.fa}.star_out \
--outFileNamePrefix ${i%.fa}.star_out \
--limitBAMsortRAM 84327490898 \
--outSAMtype BAM SortedByCoordinate \
; done
```






# used RNAseq data to remove scaffolds
- remove
    - scaffolds with no RNAseq data
    - no Haem or Tcirc CDS hits
    - no BUSCOs
```bash
cat arrow_polca_polished_plus_clean.star_outAligned.sortedByCoord.out.chr.cov | sort -k3n | awk '$5==0 {print $0}' | cut -f1 | sed 's/|arrow//g' | while read NAME; do 
    name=$(echo $NAME); count=$(show-coords -lTH ../../PRE_POLISH_COMPLETE/out.delta | grep -c "$NAME"); echo -e "$name\t$count"; 
    done | awk '$2=="0" {print $0}' | while read NAME; do 
    name=$(echo $NAME); count=$(grep -c "$NAME" ../../RECOVER_MISSING_BUSCOs/Tc_curated.busco_recovered_genome_nematoda_odb10/run_nematoda_odb10/full_table.tsv); echo -e "$name\t$count"; 
    done | cut -f1 -d " " > remove.list

awk '{print $0"\|arrow"}' remove.list > tmp; mv tmp remove.list

remove_ids=($(awk '{print $1}' arrow_polca_polished_plus_clean.fa.fai | grep -x -v -f remove.list))

samtools faidx -o ASSEMBLY_230322.fa arrow_polca_polished_plus_clean.fa "${remove_ids[@]}"

```







### Polish with Polca using RNAseq reads
```bash
bsub.py --queue basement --threads 7 20 polca_RNAseq "/nfs/users/nfs_s/sd21/lustre_link/software/GENOME_IMPROVEMENT/MaSuRCA-4.1.0/global-1/PacBio/src_reconcile/polca.sh -a tc_ASSEMBLY_230308.arrow.fa.PolcaCorrected.fa -r \"rnaseq_merged_1.fastq.gz rnaseq_merged_2.fastq.gz\" -t 7 -m 2G"


bsub.py --queue long --threads 20 60 busco_tc_polca_RNAseq_nematoda_odb10 \
    "busco --in tc_ASSEMBLY_230308.arrow.fa.PolcaCorrected.fa.PolcaCorrected.fa --out busco_tc_polca_RNAseq_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"
```


```bash
# print only contigs with duplicate BUSCOS, but no complete or fragmented BUSCOs
cat full_table.tsv | grep -v "#" |  cut -f3  | sort | uniq  | while read -r NAME; do 
    Complete=$(cat full_table.tsv | grep $NAME | grep -c "Complete" ); 
    Duplicated=$(cat full_table.tsv | grep $NAME | grep -c "Duplicated" ); 
    Fragmented=$(cat full_table.tsv | grep $NAME | grep -c "Fragmented" ); 
    echo -e "${NAME}\t${Complete}\t${Duplicated}\t${Fragmented}"; 
    done | awk '{if($2==0 && $4==0) print $1}'  > contigs.w.duplicates.list


remove_ids=($(awk '{print $1}' tc_ASSEMBLY_230322.arrow.fa.PolcaCorrected.fa.PolcaCorrected.fa.fai | grep -x -v -f remove.list))

samtools faidx -o ASSEMBLY_230402.fa tc_ASSEMBLY_230322.arrow.fa.PolcaCorrected.fa.PolcaCorrected.fa "${remove_ids[@]}"


bsub.py --queue normal --threads 20 60 busco_tc_ASSEMBLY_230402_nematoda_odb10 \
    "busco --in ASSEMBLY_230402.fa --out busco_tc_ASSEMBLY_230402_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"

```