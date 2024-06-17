# Cleaning up the genome

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk


- filtering the extra contigs not scaffolded in the genome
- identfying true contamination, ie, bacterial contigs
    - blobtools
- junk sequences not containing anything of obvious use to the assembly
    - screen against existing proteomes to id contigs without putative coding sequences


## Blobtools 
```bash
# working directory 
cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/BLOBTOOLS
```


```bash
# get taxdump
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -

#--- fix reference, and split it into 100 reads per file - this will help speed up the blast
ln -s ../FIX_MTDNA/CURRENT_ASSEMBLY_mtDNA.fa REF.fa


fastaq split_by_base_count REF.fa REFsplit 10000000



#--- run blast (vBLASTN 2.12.0+), with a separate blast job per 100 reads
for i in REFsplit*; do \
     bsub.py 10 --queue long --threads 16 blast blastn \
     -db /data/blastdb/Supported/NT/current/nt \
     -query ${i} \
     -outfmt \'6 qseqid staxids bitscore std sscinames sskingdoms stitle\' \
     -max_target_seqs 10 \
     -max_hsps 1 \
     -evalue 1e-25 \
     -num_threads 16 \
     -out blast.out_${i};
done

cat blast.out* > total_blast.out
rm x*

# diamond - didnt end up using it
# module load diamond/2.0.12 

# #------ run once to initially set up diamond
# cd /nfs/users/nfs_s/sd21/lustre_link/databases/diamond

# bsub.py 10 diamond_makedb_v2 "diamond makedb --in uniprot_ref_proteomes.fasta --taxonmap prot.accession2taxid --db uniprot_ref_proteomes"

# #------

# for i in REFsplit*; do \
# bsub.py 50 --queue small --threads 2 diamond "diamond blastx \
#  --query ${i} \
#  --db /nfs/users/nfs_s/sd21/lustre_link/databases/blobtools/uniprot/reference_proteomes.dmnd \
#  --outfmt 6 \
#  --sensitive \
#  --max-target-seqs 1 \
#  --evalue 1e-25" \
#  --salltitles \
#  --sallseqid \
#  --out diamond.out_${i};
#  done

# diamond.out

# cat diamond.out* > total_diamond.out
# rm x*


# generate coverage data
#--- use minimap to map reads to generate a bam for blobtools coverage stats

ln -s /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/DATA/illumina/42782_8_1_R1.fastq.gz
ln -s /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/DATA/illumina/42782_8_1_R2.fastq.gz

bsub.py --queue long --threads 16 20 mapping "minimap2 -ax sr -t 16 REF.fa 42782_8_1_R1.fastq.gz 42782_8_1_R2.fastq.gz \| samtools sort -@16 -O BAM -o mapped.bam - "

samtools index mapped.bam

# run blobtools
conda activate blobtools

/nfs/users/nfs_s/sd21/lustre_link/software/ASSEMBLY_QC/blobtools/blobtools taxify --hit_file total_diamond.out --taxid_mapping_file /nfs/users/nfs_s/sd21/lustre_link/databases/blobtools/uniprot/reference_proteomes.taxid_map --map_col_sseqid 0 --map_col_taxid 2

/nfs/users/nfs_s/sd21/lustre_link/software/ASSEMBLY_QC/blobtools/blobtools create -i REF.fa -b mapped.bam  --hitsfile total_diamond.out -o blobtools_out --nodes /nfs/users/nfs_s/sd21/lustre_link/databases/blobtools/uniprot/nodes.dmp --names /nfs/users/nfs_s/sd21/lustre_link/databases/blobtools/uniprot/names.dmp

/nfs/users/nfs_s/sd21/lustre_link/software/ASSEMBLY_QC/blobtools/blobtools plot -i blobtools_out.blobDB.json

/nfs/users/nfs_s/sd21/lustre_link/software/ASSEMBLY_QC/blobtools/blobtools view -i blobtools_out.blobDB.json

```

- something wrong with the nt database - dont seem to be able to get the taxonomic identifiers, which then prevents blobtools from colouring its plots.
- this is ok - gives some sense of the contamination present from the blast data
    - can use this to manually filter
- small amount of variance in coverage and GC content - need to check





## manual check of Blast data
- to used the blast data to screen for contaminants, I made a "contam.list" and a "keep.list" based on wether then had a hit to a non-worm or worm reference in the blast database
    - these lists simply had accession IDs from the blast hits, sorted and unique.
    - manually checked this list in excel

```bash
# extracting accessions and hit descripts from the blast output, sorting, and unique
cut -f5,17 total_blast.out | sort | uniq | head

# AB549207.1	Nippostrongylus brasiliensis msp1 mRNA for major sperm protein 1, partial cds
# AC087415.1	Caenorhabditis briggsae cosmid CB005O21, complete sequence
# AF000967.1	Haemonchus contortus putative glutamate dehydrogenase (HCGLDH1) mRNA, complete cds
# AF052044.1	Ostertagia ostertagi ribosomal protein L38 (rpl38) mRNA, complete cds
# AF052046.1	Ostertagia ostertagi clone jmo21 larval stage L3 mRNA sequence
# AF052591.1	Ostertagia ostertagi JMO38 larval stage L3 unknown gene
# AF079402.1	Haemonchus contortus pepsinogen (Pep1) gene, complete cds
# AF099908.1	Haemonchus contortus transposase homolog (tc1A) gene, complete cds
# AF105337.1	Ostertagia circumcincta galectin GAL-1 gene, complete cds
# AJ309529.1	Ostertagia ostertagi partial mRNA for aspartic protease (l4 asp1 gene)


wc -l contam.list
#> 175 contam.list

wc -l keep.list
# 239 keep.list

```

- used these IDs to count the number of hits per contig/scaffold
    - sequences with all comtam hits to be removed
    - sequences with contam and keep hits to be checked
    - sequences with only keep hits to be retained
- also gathered BUSCO counts per sequence, just to cross check

```bash
cat REF.fa.fai | cut -f1 | while read SEQ; do 

CONTAM_N=$(cat total_blast.out | grep -w $SEQ | grep -c -f contam.list)
KEEP_N=$(cat total_blast.out | grep -w $SEQ | grep -c -f keep.list)
BUSCO_N=$(cat ../RECOVER_MISSING_BUSCOs/Tc_curated.busco_recovered.reduced.genome_nematoda_odb10/run_nematoda_odb10/full_table.tsv | grep -c -w $SEQ )

echo -e "${SEQ}\t${KEEP_N}\t$CONTAM_N\t${BUSCO_N}";
done

```

- sorting this list identifed 30 sequences with only contamination hits to be removed.
scaffold_2842
scaffold_259
scaffold_220
scaffold_3465
scaffold_2607
scaffold_1804
scaffold_3480
scaffold_2134
scaffold_2480
scaffold_1520
scaffold_3780
scaffold_780
scaffold_2481
scaffold_408
scaffold_250
scaffold_3044
scaffold_1623
scaffold_3772
scaffold_3672
scaffold_1592
scaffold_515
scaffold_1066
scaffold_2630
scaffold_574
scaffold_591
scaffold_3313
scaffold_815
scaffold_2821
scaffold_1504
scaffold_3235

- this approach also identified 966 sequences with no blast hits at all. To see if these were junk to be removed, or contained useful hits, I downloaded the CDS sequences from Haemonchus and Teladorasagia from WBP, and used promer to map to the extra sequences. 

```bash
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS17.CDS_transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/teladorsagia_circumcincta/PRJNA72569/teladorsagia_circumcincta.PRJNA72569.WBPS17.CDS_transcripts.fa.gz

gunzip haemonchus_contortus.PRJEB506.WBPS17.CDS_transcripts.fa.gz
gunzip teladorsagia_circumcincta.PRJNA72569.WBPS17.CDS_transcripts.fa.gz

# cleanup the data
cat haemonchus_contortus.PRJEB506.WBPS17.CDS_transcripts.fa teladorsagia_circumcincta.PRJNA72569.WBPS17.CDS_transcripts.fa > hcon_tcirc_cds.fa
cut -f1 -d " " hcon_tcirc_cds.fa > tmp; mv tmp hcon_tcirc_cds.fa
fastaq enumerate_names --suffix hcon_tcirc_cds_ hcon_tcirc_cds.fa hcon_tcirc_cds.renamed.fa

# run promer
bsub.py 10 promer "promer extra.list.fa hcon_tcirc_cds.renamed.fa"


show-coords -lTH out.delta | cut -f14 | sort | uniq | wc -l
# 596 sequences with a promer hit to haem or tcirc - will be conservative and keep these

# combined the 596 sequences with the total sequence count ("junk_screen.txt"), and selected only those found once - these are the extra sequences to remove
sort junk_screen.txt | uniq -c | awk '$1==1 {print $2}' | wc -l
# 370 sequences to remove

```


## Final clean up
- removing the 30 contaminating sequences and 370 junk sequences from the assembly
- 

```bash
# extract everything except the chromosomes
remove_ids=($(awk '{print $1}' CURRENT_ASSEMBLY_mtDNA.fa.fai | grep -x -v -f remove.list))

samtools faidx -o ASSEMBLY_230308.fa CURRENT_ASSEMBLY_mtDNA.fa "${remove_ids[@]}"


```







```bash
bsub.py --queue long --threads 20 20 prothint "/nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/ProtHint/bin/prothint.py --threads 20 extra.list.fa /nfs/users/nfs_s/sd21/lustre_link/software/TRANSCRIPTOME/ProtHint/databases/odb10_metazoa/odb10_metazoa_proteins.fa"

bsub.py --queue long --threads 20 20 braker2 "/lustre/scratch118/infgen/team133/sd21/software/TRANSCRIPTOME/BRAKER_v2.0/braker.pl --genome=cjohnstoni_genome_200917.fa --hints=prothint_augustus.gff3 --gff3 --cores 20 --extrinsicCfgFile=extrinsic.M.RM.E.W.P.C.cfg"

```




```bash
cut -f3 full_table.tsv | sort | uniq 


# post polishing contam screen
/nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/ASSEMBLY/POST_CANU_IMPROVEMENT/POLISH/pb_arrow_polish_tc_ASSEMBLY_230308_out/POLCA

# print only contigs with duplicate BUSCOS, but no complete or fragmented BUSCOs
cut -f3 BUSCO_nematoda/run_nematoda_odb10/full_table.tsv | sort | uniq | grep -v "#" | while read -r NAME; do 
    Complete=$(cat BUSCO_nematoda/run_nematoda_odb10/full_table.tsv | grep $NAME | grep -c "Complete" ); 
    Duplicated=$(cat BUSCO_nematoda/run_nematoda_odb10/full_table.tsv | grep $NAME |grep -c "Duplicated" ); 
    Fragmented=$(cat BUSCO_nematoda/run_nematoda_odb10/full_table.tsv | grep $NAME | grep -c "Fragmented" ); 
    echo -e "${NAME}\t${Complete}\t${Duplicated}\t${Fragmented}"; 
    done | awk '{if($2==0 && $4==0) print $1}'  > contigs.w.duplicates.list


while read -r NAME; do 
    samtools faidx tc_ASSEMBLY_230308.arrow.fa.PolcaCorrected.fa ${NAME}; 
    done < contigs.w.duplicates.list > contigs.w.duplicates.fa


nucmer -p contigs.w.duplicates_v_TcHc contigs.w.duplicates.fa ./../../../PRE_POLISH_COMPLETE/hcon_tcirc_cds.renamed.fa




# extract everything except the chromosomes
ln -s contigs.w.duplicates.list remove.list

remove_ids=($(awk '{print $1}' tc_ASSEMBLY_230308.arrow.fa.PolcaCorrected.fa.fai | grep -x -v -f remove.list))

samtools faidx -o ASSEMBLY_230317.fa tc_ASSEMBLY_230308.arrow.fa.PolcaCorrected.fa "${remove_ids[@]}"
```