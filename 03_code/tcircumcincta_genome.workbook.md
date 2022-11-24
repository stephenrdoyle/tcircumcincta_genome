
# purge duplicates

cp /nfs/users/nfs_s/sd21/lustre118_link/tc//PURGEDUPS/pb.fofn .

cat pb.fofn
/lustre/scratch118/infgen/team133/sd21/tc/PACBIO_CANU_NEW/all_subreads.fastq


pd_config.py tc_canu1.9.contigs.fasta pb.fofn


# run purge dups - note, does not need to be bsubbed, as it automatically submits jobs
run_purge_dups.py config.json ~sd21/lustre118_link/software/GENOME_IMPROVEMENT/purge_dups/bin TC


# wrapper
!#/bin/bash

# run_purge_dups

PREFIX=$1
REFERENCE=$2
PB_FOFNs=$3

pd_config.py ${REFERENCE} ${PB_FOFNs}

run_purge_dups.py config.json ~sd21/lustre118_link/software/GENOME_IMPROVEMENT/purge_dups/bin ${PREFIX}
