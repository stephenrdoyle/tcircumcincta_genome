# Teladorsagia circumcincta genome analysis: repeats

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk





```bash

cd /nfs/users/nfs_s/sd21/lustre_link/teladorsagia_circumcincta/GENOME/REPEATS

ln -s ../REFERENCES/DNAZOO/Tcircumcincta.DNAzoo.fa DNAZOO.fa

bsub.py --queue long --threads 20 50 repeat_masker_dnazoo "~sd21/bash_scripts/run_genome_repeatmasker DNAZOO DNAZOO.fa"


# WASHU
ln -s ../REFERENCES/WASHU/teladorsagia_circumcincta.PRJNA72569.WBPS17.genomic.fa WASHU.fa

bsub.py --queue long --threads 20 50 repeat_masker_washu "~sd21/bash_scripts/run_genome_repeatmasker WASHU WASHU.fa"

```