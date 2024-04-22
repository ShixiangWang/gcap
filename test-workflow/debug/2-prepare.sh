# Set up conda env
# mamba create -n cancerit -c bioconda cancerit-allelecount 

cd /data3/wsx/share/gcap_reference
wget -c https://zenodo.org/records/6524005/files/1000G_loci_hg38.tar.gz
wget -c https://zenodo.org/records/6524005/files/GC_correction_hg38.txt.gz
wget -c https://zenodo.org/records/6524005/files/RT_correction_hg38.txt.gz

tar zxvf 1000G_loci_hg38.tar.gz
gunzip GC_correction_hg38.txt.gz
gunzip RT_correction_hg38.txt.gz