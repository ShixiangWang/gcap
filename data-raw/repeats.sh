# https://www.biostars.org/p/278914/

wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.txt.gz | gunzip -c - | awk -v OFS="\t" '{ print $6, $7, $8, $12, $11, $10 }' - | sort -k1,1 -k2,2n > rmsk_hg38.bed

wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/rmsk.txt.gz | gunzip -c - | awk -v OFS="\t" '{ print $6, $7, $8, $12, $11, $10 }' - | sort -k 1,1 -k2,2n > rmsk_hg19.bed

