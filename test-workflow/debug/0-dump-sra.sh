#!/bin/bash
mkdir -p /data3/wsx/share/gcap_debug
cd /data3/wsx/share

export PATH=$HOME/soft/sratoolkit/bin:$PATH

for i in ERR5242993 ERR5243012
do
  echo handling $i
  parallel-fastq-dump -t 20 -O gcap_debug/ --split-3  --gzip -s $i
done
