#! /bin/sh
wget -i ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_HUVEC_combined_30.hic.gz
rm wget-log
gunzip GSE63525_HUVEC_combined_30.hic.gz
