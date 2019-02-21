#! /bin/sh
xargs -n 1 curl -O -L < links.txt

rm nohup.out

mv ENCFF000BSM.bam HUVEC_H3K4me1_rep1.bam
mv ENCFF000BSN.bam HUVEC_H3K4me1_rep2.bam
mv ENCFF000BSO.bam HUVEC_H3K4me1_rep3.bam
mv ENCFF353JMR.bam HUVEC_H3K4me3_rep1.bam
mv ENCFF456FFM.bam HUVEC_H3K4me3_rep2.bam
mv ENCFF977XML.bam HUVEC_H3K27ac_rep1.bam
mv ENCFF757GIT.bam HUVEC_H3K27ac_rep2.bam
mv ENCFF921ASC.bam HUVEC_H3K27me3_rep1.bam
mv ENCFF028HLI.bam HUVEC_H3K27me3_rep2.bam

samtools index HUVEC_H3K4me1_rep1.bam
samtools index HUVEC_H3K4me1_rep2.bam
samtools index HUVEC_H3K4me1_rep3.bam
samtools index HUVEC_H3K4me3_rep1.bam
samtools index HUVEC_H3K4me3_rep2.bam
samtools index HUVEC_H3K27ac_rep1.bam
samtools index HUVEC_H3K27ac_rep2.bam
samtools index HUVEC_H3K27me3_rep1.bam
samtools index HUVEC_H3K27me3_rep2.bam

samtools merge HUVEC_H3K4me1.bam HUVEC_H3K4me1_rep1.bam HUVEC_H3K4me1_rep2.bam HUVEC_H3K4me1_rep3.bam
samtools merge HUVEC_H3K4me3.bam HUVEC_H3K4me3_rep1.bam HUVEC_H3K4me3_rep2.bam
samtools merge HUVEC_H3K27ac.bam HUVEC_H3K27ac_rep1.bam HUVEC_H3K27ac_rep2.bam
samtools merge HUVEC_H3K27me3.bam HUVEC_H3K27me3_rep1.bam HUVEC_H3K27me3_rep2.bam

samtools index HUVEC_H3K4me1.bam
samtools index HUVEC_H3K4me3.bam
samtools index HUVEC_H3K27ac.bam
samtools index HUVEC_H3K27me3.bam

rm HUVEC_H3K4me1_rep1.bam HUVEC_H3K4me1_rep1.bam.bai HUVEC_H3K4me1_rep2.bam HUVEC_H3K4me1_rep2.bam.bai HUVEC_H3K4me1_rep3.bam HUVEC_H3K4me1_rep3.bam.bai HUVEC_H3K4me3_rep1.bam HUVEC_H3K4me3_rep1.bam.bai HUVEC_H3K4me3_rep2.bam HUVEC_H3K4me3_rep2.bam.bai HUVEC_H3K27ac_rep1.bam HUVEC_H3K27ac_rep1.bam.bai HUVEC_H3K27ac_rep2.bam HUVEC_H3K27ac_rep2.bam.bai HUVEC_H3K27me3_rep1.bam HUVEC_H3K27me3_rep1.bam.bai HUVEC_H3K27me3_rep2.bam HUVEC_H3K27me3_rep2.bam.bai

