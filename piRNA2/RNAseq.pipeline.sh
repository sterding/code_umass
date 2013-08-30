################################################################
# The wrapper script for the RNAseq analysis pipeline
# Usage: TOADD
# version: 1.0
# date: 2013-02-28
# Author: Xianjun Dong
################################################################


#!/bin/sh
#$ -V
#$ -pe single 8
#$ -cwd
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=4G

# include
source $HOME/projects/mylib.sh

###########################################
echo "############## 0. Configuration"
###########################################

samplename=$1 #mouse_adult_wt_RNAseq_PE50nt_strand_R2
read1=$2 # R1.fastq
read2=$3 # R2.fastq
smallRNAreads=$4 #$HOME/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.ox.6wk.testes.raw.uniq.reads.fa.gz

trinity_output=$HOME/nearline/Xin/RNAseq/$samplename/trinity.introns.bed
#link to /home/wangj2/wnearline/xin_rnaseq/junc_pirna_cluster/Phil.Rnaseq.mouse_adult_testis.npa.trinity.blat.maplen95p.mm1p.allmap.introns.bed # intron in bed6

cpu=8

export BOWTIE2_INDEXES=$GENOME/mm9/Sequence/Bowtie2Index/
export ANNOTATION=$GENOME/mm9/Annotation/Genes
[ -d $HOME/scratch/ucsc ] || mkdir $HOME/scratch/ucsc
[ -d $HOME/scratch/$samplename ] || mkdir -p $HOME/scratch/$samplename

ln -fs $HOME/sge_jobs_output/sge_job.$JOB_ID.out $HOME/scratch/$samplename/sge2.log

cd $HOME/nearline/Xin/RNAseq/$samplename

phred=`getphred $read1`; tophat_scoreoption=""; # default
[ "$phred" == "Phred+64" ] && tophat_scoreoption="--solexa1.3-quals";
[ "$phred" == "Solexa+64" ] && tophat_scoreoption="--solexa-quals";
echo "tophat_scoreoption: $tophat_scoreoption"

###########################################
echo "############## 1. QC"
###########################################
run_QC R1.fastq R2.fastq

###########################################
echo "############## 2. Mapping"
###########################################
run_mapping R1.fastq R2.fastq

# mapping stat
get_mapping_stat mapping.sam

###########################################
echo "############## 3. Assembly"
###########################################
run_assembly mapping.bam

###########################################
echo "############## 3.2 Refine with CAGE/PAS"
###########################################
run_5refinement annotation.gtf CAGE_narrowpeak
run_3refinement annotation.gtf PAS_narrowpeak

###########################################
echo "############## 4. quantify mRNA"
###########################################
run_quanfiy_mRNA annotation.gtf mapping.bam

###########################################
echo "############## 5. quantify smallRNA"
###########################################
run_quanfiy_smallRNA transcripts.gtf $HOME/scratch/Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA/map2piRNA.bed
run_quanfiy_smallRNA transcripts.gtf $HOME/scratch/Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA/map2piRNA_u.bed
