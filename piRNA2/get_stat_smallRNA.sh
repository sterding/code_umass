#!/bin/sh
#$ -V
#$ -pe openmpi 8
#$ -cwd
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=2G

samplename=$1
output=$2

cd $HOME/scratch/$samplename;
echo $samplename > $output
echo "totalreads:" `wc -l $samplename.fq | awk '{print $1/4}'` >> $output
echo "totalreads_after_filtered:" `fgrep '>' $samplename.fa | cut -f2 -d'_' | addCols stdin` >> $output
cat allmap.bed.mapping.statistics.txt >> $output
cat uniqmap.bed.mapping.statistics.txt >> $output


