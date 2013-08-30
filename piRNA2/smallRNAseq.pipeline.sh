################################################################
# The wrapper script for the smallRNAseq analysis pipeline
# Usage:
# qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/fromPublication/Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA.fastq Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA 0 /home/dongx/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks.trinity.transcript.bed12
# qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_10dpp_wt_smallRNA_ox_rep1.fastq mouse_10dpp_wt_smallRNA_ox_rep1 0 /home/dongx/scratch/mouse_10.5dpp_wt_RNAseq_PE50nt_strand/cufflinks.trinity.transcript.bed12
# version: 1.0
# date: 2013-02-28
# Author: Xianjun Dong
# Limit: this only works for single-end smallRNA fastq input
################################################################

#!/bin/sh
#$ -V
#$ -pe openmpi 24
#$ -cwd
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=10G

export BOWTIE_INDEXES=$GENOME/mm9/Sequence/BowtieIndex/
export ANNOTATION=$GENOME/mm9/Annotation
#ncRNA=$ANNOTATION/ncRNA/ncRNA.wt_miRNA_piRNA.bed12  # ncRNA from fRNAdb (http://www.ncrna.org/frnadb/download)
ncRNA=$ANNOTATION/ncRNA/biomart.mm9.snRNA_snoRNA_tRNA_rRNA.bed
miRNA=$ANNOTATION/SmallRNA/mm10.liftover.to.mm9.bed  # precursor from miRBase
repeats=$ANNOTATION/Variation/repeatmasker.mm9.ucsc.sorted.bed  # repeat from repeatmasker
genes=$ANNOTATION/Genes/genes.gtf  # mRNA
exons=$ANNOTATION/Genes/genes.exons.bed12  # mRNA
cds=$ANNOTATION/Genes/genes.cds.bed12  # mRNA
introns=$ANNOTATION/Genes/genes.introns.bed12
three_prime_utr=$ANNOTATION/Genes/genes.3utr.bed12
five_prime_utr=$ANNOTATION/Genes/genes.5utr.bed12

cpu=24

###########################################
echo "############## 0. Configuration"
###########################################
# include
source $HOME/projects/mylib.sh
export PATH=$PATH:$HOME/projects/piRNA/src

fullfilepath=$1
samplename=$2
mm=$3
annotation=$4 # bed12 format
# gtf2bed /home/dongx/scratch/mouse_04dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf | cat - ~/nearline/Xin/RNAseq/mouse_04dpp_wt_RNAseq_PE100nt_strand/Phil.Rnaseq.nonPolyA.B6.4dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed | sort -k1,1 -k2,2n > /home/dongx/scratch/mouse_04dpp_wt_RNAseq_PE100nt_strand/cufflinks.trinity.transcript.bed12
# gtf2bed < /home/dongx/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf | cat - ~/nearline/Xin/RNAseq/mouse_07dpp_wt_RNAseq_PE100nt_strand/Phil.Rnaseq.nonPolyA.B6.7dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed | sort -k1,1 -k2,2n > /home/dongx/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks.trinity.transcript.bed12

adaptersequence="TCGTATGCCGTCTTCTGCTTG";

ext="${fullfilepath##*.}"
[[ "$samplename" == "" ]] && { samplename=$(basename "$fullfilepath"); samplename="${samplename%.*}";}
[[ "$mm" == "" ]] && mm=0;

# check the file format by file extension
[[ "$ext" =~ "fq|fastq|FASTQ" ]] || die "The input is not in fastq format. You have to convert it into fastq in order to run the pipeline.";

Min=16; Max=40

[ -d $HOME/scratch/ucsc ] || mkdir $HOME/scratch/ucsc
[ -d $HOME/scratch/$samplename ] || mkdir -p $HOME/scratch/$samplename

cd $HOME/scratch/$samplename
[ -e "adaptor.fa" ] || echo -e ">adapter1\n$adaptersequence" > adaptor.fa
ln -fs $fullfilepath $samplename.fq
ln -fs $HOME/sge_jobs_output/sge_job.$JOB_ID.out run.log

#############################################
##echo "############## 1. QC"
#############################################
#run_QC_smallRNA $samplename.fq
#
############################################
#echo "############## 1.2 Collapse reads into uniq species [optional for smallRNA only]"
############################################
#_insert2uniqreads filtered/$samplename.fq $Min $Max > $samplename.fa
#
############################################
#echo "############## 2. Mapping"
############################################
### map to genome
#[ "$mm" == "0" ] && bowtie -v 0 -a --un unmapped.fa -p 8 -f genome $samplename.fa 2>mapping.log | awk '{OFS="\t";print $3,$4,$4+length($5),$1,$7,$2,$4,$4+length($5),0,1,length($5),0;}' > map2genome.bed
#[ "$mm" != "0" ] && bowtie -v $mm -a --best --un unmapped.fa -p 8 -f genome $samplename.fa 2>mapping.log | awk '{OFS="\t";print $3,$4,$4+length($5),$1,$7,$2,$4,$4+length($5),0,1,length($5),0;}' > map2genome.bed

#echo "### map to transcriptome"
#mapping_junction_reads2 map2junctions -f $samplename.fa $mm $annotation
## merge reads mapped to genome and transcriptome (using sort-merge)
#sort -m -k4,4 <(sort -k4,4 map2junctions/junctions.bed) map2genome.bed > allmap.bed
#
## add mapping count (like NH tag)
#join -1 4 -2 1 -a 1 -a 2 allmap.bed <(groupBy -g 4 -c 4 -o count -i allmap.bed) | awk 'BEGIN{OFS="\t";}{print $2,$3,$4,$1,$13,$6,$7,$8,$9,$10,$11,$12}' > allmap.bed2
#LC_ALL=C sort --parallel=$cpu --buffer-size=2G -k1,1 -k2,2n allmap.bed2 > allmap.bed  # parallel version
## sort -k1,1 -k2,2n --parallel=$cpu allmap.bed2 > allmap.bed # single-processor version
#rm allmap.bed2
#
## unique mapper only
#awk '{if($5==1) print}' allmap.bed > uniqmap.bed
#
###########################################
echo "############## 3. discard reads mapped to miRNA(precursor) and tRNA/rRNA/snRNA/snoRNA, and also reads length not in (23,33)..."
###########################################
echo "### handle uniq-mappers"

# is known ncRNA? miRNA and other ncRNA (tRNA/rRNA/snRNA/snoRNA)
intersectBed -a uniqmap.bed -b $miRNA -s -wo -sorted -f 0.5> map2miRNA.uniqmap.bed
intersectBed -a uniqmap.bed -b $ncRNA -s -wo -sorted -f 0.5> map2ncRNA.uniqmap.bed
## Alternative1: remained ones with lenght in range of (23, 35)
#intersectBed -a allmap.bed -b <(cat $miRNA $ncRNA | sort -k1,1 -k2,2n -m) -s -wa -sorted -v -f 0.5| awk '{if(($3-$2)>23 && ($3-$2)<35) print}' > map2piRNA.bed2
## Alternative2: if a read maps to miRNA/ncRNA annotation, it won't go to the next step for calculating piRNA coveragex
rm map2piRNA.uniqmap.bed non24_40nt.uniqmap.bed
fgrep -f <(cat map2miRNA.uniqmap.bed map2ncRNA.uniqmap.bed | awk '{a[$4]=1;}END{for(i in a) print i;}') -vw uniqmap.bed | awk '{if(($3-$2)>=24 && ($3-$2)<=40) print >> "map2piRNA.uniqmap.bed"; else print >> "non24_40nt.uniqmap.bed";}'

## futher divide the putative piRNA reads to repeats/exons/UTR/intron/intergenic
# reads2repeats in individual class
rm map2repeat.uniqmap.*.bed
intersectBed -a map2piRNA.uniqmap.bed -b $repeats -s -wo -sorted -split -f 0.5| awk '{split($16, n, "__"); print >> "map2repeat.uniqmap."n[3]".bed";}'
intersectBed -a map2piRNA.uniqmap.bed -b $repeats -s -wa -sorted -split -v -f 0.5 > map2piRNA.uniqmap.nonrepeats.bed
# reads2gene
intersectBed -a map2piRNA.uniqmap.nonrepeats.bed -b $five_prime_utr -s -wa -split -u -f 0.5> map2five_prime_utr.uniqmap.bed
intersectBed -a map2piRNA.uniqmap.nonrepeats.bed -b $cds -s -wa -split -u -f 0.5> map2cds.uniqmap.bed
intersectBed -a map2piRNA.uniqmap.nonrepeats.bed -b $introns -s -wa -split -u -f 0.5> map2introns.uniqmap.bed
intersectBed -a map2piRNA.uniqmap.nonrepeats.bed -b $three_prime_utr -s -wa -split -u -f 0.5> map2three_prime_utr.uniqmap.bed
#read2intergenic
intersectBed -a map2piRNA.uniqmap.nonrepeats.bed -b $genes -s -v -f 0.5> map2intergenic.uniqmap.bed

[ -e mapping.statistics.normalized.txt ] && rm mapping.statistics.normalized.txt;
for i in uniqmap.bed map2miRNA.uniqmap.bed map2ncRNA.uniqmap.bed non24_40nt.uniqmap.bed map2repeat.uniqmap.LTR.bed map2repeat.uniqmap.LINE.bed map2repeat.uniqmap.SINE.bed map2repeat.uniqmap.DNA.bed map2repeat.uniqmap.Simple_repeat.bed map2repeat.uniqmap.Low_complexity.bed map2repeat.uniqmap.Satellite.bed map2repeat.uniqmap.RNA.bed map2repeat.uniqmap.rRNA.bed map2repeat.uniqmap.scRNA.bed  map2repeat.uniqmap.snRNA.bed map2repeat.uniqmap.srpRNA.bed map2repeat.uniqmap.tRNA.bed map2repeat.uniqmap.Other.bed  map2repeat.uniqmap.RC.bed map2repeat.uniqmap.Unknown.bed map2five_prime_utr.uniqmap.bed map2cds.uniqmap.bed map2introns.uniqmap.bed map2three_prime_utr.uniqmap.bed map2intergenic.uniqmap.bed; do type=`echo $i | sed 's/.bed//g'`; awk -v type=$type '{OFS="\t"; split($4,a,"_"); s0=s0+a[2]; s=s+(a[2]/$5);}END{print type, s0, s;}' $i >> mapping.statistics.normalized.txt; done

echo "### handle all-mappers (which takes time)"

# is known ncRNA? miRNA and other ncRNA (tRNA/rRNA/snRNA/snoRNA)
intersectBed -a allmap.bed -b $miRNA -s -wo -sorted -f 0.5> map2miRNA.allmap.bed
intersectBed -a allmap.bed -b $ncRNA -s -wo -sorted -f 0.5> map2ncRNA.allmap.bed
## remained ones with lenght in range of (23, 35)
#intersectBed -a allmap.bed -b <(cat $miRNA $ncRNA | sort -k1,1 -k2,2n -m) -s -wa -sorted -v -f 0.5| awk '{if(($3-$2)>23 && ($3-$2)<40) print}' > map2piRNA.bed2
## if a read maps to miRNA/ncRNA annotation, it won't go to the next step for calculating piRNA coveragex
rm map2piRNA.allmap.bed non24_40nt.allmap.bed
fgrep -f <(cat map2miRNA.allmap.bed map2ncRNA.allmap.bed | awk '{a[$4]=1;}END{for(i in a) print i;}') -vw allmap.bed | awk '{if(($3-$2)>=24 && ($3-$2)<=40) print >> "map2piRNA.allmap.bed"; else print >> "non24_40nt.allmap.bed";}'

## futher divide the putative piRNA reads to repeats/exons/UTR/intron/intergenic
# reads2repeats in individual class
rm map2repeat.allmap.*.bed
intersectBed -a map2piRNA.allmap.bed -b $repeats -s -wo -sorted -split -f 0.5| awk '{split($16, n, "__"); print >> "map2repeat.allmap."n[3]".bed";}'
intersectBed -a map2piRNA.allmap.bed -b $repeats -s -wa -sorted -split -v -f 0.5 > map2piRNA.allmap.nonrepeats.bed
# reads2exons
intersectBed -a map2piRNA.allmap.nonrepeats.bed -b $five_prime_utr -s -wa -split -u -f 0.5> map2five_prime_utr.allmap.bed
intersectBed -a map2piRNA.allmap.nonrepeats.bed -b $cds -s -wa -split -u -f 0.5> map2cds.allmap.bed
intersectBed -a map2piRNA.allmap.nonrepeats.bed -b $introns -s -wa -split -u -f 0.5> map2introns.allmap.bed
intersectBed -a map2piRNA.allmap.nonrepeats.bed -b $three_prime_utr -s -wa -split -u -f 0.5> map2three_prime_utr.allmap.bed
#read2intergenic
intersectBed -a map2piRNA.allmap.nonrepeats.bed -b $genes -s -v -f 0.5> map2intergenic.allmap.bed

## get stat
for i in allmap.bed map2miRNA.allmap.bed map2ncRNA.allmap.bed non24_40nt.allmap.bed map2repeat.allmap.LTR.bed map2repeat.allmap.LINE.bed map2repeat.allmap.SINE.bed map2repeat.allmap.DNA.bed map2repeat.allmap.Simple_repeat.bed map2repeat.allmap.Low_complexity.bed map2repeat.allmap.Satellite.bed map2repeat.allmap.RNA.bed map2repeat.allmap.rRNA.bed map2repeat.allmap.scRNA.bed  map2repeat.allmap.snRNA.bed map2repeat.allmap.srpRNA.bed map2repeat.allmap.tRNA.bed map2repeat.allmap.Other.bed  map2repeat.allmap.RC.bed map2repeat.allmap.Unknown.bed map2five_prime_utr.allmap.bed map2cds.allmap.bed map2introns.allmap.bed map2three_prime_utr.allmap.bed map2intergenic.allmap.bed; do type=`echo $i | sed 's/.bed//g'`; awk -v type=$type '{OFS="\t"; split($4,a,"_"); s0=s0+a[2]; s=s+(a[2]/$5);}END{print type, s0, s;}' $i >> mapping.statistics.normalized.txt; done


echo "Job done!"; exit

echo "
df=read.table('mapping.statistics.normalized.txt', header=F)
rownames(df)=df[,1]; df=df[,-1]; colnames(df)=c('rawCount','normCount')
uniqmap=df[grep('.uniqmap', rownames(df)),2]; names(uniqmap)=sub('.uniqmap', '', rownames(df)[grep('.uniqmap', rownames(df))]); names(uniqmap)=sub('map2','',names(uniqmap))
allmap=df[grep('.allmap', rownames(df)),2]; names(allmap)=sub('.allmap', '', rownames(df)[grep('.allmap', rownames(df))]); names(allmap)=sub('map2','',names(allmap))

pdf('$samplename.mapping.statistics.normalized.pdf')
par(mfrow=c(2,1))
#pie(uniqmap)
par(mar=c(8,5,2,2)); barplot(uniqmap, las=2, cex.name=0.8, main='unique mapping') #, log='y')
#pie(allmap)
par(mar=c(9,5,2,2)); barplot(allmap, las=2, cex.name=0.8, main='all mapping')
dev.off()
" > /tmp/stat.R
Rscript /tmp/stat.R
