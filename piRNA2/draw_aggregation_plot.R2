# script to draw aggregation plot of CAGE/Degredome/PAS signal for meta-gene
# input files: genes annotation file
# Usage: Rscript draw_aggregation_plot.R control2
# Rscript draw_aggregation_plot.R piRNA
# requirement: run bigWigAverageOverBed_81bins.sh before (see detail there)

source('http://zlab.umassmed.edu/~dongx/mylib/mylib.R')

args <- commandArgs(TRUE)
genetype=args[1]
#source('http://zlab.umassmed.edu/~dongx/mylib/mylib.R'); genetype='control2'; type='EncodeCshlLongRnaSeq'; species='mouse'

species='mouse';

# gene file
if(genetype=='control2')
{
             ## use controlset
             g1=read.table("../data/controlSet2/piRNA.genic.bed"); rownames(g1)=paste(g1[,4], g1[,5], sep="__"); colnames(g1)=c('chr','start','end','name','score','strand','txstart','txend','rgb','exoncount','exonsize','exonstart')
             g2=read.table("../data/controlSet2/piRNA.intergenic.bed"); rownames(g2)=paste(g2[,4], g2[,5], sep="__"); colnames(g2)=c('chr','start','end','name','score','strand','txstart','txend','rgb','exoncount','exonsize','exonstart')
             g3=read.table("../data/controlSet2/x100910.rnaseq.transcripts.NM.all.final.bed.c"); rownames(g3)=paste(g3[,4], g3[,5], sep="__"); colnames(g3)=c('chr','start','end','name','score','strand','txstart','txend','rgb','exoncount','exonsize','exonstart')
             g4=read.table("../data/controlSet2/x100910.rnaseq.transcripts.NR.all.final.bed.c"); rownames(g4)=paste(g4[,4], g4[,5], sep="__"); colnames(g4)=c('chr','start','end','name','score','strand','txstart','txend','rgb','exoncount','exonsize','exonstart')
             #g5=read.table("../data/controlSet2/x100910.rnaseq.transcripts.NR.mir.bed"); rownames(g5)=paste(g5[,4], g5[,5], sep="__"); colnames(g5)=c('chr','start','end','name','score','strand','txstart','txend','rgb','exoncount','exonsize','exonstart')
             pattern="output_controlSet2.*.[bw|bigWig]"
             nc=4
             marks=c('CAGE', 'PAS','H3K4me3', 'Pol2', 'smallRNAseq_76nt_strand_ox', 'EncodeCshlLongRnaSeq', 'RNAseq_PE50nt_strand_R1');
}
if(genetype=='BIP')
{
             ## use controlset
             g=read.table("../data/BIP.all.piNMNR.cpg.bed"); rownames(g)=g[,4]; colnames(g)=c('chr','start','end','name','score','strand','location', 'category', 'normalizedCpG')
             g1=subset(g, category=="pi:pi")
             g2=subset(g, category=="pi:nc")
             g3=subset(g, category=="pi:pc")
             g4=subset(g, location=="NM")
             g5=subset(g, location=="NR")
             pattern="output_BIP.all.piNMNR.*.[bw|bigWig]"
             nc=5
             marks=c('CAGE', 'Degradome.polyA_', 'Degradome.polyA-','H3K4me3', 'Pol2', 'PolIII', 'Amyb', 'RNAseq_PE50nt', 'smallRNAseq');
}
if(genetype=='piRNA'){
             #genes=read.table("../data/10092012.annotation.TSSregion.table", header=T, fill=T)
             genes=read.table("../data/TSS_table_piRNA_adult_mm9.tab", header=T, fill=T)
             rownames(genes)=genes[,4] # TransID
             pattern="output_100312.*.[bw|bigWig]"
             marks=c('CAGE_', 'PAS','H3K4me3', 'Pol2', 'PolIII', 'Amyb', 'EncodeCshlLongRnaSeq', 'RNAseq_PE50nt_strand_R1', 'smallRNAseq_76nt_strand_ox');
}
if(genetype=='uncorrectedpiRNA'){
             #genes=read.table("../data/10092012.annotation.TSSregion.table", header=T, fill=T)
             genes=read.table("../data/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks_o250.transcripts.piRNAuncorrected.bed")
             rownames(genes)=genes[,4]
             colnames(genes)=c('chr','start','end','name','score','strand','txstart','txend','rgb','exoncount','exonsize','exonstart')
             pattern="output_cufflinks_o250.*.[bw|bigWig]"
             marks=c('CAGE_', 'PAS','H3K4me3', 'Pol2', 'PolIII', 'Amyb', 'EncodeCshlLongRnaSeq', 'RNAseq_PE50nt_strand_R1', 'smallRNAseq_76nt_strand_ox');
}
if(genetype=='control')
{
             ## use controlset
             g1=read.table("../data/controlSet/A-Myb.target.bed"); rownames(g1)=g1[,4]; colnames(g1)=c('chr','start','end','name','score','strand','txstart','txend','rgb','exoncount','exonsize','exonstart')
             g2=read.table("../data/controlSet/germline.silenced.bed"); rownames(g2)=g2[,4]; colnames(g2)=c('chr','start','end','name','score','strand','txstart','txend','rgb','exoncount','exonsize','exonstart')
             g3=read.table("../data/controlSet/germline.specifc.nonAMyb.bed"); rownames(g3)=g3[,4]; colnames(g3)=c('chr','start','end','name','score','strand','txstart','txend','rgb','exoncount','exonsize','exonstart')
             g4=read.table("../data/controlSet/housekeeping.bed"); rownames(g4)=g4[,4]; colnames(g4)=c('chr','start','end','name','score','strand','txstart','txend','rgb','exoncount','exonsize','exonstart')
             nc=4
             pattern="output_controlSet.*.[bw|bigWig]"
             marks=c('CAGE_', 'PAS','H3K4me3', 'Pol2', 'PolIII', 'Amyb', 'EncodeCshlLongRnaSeq', 'RNAseq_PE50nt_strand_R1', 'smallRNAseq_76nt_strand_ox');
}

if(genetype=='ensembl')
{
             ## use ensembl genes
             genes=read.table("/home/dongx/nearline/genomes/mm9/Annotation/Genes/NCBIM37.biomart67.transcripts.cpg.tab", header=T)
             rownames(genes)=genes[,6] # TransID
             pattern="output_NCBIM37.*.[bw|bigWig]"
             g1=genes[genes$gene_type=='protein_coding',]
             g2=genes[genes$gene_type=='pseudogene',]
             g3=genes[genes$gene_type=='lincRNA',]
             g4=genes
             marks=c('CAGE_', 'PAS','H3K4me3', 'Pol2', 'PolIII', 'Amyb', 'EncodeCshlLongRnaSeq', 'RNAseq_PE50nt_strand_R1', 'smallRNAseq_76nt_strand_ox');
}


# check outliner!! (which can be interesting biologically)
#`%ni%` <- Negate(`%in%`)
#genes=subset(genes, trans_id %ni% c('ENSMUST00000175219','ENSMUST00000175544','ENSMUST00000175227'))

#hcp=rownames(genes)[genes$normalizedCpG>0.4]
#lcp=rownames(genes)[genes$normalizedCpG<=0.4]

QSUBOUPUT="../data/for_aggregationplot"
filenames = paste(QSUBOUPUT, dir(QSUBOUPUT, pattern=pattern), sep="/")  # for de novo piRNA genes

#png(paste("../result/aggregation",species, paste(marks,collapse="_"), genetype,"png",sep="."), width=1200,height=400*length(marks))
pdf(paste("../result/aggregation",species, genetype,"pdf",sep="."), height=3*length(marks), width=nc*3, colormodel='cmyk')
par(mfrow=c(length(marks),nc))
for(type in marks)
{
             fnames=filenames[grep(paste(type,sep=".*"), filenames, ignore.case = T)]
             print(type)
             print(length(fnames));
             if(length(fnames)==2)
             {
                          ##  ================= all

                          plus=fnames[grep('plus|\\.\\+\\.', fnames, ignore.case =T)][1]
                          minus=fnames[grep('minus|\\.-\\.', fnames, ignore.case =T)][1]

                          if(genetype=='control2') {
                                       mp1=getAGGREGATION2(plus, minus, g1, 0.05);
                                       mp2=getAGGREGATION2(plus, minus, g2, 0.05);
                                       mp3=getAGGREGATION2(plus, minus, g3, 0.05);
                                       mp4=getAGGREGATION2(plus, minus, g4, 0.05);
                                       if(type=='EncodeCshlLongRnaSeq'){
                                              mp1=getAGGREGATION2.cshl(plus, minus, g1, 0.05); mp2=getAGGREGATION2.cshl(plus, minus, g2, 0.05); mp3=getAGGREGATION2.cshl(plus, minus, g3, 0.05); mp4=getAGGREGATION2.cshl(plus, minus, g4, 0.05);
                                       }
                                       #mp5=getAGGREGATION2(plus, minus, g5, 0.05);
                                       draw.plot2(mp1[[1]],mp1[[2]],range(mp1,mp2,mp3,mp4,na.rm=T), ifelse(grepl("CAGE", type), paste("piRNA.genic: n=",formatC(nrow(g1),format="d",big.mark=","), sep=""),""), grepl("CAGE", type), "Mean signal")
                                       draw.plot2(mp2[[1]],mp2[[2]],range(mp1,mp2,mp3,mp4,na.rm=T), ifelse(grepl("CAGE", type), paste("piRNA.intergenic: n=",formatC(nrow(g2),format="d",big.mark=","), sep=""),""))
                                       draw.plot2(mp3[[1]],mp3[[2]],range(mp1,mp2,mp3,mp4,na.rm=T), ifelse(grepl("CAGE", type), paste("protein-coding genes: n=",formatC(nrow(g3),format="d",big.mark=","), sep=""),""))
                                       draw.plot2(mp4[[1]],mp4[[2]],range(mp1,mp2,mp3,mp4,na.rm=T), ifelse(grepl("CAGE", type), paste("ncRNA: n=",formatC(nrow(g4),format="d",big.mark=","), sep=""),""))
                                       #draw.plot2(mp5[[1]],mp5[[2]],range(mp1,mp2,mp3,mp4,mp5,na.rm=T), ifelse(grepl("CAGE", type), paste(type," (pri-miRNA: n=",formatC(nrow(g5),format="d",big.mark=","),")", sep=""),""))
                          }
                          if(genetype=='BIP') {
                                       mp1=getAGGREGATION2(plus, minus, g1, 0.1);
                                       mp2=getAGGREGATION2(plus, minus, g2, 0.1);
                                       mp3=getAGGREGATION2(plus, minus, g3, 0.1);
                                       mp4=getAGGREGATION2(plus, minus, g4, 0.1);
                                       mp5=getAGGREGATION2(plus, minus, g5, 0.1);
                                       draw.plot2(mp1[[1]],mp1[[2]],range(mp1,mp2,mp3,mp4,mp5,na.rm=T), '', grepl("CAGE", type), "Mean signal")
                                       draw.plot2(mp2[[1]],mp2[[2]],range(mp1,mp2,mp3,mp4,mp5,na.rm=T))
                                       draw.plot2(mp3[[1]],mp3[[2]],range(mp1,mp2,mp3,mp4,mp5,na.rm=T))
                                       draw.plot2(mp4[[1]],mp4[[2]],range(mp1,mp2,mp3,mp4,mp5,na.rm=T))
                                       draw.plot2(mp5[[1]],mp5[[2]],range(mp1,mp2,mp3,mp4,mp5,na.rm=T))
                          }
                          if(genetype=='ensembl') {
                                       mp1=getAGGREGATION2(plus, minus, g1);
                                       mp2=getAGGREGATION2(plus, minus, g2);
                                       mp3=getAGGREGATION2(plus, minus, g3);
                                       mp4=getAGGREGATION2(plus, minus, g4);
                                       draw.plot2(mp1[[1]],mp1[[2]],range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (protein_coding: n=",formatC(nrow(g1),format="d",big.mark=","),")", sep=""))
                                       draw.plot2(mp2[[1]],mp2[[2]],range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (pseudogene: n=",formatC(nrow(g2),format="d",big.mark=","),")", sep=""), T)
                                       draw.plot2(mp3[[1]],mp3[[2]],range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (lincRNA: n=",formatC(nrow(g3),format="d",big.mark=","),")", sep=""), T)
                                       draw.plot2(mp4[[1]],mp4[[2]],range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (all: n=",formatC(nrow(g4),format="d",big.mark=","),")", sep=""), T)
                          }
                          if(genetype=='control') {
                                       mp1=getAGGREGATION2(plus, minus, g1);
                                       mp2=getAGGREGATION2(plus, minus, g2);
                                       mp3=getAGGREGATION2(plus, minus, g3);
                                       mp4=getAGGREGATION2(plus, minus, g4);
                                       draw.plot2(mp1[[1]],mp1[[2]],range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (A-Myb.target: n=",formatC(nrow(g1),format="d",big.mark=","),")", sep=""))
                                       draw.plot2(mp2[[1]],mp2[[2]],range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (germline.silenced: n=",formatC(nrow(g2),format="d",big.mark=","),")", sep=""), T)
                                       draw.plot2(mp3[[1]],mp3[[2]],range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (germline.specifc.nonAMyb: n=",formatC(nrow(g3),format="d",big.mark=","),")", sep=""), T)
                                       draw.plot2(mp4[[1]],mp4[[2]],range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (housekeeping: n=",formatC(nrow(g4),format="d",big.mark=","),")", sep=""), T)
                          }
                          if(genetype=='piRNA'){
                                       mp1=getAGGREGATION23(plus, minus, genes)
                                       mp2=getAGGREGATION23(plus, minus, subset(genes, intergenic_genic=='intergenic'))
                                       mp3=getAGGREGATION23(plus, minus, subset(genes, intergenic_genic=='genic'))
                                       draw.plot2(mp1[[1]],mp1[[2]],range(mp1,mp2,mp3,na.rm=T), paste(type," (all: n=",formatC(nrow(genes),format="d",big.mark=","),")", sep=""))
                                       draw.plot2(mp2[[1]],mp2[[2]],range(mp1,mp2,mp3,na.rm=T), paste(type," (intergenic: n=",formatC(sum(genes$intergenic_genic=='intergenic'),format="d",big.mark=","),")", sep=""), T)
                                       draw.plot2(mp3[[1]],mp3[[2]],range(mp1,mp2,mp3,na.rm=T), paste(type," (genic: n=",formatC(sum(genes$intergenic_genic=='genic'),format="d",big.mark=","),")", sep=""), T)
                                       #PLOTAGGREGATION23(plus, minus, genes)
                          }
                          if(genetype=='uncorrectedpiRNA'){
                                       mp1=getAGGREGATION23(plus, minus, genes)
                                       draw.plot2(mp1[[1]],mp1[[2]],range(mp1,na.rm=T), paste(type," (uncorrectedpiRNA: n=",formatC(nrow(genes),format="d",big.mark=","),")", sep=""))
                                       #PLOTAGGREGATION23(plus, minus, genes)
                          }
                          #PLOTAGGREGATION2(plus, minus, genes, 'All : ')
                          #if(genetype=='piRNA'){
                          #             PLOTAGGREGATION2(plus, minus, subset(genes, intergenic_genic=='intergenic'), 'Intergenic : ')
                          #             PLOTAGGREGATION2(plus, minus, subset(genes, intergenic_genic=='genic'), 'Genic : ')
                          #}
             }
             if(length(fnames)==1)
             {
                          if(genetype=='control2') {
                                       mp1=getAGGREGATION(fnames, g1, 0.05);
                                       mp2=getAGGREGATION(fnames, g2, 0.05);
                                       mp3=getAGGREGATION(fnames, g3, 0.05);
                                       mp4=getAGGREGATION(fnames, g4, 0.05);
                                       #mp5=getAGGREGATION(fnames, g5, 0.05);
                                       draw.plot(mp1,range(mp1,mp2,mp3,mp4,0,1,na.rm=T), '',0, "Mean signal")
                                       draw.plot(mp2,range(mp1,mp2,mp3,mp4,0,1,na.rm=T), '')
                                       draw.plot(mp3,range(mp1,mp2,mp3,mp4,0,1,na.rm=T), '')
                                       draw.plot(mp4,range(mp1,mp2,mp3,mp4,0,1,na.rm=T), '')
                                       #draw.plot(mp5,range(mp1,mp2,mp3,mp4,mp5,1,na.rm=T), '')
                          }
                          if(genetype=='BIP') {
                                       mp1=getAGGREGATION(fnames, g1, 0.1);
                                       mp2=getAGGREGATION(fnames, g2, 0.1);
                                       mp3=getAGGREGATION(fnames, g3, 0.1);
                                       mp4=getAGGREGATION(fnames, g4, 0.1);
                                       mp5=getAGGREGATION(fnames, g5, 0.1);
                                       draw.plot(mp1,range(mp1,mp2,mp3,mp4,mp5,0,1,na.rm=T), '', 0, "Mean signal")
                                       draw.plot(mp2,range(mp1,mp2,mp3,mp4,mp5,0,1,na.rm=T))
                                       draw.plot(mp3,range(mp1,mp2,mp3,mp4,mp5,0,1,na.rm=T))
                                       draw.plot(mp4,range(mp1,mp2,mp3,mp4,mp5,0,1,na.rm=T))
                                       draw.plot(mp5,range(mp1,mp2,mp3,mp4,mp5,0,1,na.rm=T))
                          }
                          if(genetype=='control') {
                                       mp1=getAGGREGATION(fnames, g1);
                                       mp2=getAGGREGATION(fnames, g2);
                                       mp3=getAGGREGATION(fnames, g3);
                                       mp4=getAGGREGATION(fnames, g4);
                                       draw.plot(mp1,range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (A-Myb.target: n=",formatC(nrow(g1),format="d",big.mark=","),")", sep=""))
                                       draw.plot(mp2,range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (germline.silenced: n=",formatC(nrow(g2),format="d",big.mark=","),")", sep=""), T)
                                       draw.plot(mp3,range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (germline.specifc.nonAMyb: n=",formatC(nrow(g3),format="d",big.mark=","),")", sep=""), T)
                                       draw.plot(mp4,range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (housekeeping: n=",formatC(nrow(g4),format="d",big.mark=","),")", sep=""), T)
                          }
                          if(genetype=='ensembl') {
                                       mp1=getAGGREGATION(fnames, g1);
                                       mp2=getAGGREGATION(fnames, g2);
                                       mp3=getAGGREGATION(fnames, g3);
                                       mp4=getAGGREGATION(fnames, g4);
                                       draw.plot(mp1,range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (protein_coding: n=",formatC(nrow(g1),format="d",big.mark=","),")", sep=""))
                                       draw.plot(mp2,range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (pseudogene: n=",formatC(nrow(g2),format="d",big.mark=","),")", sep=""), T)
                                       draw.plot(mp3,range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (lincRNA: n=",formatC(nrow(g3),format="d",big.mark=","),")", sep=""), T)
                                       draw.plot(mp4,range(mp1,mp2,mp3,mp4,na.rm=T), paste(type," (all: n=",formatC(nrow(g4),format="d",big.mark=","),")", sep=""), T)
                          }
                          if(genetype=='piRNA'){
                                       mp1=getAGGREGATION3(fnames, genes)
                                       mp2=getAGGREGATION3(fnames, subset(genes, intergenic_genic=='intergenic'))
                                       mp3=getAGGREGATION3(fnames, subset(genes, intergenic_genic=='genic'))
                                       draw.plot(mp1,range(mp1,mp2,mp3,na.rm=T), paste(type," (all: n=",formatC(nrow(genes),format="d",big.mark=","),")", sep=""))
                                       draw.plot(mp2,range(mp1,mp2,mp3,na.rm=T), paste(type," (intergenic: n=",formatC(sum(genes$intergenic_genic=='intergenic'),format="d",big.mark=","),")", sep=""), T)
                                       draw.plot(mp3,range(mp1,mp2,mp3,na.rm=T), paste(type," (genic: n=",formatC(sum(genes$intergenic_genic=='genic'),format="d",big.mark=","),")", sep=""), T)
                                       #PLOTAGGREGATION3(fnames, genes,  0)
                                       #PLOTAGGREGATION(fnames, subset(genes, intergenic_genic=='intergenic'), 'Intergenic : ')
                                       #PLOTAGGREGATION(fnames, subset(genes, intergenic_genic=='genic'), 'Genic : ')
                          }
                          if(genetype=='uncorrectedpiRNA'){
                                       mp1=getAGGREGATION3(fnames, genes)
                                       draw.plot(mp1,range(mp1,na.rm=T), paste(type," (uncorrectedpiRNA: n=",formatC(nrow(genes),format="d",big.mark=","),")", sep=""))
                                       #PLOTAGGREGATION23(plus, minus, genes)
                          }
             }
}
dev.off()

quit("no")

## OLD VERSION [ DEPRECATED ]
QSUBOUPUT="../data"
filenames = paste(QSUBOUPUT, dir(QSUBOUPUT, pattern=pattern), sep="/")  # for de novo piRNA genes

png(paste("../result/aggregation.TSS.TTS",genetype,"png",sep="."), width=1200,height=1200)
par(mfrow=c(3,3))
for(type in c('CAGE','Degradome','PAS'))
{
             plus=filenames[grep(paste(type,'plus',sep=".*"), filenames)]
             minus=filenames[grep(paste(type,'minus',sep=".*"), filenames)]

             histdata=read.table(plus)
             histdata=data.frame(histdata[,-1], row.names=histdata[,1])

             p=histdata[rownames(genes)[genes$strand=='+'], ]
             m=histdata[rownames(genes)[genes$strand=='-'], ]*-1

             histdata=read.table(minus)
             histdata=data.frame(histdata[,-1], row.names=histdata[,1])

             m0=rbind(m,histdata[rownames(genes)[genes$strand=='+'], ])
             p0=rbind(p,histdata[rownames(genes)[genes$strand=='-'], ]*-1)

             # remove trans with nan values from bigWigAverageBed (e.g. chr1 1 1)
             m0=m0[!is.na(rowSums(m0)),] # 79608 out of 79652 remained (for protein-cding gene)
             p0=p0[!is.na(rowSums(p0)),]


             # all genes
             m1=m0; p1=p0;
             # remove gene body bin, and 1% outlier
             n1=nrow(m1)
             m1=apply(m1, 2, function(x) mean(x, trim=0.01)); m1[41]=NA
             p1=apply(p1, 2, function(x) mean(x, trim=0.01)); p1[41]=NA

             # HCP genes
             m2=m0[intersect(rownames(m0),hcp),]; p2=p0[intersect(rownames(m0),hcp),];
             # remove gene body bin, and 1% outlier
             n2=nrow(m2)
             m2=apply(m2, 2, function(x) mean(x, trim=0.01)); m2[41]=NA
             p2=apply(p2, 2, function(x) mean(x, trim=0.01)); p2[41]=NA

             # LCP genes
             m3=m0[intersect(rownames(m0),lcp),]; p3=p0[intersect(rownames(m0),lcp),];
             # remove gene body bin, and 1% outlier
             n3=nrow(m3)
             m3=apply(m3, 2, function(x) mean(x, trim=0.01)); m3[41]=NA
             p3=apply(p3, 2, function(x) mean(x, trim=0.01)); p3[41]=NA

             #plot
             draw.plot2(p1,m1,range(m1,m2,m3, p1,p2,p3, na.rm=T), paste(type," (all genes; n=",formatC(n1,format="d",big.mark=","),")", sep=""))
             draw.plot2(p2,m2,range(m1,m2,m3, p1,p2,p3, na.rm=T), paste(type," (HCP genes; n=",formatC(n2,format="d",big.mark=","),")", sep=""))
             draw.plot2(p3,m3,range(m1,m2,m3, p1,p2,p3, na.rm=T), paste(type," (LCP genes; n=",formatC(n3,format="d",big.mark=","),")", sep=""))

             print(type)
}
dev.off()
