#!/usr/bin/python
import os,sys
import json
import argparse

parser = argparse.ArgumentParser(description="Parse fastp output")
parser.add_argument('--fastp-json','-i',help='Fastp json file',required=True)
parser.add_argument('--outdir','-o',help='Output directory',required=True)
parser.add_argument('--samplename','-s',help='sample name',required=True)
argv = parser.parse_args()

json_file = argv.fastp_json
sample = argv.samplename
outdir = argv.outdir

prefix = sample

def qs2error(qual_score):
	return 10**(-qual_score/10)

def mean(list_data):
	if len(list_data) == 0:
		return 0
	return sum(list_data)/float(len(list_data))

data = json.load(open(json_file))

raw_gc = open(os.path.join(outdir,"raw_%s.GC"%prefix),'w')
raw_qm = open(os.path.join(outdir,"raw_%s.QM"%prefix),'w')
clean_gc = open(os.path.join(outdir,"clean_%s.GC"%prefix),'w')
clean_qm = open(os.path.join(outdir,"clean_%s.QM"%prefix),'w')
stat = open(os.path.join(outdir,"%s.stat"%prefix),'w')

raw_len1 =data['read1_before_filtering']['total_cycles']
raw_len2 =data['read2_before_filtering']['total_cycles']
clean_len1 =data['read1_after_filtering']['total_cycles']
clean_len2 =data['read2_after_filtering']['total_cycles']

## raw left read
for i in range(raw_len1):
	qs = data['read1_before_filtering']['quality_curves']['mean'][i]
	raw_qm.write("%d\t%0.6f\t%0.6f\n" % (i,qs,qs2error(qs)))
	gc_content = data['read1_before_filtering']['content_curves']
	raw_gc.write("%d\tA\t%f\tT\t%f\tG\t%f\tC\t%f\tN\t%f\n" % (i,gc_content['A'][i],gc_content['T'][i], \
		gc_content['G'][i],gc_content['C'][i],gc_content['N'][i]))
## raw right read
for i in range(raw_len2):
	qs = data['read2_before_filtering']['quality_curves']['mean'][i]
	raw_qm.write("%d\t%0.6f\t%0.6f\n" % (raw_len1+i,qs,qs2error(qs)))
	gc_content = data['read2_before_filtering']['content_curves']
	raw_gc.write("%d\tA\t%f\tT\t%f\tG\t%f\tC\t%f\tN\t%f\n" % (raw_len1+i,gc_content['A'][i],gc_content['T'][i], \
		gc_content['G'][i],gc_content['C'][i],gc_content['N'][i]))
## clean left read
for i in range(clean_len1):
	qs = data['read1_after_filtering']['quality_curves']['mean'][i]
	clean_qm.write("%d\t%0.6f\t%0.6f\n" % (i,qs,qs2error(qs)))
	gc_content = data['read1_after_filtering']['content_curves']
	clean_gc.write("%d\tA\t%f\tT\t%f\tG\t%f\tC\t%f\tN\t%f\n" % (i,gc_content['A'][i],gc_content['T'][i], \
		gc_content['G'][i],gc_content['C'][i],gc_content['N'][i]))
## clean right read
for i in range(clean_len2):
	qs = data['read2_after_filtering']['quality_curves']['mean'][i]
	clean_qm.write("%d\t%0.6f\t%0.6f\n" % (clean_len1+i,qs,qs2error(qs)))
	gc_content = data['read2_after_filtering']['content_curves']
	clean_gc.write("%d\tA\t%f\tT\t%f\tG\t%f\tC\t%f\tN\t%f\n" % (clean_len1+i,gc_content['A'][i],gc_content['T'][i], \
		gc_content['G'][i],gc_content['C'][i],gc_content['N'][i]))

q20_raw1 = float(data['read1_before_filtering']['q20_bases'])/(data['read1_before_filtering']['total_bases'])*100
q20_raw2 = float(data['read2_before_filtering']['q20_bases'])/(data['read2_before_filtering']['total_bases'])*100
q20_clean1 = float(data['read1_after_filtering']['q20_bases'])/(data['read1_after_filtering']['total_bases'])*100
q20_clean2 = float(data['read2_after_filtering']['q20_bases'])/(data['read2_after_filtering']['total_bases'])*100
q30_raw1 = float(data['read1_before_filtering']['q30_bases'])/(data['read1_before_filtering']['total_bases'])*100
q30_raw2 = float(data['read2_before_filtering']['q30_bases'])/(data['read2_before_filtering']['total_bases'])*100
q30_clean1 = float(data['read1_after_filtering']['q30_bases'])/(data['read1_after_filtering']['total_bases'])*100
q30_clean2 = float(data['read2_after_filtering']['q30_bases'])/(data['read2_after_filtering']['total_bases'])*100
gc_raw1 = mean(data['read1_before_filtering']['content_curves']['GC'])*100
gc_raw2 = mean(data['read2_before_filtering']['content_curves']['GC'])*100
gc_clean1 = mean(data['read1_after_filtering']['content_curves']['GC'])*100
gc_clean2 = mean(data['read2_after_filtering']['content_curves']['GC'])*100
qs_raw1 = mean(data['read1_before_filtering']['quality_curves']['mean'])
qs_raw2 = mean(data['read2_before_filtering']['quality_curves']['mean'])
qs_clean1 = mean(data['read1_after_filtering']['quality_curves']['mean'])
qs_clean2 = mean(data['read2_after_filtering']['quality_curves']['mean'])

total_reads = data['summary']['before_filtering']['total_reads']
clean_reads = data['summary']['after_filtering']['total_reads']
low_quality_reads = data['filtering_result']['low_quality_reads']
too_many_N_reads = data['filtering_result']['too_many_N_reads']
adapter_reads = data['filtering_result']['too_short_reads']

stat.write("Type\tRaw data\tClean data\n")
stat.write("Number of Reads:\t%d\t%d\n"%(total_reads, clean_reads))
stat.write("Data Size:\t%d\t%d\n"%(data['summary']['before_filtering']['total_bases'],data['summary']['after_filtering']['total_bases']))
stat.write("Number of fq1 Bases:\t%d\t%d\n"%(data['read1_before_filtering']['total_bases'],data['read1_after_filtering']['total_bases']))
stat.write("Number of fq2 Bases:\t%d\t%d\n"%(data['read2_before_filtering']['total_bases'],data['read2_after_filtering']['total_bases']))
stat.write("N of fq1:\t0.00%\t0.00%\n")
stat.write("N of fq2:\t0.00%\t0.00%\n")
stat.write("Low qual base of fq1:(<=5)\t0.00%\t0.00%\n")
stat.write("Low qual base of fq2:(<=5)\t0.00%\t0.00%\n")
stat.write("Q20 of fq1:\t%0.2f%%\t%0.2f%%\n"%(q20_raw1,q20_clean1))
stat.write("Q20 of fq2:\t%0.2f%%\t%0.2f%%\n"%(q20_raw2,q20_clean2))
stat.write("Q30 of fq1:\t%0.2f%%\t%0.2f%%\n"%(q30_raw1,q30_clean1))
stat.write("Q30 of fq2:\t%0.2f%%\t%0.2f%%\n"%(q30_raw2,q30_clean2))
stat.write("GC of fq1:\t%0.2f%%\t%0.2f%%\n"%(gc_raw1,gc_clean1))
stat.write("GC of fq2:\t%0.2f%%\t%0.2f%%\n"%(gc_raw2,gc_clean2))
stat.write("Error of fq1:\t%0.2f%%\t%0.2f%%\n"%(qs2error(qs_raw1)*100,qs2error(qs_clean1)*100))
stat.write("Error of fq2:\t%0.2f%%\t%0.2f%%\n"%(qs2error(qs_raw2)*100,qs2error(qs_clean2)*100))
stat.write("Discard Reads related to N:\t%0.2f%%\n" % (float(too_many_N_reads)/total_reads*100))
stat.write("Discard Reads related to low qual:\t%0.2f%%\n" % (float(low_quality_reads)/total_reads*100))
stat.write("Discard Reads related to Adapter:\t%0.2f%%\n" % (float(adapter_reads)/total_reads*100))
stat.write("Reads Classification:\t%d\t%d\t%d\t%d\t%d\n"%(total_reads,clean_reads,too_many_N_reads,low_quality_reads,adapter_reads))

raw_gc.close()
raw_qm.close()
clean_gc.close()
clean_qm.close()
stat.close()

datasize = float(data['summary']['before_filtering']['total_bases'])/1e9
effective = float(clean_reads)/total_reads*100
err_rate = qs2error(mean([qs_clean1,qs_clean2]))
q20 = float(data['summary']['after_filtering']['q20_bases'])/data['summary']['after_filtering']['total_bases']*100
q30 = float(data['summary']['after_filtering']['q30_bases'])/data['summary']['after_filtering']['total_bases']*100
gc_content = data['summary']['after_filtering']['gc_content']*100
summary = open(os.path.join(outdir,"%s.summary"%prefix),'w')
summary.write("Sample\tRaw reads\tRaw data(G)\tEffective(%s)\tError(%)\tQ20(%)\tQ30(%)\tGC(%)\n")
summary.write("%s\t%d\t%0.1f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n"%(prefix,total_reads,datasize,effective,err_rate,q20,q30,gc_content))
summary.close()

plot_rscript='''setwd("%s")
library(ggplot2)
library(reshape2)
# GC
data<-read.table("clean_%s.GC",sep="\\t",head=F)
df<-as.data.frame(cbind(data[,3],data[,5],data[,7],data[,9],data[,11]))
colnames(df)<-c("A","T","G","C","N")
pos<-data[,1]+1
pos<-as.data.frame(rep(pos,5))
middle=(max(data[,1]+1))/2
mdf<-melt(df,measure=colnames(df))
mdf<-cbind(mdf,pos)
colnames(mdf)<-c("type","percent","pos")
p <- ggplot(mdf, aes(x=pos, y=percent, group=type))+
geom_line(aes(colour = type),size=0.5)+
scale_x_discrete(breaks=seq(0,max(data[,1]+1),20),labels=seq(0,max(data[,1]+1),20))+
ylim(0,50)+
xlab("Position along reads") + ylab("Percent of bases")+
labs(title="Bases content along reads \\n(%s)")+
geom_vline(xintercept = middle,colour="#377EB8",linetype="dashed")
ggsave(filename="%s.GC.pdf",plot=p)
ggsave(filename="%s.GC.png",type="cairo-png",plot=p)
# QM
df<-read.table("clean_%s.QM",sep="\\t",head=F)
colnames(df)=c("pos","Q","E")
df$pos<-df$pos+1
middle=(max(df$pos))/2
p <- ggplot(df,aes(x=pos,y=Q))+
geom_point(size = I(3),colour="#66C2A5")+
geom_line(colour="#66C2A5",size=1)+
xlab("Position along reads") + ylab("Quality score")+
ylim(0,45)+
labs(title="Quality score distribution along reads \\n(%s)")+
scale_x_discrete(breaks=seq(0,max(df$pos),20),labels=seq(0,max(df$pos),20))+
geom_vline(xintercept = middle,colour="#377EB8",linetype="dashed")
ggsave(filename="%s.QM.pdf",plot=p)
ggsave(filename="%s.QM.png",type="cairo-png",plot=p)
# Error
p <- ggplot(df,aes(x=pos,y=E,ymin=0,ymax=E))+
geom_linerange(colour="#66C2A5",size=0.5)+
xlab("Position along reads") + ylab("%% Error rate")+
labs(title="Error rate distribution along reads \\n(%s)")+
geom_vline(xintercept = middle,colour="#377EB8",linetype="dashed")
ggsave(filename="%s.Error.pdf",plot=p)
ggsave(filename="%s.Error.png",type="cairo-png",plot=p)
#pie figure
names<-c("Clean Reads ","Containing N ","Low Quality ","Adapter Related ")
values<-c(%d, %d, %d, %d)
ratio<-sprintf("%%.2f",100*values/sum(values))
regions=c("clean","n","low","adapter")
labels <- paste(names,"(",values,", ",ratio,"%%)",sep="\\t")
df <- data.frame(values = values, type = regions)
p<-ggplot(df,aes(x=factor(0),y = values, fill = regions)) +
geom_bar(width = 1,stat="identity") +
coord_polar(theta="y") +
xlab("") + ylab("") + labs(fill="regions",title = "Classification of Raw Reads \\n(%s)") +
theme(axis.ticks = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_blank(),
panel.grid.minor=element_blank(),
panel.grid.major=element_blank(),
panel.background=element_blank())+
scale_fill_brewer(name="",breaks=regions,labels=labels,palette="Accent")
ggsave(filename="%s.pie3d.pdf", plot=p)
ggsave(filename="%s.pie3d.png", type="cairo-png", plot=p)
'''%(outdir,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,clean_reads,too_many_N_reads,low_quality_reads,adapter_reads,prefix,prefix,prefix)

open(os.path.join(outdir,prefix+'.R'),'w').write(plot_rscript)
assert not os.system("Rscript %s" % os.path.join(outdir,prefix+'.R'))
