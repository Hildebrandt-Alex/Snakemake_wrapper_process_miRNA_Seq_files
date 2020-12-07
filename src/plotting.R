setwd('/naslx/projects/pr62lo/di52dik/projects/exosome_reference_gene_project/data_and_snakefile/snakemake_files/Seq_Data_and_Output')

library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

dir.create("Trimming")
dir.create("Trimming/Untrimmed")
dir.create("Trimming/Trimming_Output")
dir.create("Unmapped")
dir.create("Plots/")
dir.create("Plots/Additional_Files")
dir.create("Readcounts")
dir.create("Readcounts/rRNA")
dir.create("Readcounts/snRNA")
dir.create("Readcounts/snoRNA")
dir.create("Readcounts/tRNA")
dir.create("Readcounts/isomir")


groupnames <- as.character(read.delim("Groupnames.txt",header=F, as.is=T))
sample_names <- gsub(".fastq","",list.files(pattern="*.fastq"))
my_files <- list.files(pattern="*rRNA_readcount.txt")
my_data <- lapply(my_files, function(x) {dat=read.table(x, row.names=NULL, header=F, quote= "",col.names=c("gene",x)); return (dat)})
DT <- lapply(my_data, data.table)
rRNA_readcount <- Reduce(function(...) merge.data.frame(..., all = T, by = c("gene")), DT)
rRNA_readcount <- setnames(rRNA_readcount,c("gene", sample_names))
rRNA_readcount[is.na(rRNA_readcount)] <- 0
write.table(rRNA_readcount,"rRNA_Readcount.txt",row.names=F,col.names=T, quote=F, sep = '\t')
meanrRNA_readcount <- cbind.data.frame(as.data.frame (rRNA_readcount [, grep(paste(groupnames, collapse="|"), colnames(rRNA_readcount), invert=T)]),sapply(groupnames, function(x) rowMeans(as.data.frame (rRNA_readcount) [, grep(x, colnames(rRNA_readcount))] )  ))
meanrRNA_readcount <- setnames(meanrRNA_readcount,1,"gene")
write.table(meanrRNA_readcount,"Mean_rRNA_Readcount.txt",row.names=F,col.names=T, quote=F, sep = '\t')


my_files <- list.files(pattern="*snRNA_readcount.txt")
my_data <- lapply(my_files, function(x) {dat=read.table(x, row.names=NULL, header=F, quote= "",col.names=c("gene",x)); return (dat)})
DT <- lapply(my_data, data.table)
snRNA_readcount <- Reduce(function(...) merge.data.frame(..., all = T, by = c("gene")), DT)
snRNA_readcount <- setnames(snRNA_readcount,c("gene", sample_names))
snRNA_readcount[is.na(snRNA_readcount)] <- 0
write.table(snRNA_readcount,"snRNA_Readcount.txt",row.names=F,col.names=T, quote=F, sep = '\t')
meansnRNA_readcount <- cbind.data.frame(as.data.frame (snRNA_readcount [, grep(paste(groupnames, collapse="|"), colnames(snRNA_readcount), invert=T)]),sapply(groupnames, function(x) rowMeans(as.data.frame (snRNA_readcount) [, grep(x, colnames(snRNA_readcount))] )  ))
meansnRNA_readcount <- setnames(meansnRNA_readcount,1,"gene")
write.table(meansnRNA_readcount,"Mean_snRNA_Readcount.txt",row.names=F,col.names=T, quote=F, sep = '\t')

my_files <- list.files(pattern="*snoRNA_readcount.txt")
my_data <- lapply(my_files, function(x) {dat=read.table(x, row.names=NULL, header=F, quote= "",col.names=c("gene",x)); return (dat)})
DT <- lapply(my_data, data.table)
snoRNA_readcount <- Reduce(function(...) merge.data.frame(..., all = T, by = c("gene")), DT)
snoRNA_readcount <- setnames(snoRNA_readcount,c("gene", sample_names))
snoRNA_readcount[is.na(snoRNA_readcount)] <- 0
write.table(snoRNA_readcount,"snoRNA_Readcount.txt",row.names=F,col.names=T, quote=F, sep = '\t')
meansnoRNA_readcount <- cbind.data.frame(as.data.frame (snoRNA_readcount [, grep(paste(groupnames, collapse="|"), colnames(snoRNA_readcount), invert=T)]),sapply(groupnames, function(x) rowMeans(as.data.frame (snoRNA_readcount) [, grep(x, colnames(snoRNA_readcount))] )  ))
meansnoRNA_readcount <- setnames(meansnoRNA_readcount,1,"gene")
write.table(meansnoRNA_readcount,"Mean_snoRNA_Readcount.txt",row.names=F,col.names=T, quote=F, sep = '\t')

my_files <- list.files(pattern="*tRNA_readcount.txt")
my_data <- lapply(my_files, function(x) {dat=read.table(x, row.names=NULL, header=F, quote= "",col.names=c("gene",x)); return (dat)})
DT <- lapply(my_data, data.table)
tRNA_readcount <- Reduce(function(...) merge.data.frame(..., all = T, by = c("gene")), DT)
tRNA_readcount <- setnames(tRNA_readcount,c("gene", sample_names))
tRNA_readcount[is.na(tRNA_readcount)] <- 0
write.table(tRNA_readcount,"tRNA_Readcount.txt",row.names=F,col.names=T, quote=F, sep = '\t')
meantRNA_readcount <- cbind.data.frame(as.data.frame (tRNA_readcount [, grep(paste(groupnames, collapse="|"), colnames(tRNA_readcount), invert=T)]),sapply(groupnames, function(x) rowMeans(as.data.frame (tRNA_readcount) [, grep(x, colnames(tRNA_readcount))] )  ))
meantRNA_readcount <- setnames(meantRNA_readcount,1,"gene")
write.table(meantRNA_readcount,"Mean_tRNA_Readcount.txt",row.names=F,col.names=T, quote=F, sep = '\t')

my_files <- list.files(pattern="*short_lengths.txt")
my_data <- lapply(my_files, read.table, row.names=NULL, header=F)
DT <- lapply(my_data, data.table)
DT <- lapply(DT, setnames, c("Length","Reads","non-identical sequences")) 
DT <- lapply(DT, function(x) melt(x,id="Length"))
shortlengths <- Reduce(function(...) merge.data.frame(..., all = T, by = c("Length","variable")), DT)
shortlengths[is.na(shortlengths)] <- 0
shortlengths <- setnames(shortlengths,c("Length","Reads", sample_names))
write.table(shortlengths,"Short_Length_Distribution.txt",row.names=F,col.names=T, quote=F, sep = '\t')

meanshortlengths <- cbind.data.frame(as.data.frame (shortlengths) [, grep(paste(groupnames, collapse="|"), colnames(shortlengths), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (shortlengths) [, grep(x, colnames(shortlengths))] )  ))
write.table(meanshortlengths,"Mean_Short_Length_Distribution.txt",row.names=F,col.names=T, quote=F, sep = '\t')

ggplot(melt(meanshortlengths,id.vars=c("Length","Reads")), aes_string(x= "Length", y= "value", colour= "Reads")) +
  geom_line() + labs(title="Short Reads Length Distribution", x = "Sequence Length [nt]", y= "Reads", colour = "Total/ Unique Sequences") + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 20)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10)) + facet_wrap(~variable, nrow=3, scales='free')
ggsave(filename="Short_Length_Distribution.png", width= 14, height = 8, dpi = 600)

ggplot(melt(meanshortlengths,id.vars=c("Length","Reads")), aes_string(x= "Length", y= "value", colour= "Reads")) +
  geom_line() + labs(title="Short Reads Length Distribution", x = "Sequence Length [nt]", y= "Reads", colour = "Total/ Unique Sequences") + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 20)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10)) + facet_wrap(~variable, nrow=3)
ggsave(filename="Short_Length_Distribution_fixed_axis.png", width= 14, height = 8, dpi = 600)

my_files <- list.files(pattern="*short_noadapter.txt")
my_data <- lapply(my_files, read.table, row.names=NULL, header=F, sep=":")
DT <- lapply(my_data, data.table)
short_adapter <- Reduce(function(...) merge.data.frame(..., all = T, by = c("V1")), DT)
short_adapter <- setnames(short_adapter,c("Mapping", sample_names))

my_files <- list.files(pattern="*trimmedreads_lengths.txt")
my_data <- lapply(my_files, read.table, row.names=NULL, header=T, sep="\t")
DT <- lapply(my_data, data.table)
DT <- lapply(DT, function(x) melt(x,id="length"))
trimmedlengths <- Reduce(function(...) merge.data.frame(..., all = T, by = c("length","variable")), DT)
trimmedlengths <- setnames(trimmedlengths,c("Length", "Reads", sample_names))
trimmedlengths[is.na(trimmedlengths)] <- 0

my_files <- list.files(pattern="*rRNA_unmapped_lengths.txt")
my_data <- lapply(my_files, read.table, row.names=NULL, header=T, sep="\t")
DT <- lapply(my_data, data.table)
DT <- lapply(DT, function(x) melt(x,id="length"))
rRNAunmapped <- Reduce(function(...) merge.data.frame(..., all = T, by = c("length","variable")), DT)
rRNAunmapped <- setnames(rRNAunmapped,c("Length", "Reads", sample_names))
rRNAunmapped[is.na(rRNAunmapped)] <- 0
Diff <- merge(trimmedlengths, rRNAunmapped, by=c("Length", "Reads"))
rRNAlengths <- cbind(Diff[,1:2], Mapping="rRNA", Diff[,grepl("*\\.x$",names(Diff))] - Diff[,grepl("*\\.y$",names(Diff))])
rRNAlengths <- setnames(rRNAlengths,c("Length", "Reads", "Mapping", sample_names))

my_files <- list.files(pattern="*tRNA_unmapped_lengths.txt")
my_data <- lapply(my_files, read.table, row.names=NULL, header=T, sep="\t")
DT <- lapply(my_data, data.table)
DT <- lapply(DT, function(x) melt(x,id="length"))
tRNAunmapped <- Reduce(function(...) merge.data.frame(..., all = T, by = c("length","variable")), DT)
tRNAunmapped <- setnames(tRNAunmapped,c("Length", "Reads", sample_names))
tRNAunmapped[is.na(tRNAunmapped)] <- 0
Diff <- merge(rRNAunmapped, tRNAunmapped, by=c("Length", "Reads"))
tRNAlengths <- cbind(Diff[,1:2], Mapping="tRNA", Diff[,grepl("*\\.x$",names(Diff))] - Diff[,grepl("*\\.y$",names(Diff))])
tRNAlengths <- setnames(tRNAlengths,c("Length", "Reads", "Mapping", sample_names))

my_files <- list.files(pattern="*isomir_unmapped_lengths.txt")
my_data <- lapply(my_files, read.table, row.names=NULL, header=T, sep="\t")
DT <- lapply(my_data, data.table)
DT <- lapply(DT, function(x) melt(x,id="length"))
isomirunmapped <- Reduce(function(...) merge.data.frame(..., all = T, by = c("length","variable")), DT)
isomirunmapped <- setnames(isomirunmapped,c("Length", "Reads", sample_names))
isomirunmapped[is.na(isomirunmapped)] <- 0
Diff <- merge(tRNAunmapped, isomirunmapped, by=c("Length", "Reads"))
isomirlengths <- cbind(Diff[,1:2], Mapping="isomir", Diff[,grepl("*\\.x$",names(Diff))] - Diff[,grepl("*\\.y$",names(Diff))])
isomirlengths <- setnames(isomirlengths,c("Length", "Reads", "Mapping", sample_names))

my_files <- list.files(pattern="*snRNA_unmapped_lengths.txt")
my_data <- lapply(my_files, read.table, row.names=NULL, header=T, sep="\t")
DT <- lapply(my_data, data.table)
DT <- lapply(DT, function(x) melt(x,id="length"))
snRNAunmapped <- Reduce(function(...) merge.data.frame(..., all = T, by = c("length","variable")), DT)
snRNAunmapped <- setnames(snRNAunmapped,c("Length", "Reads", sample_names))
snRNAunmapped[is.na(snRNAunmapped)] <- 0
Diff <- merge(isomirunmapped, snRNAunmapped, by=c("Length", "Reads"))
snRNAlengths <- cbind(Diff[,1:2], Mapping="snRNA", Diff[,grepl("*\\.x$",names(Diff))] - Diff[,grepl("*\\.y$",names(Diff))])
snRNAlengths <- setnames(snRNAlengths,c("Length", "Reads", "Mapping", sample_names))

my_files <- list.files(pattern="*snoRNA_unmapped_lengths.txt")
my_data <- lapply(my_files, read.table, row.names=NULL, header=T, sep="\t")
DT <- lapply(my_data, data.table)
DT <- lapply(DT, function(x) melt(x,id="length"))
snoRNAunmapped <- Reduce(function(...) merge.data.frame(..., all = T, by = c("length","variable")), DT)
snoRNAunmapped <- setnames(snoRNAunmapped,c("Length", "Reads", sample_names))
snoRNAunmapped[is.na(snoRNAunmapped)] <- 0
Diff <- merge(snRNAunmapped, snoRNAunmapped, by=c("Length", "Reads"))
snoRNAlengths <- cbind(Diff[,1:2], Mapping="snoRNA", Diff[,grepl("*\\.x$",names(Diff))] - Diff[,grepl("*\\.y$",names(Diff))])
snoRNAlengths <- setnames(snoRNAlengths,c("Length", "Reads", "Mapping", sample_names))

unmappedlengths <- cbind(snoRNAunmapped[,1:2], Mapping="Unmapped",snoRNAunmapped[,-(1:2)])

mappinglengths <- rbind (unmappedlengths, rRNAlengths, snRNAlengths, snoRNAlengths, tRNAlengths, isomirlengths)
write.table(mappinglengths,"Mapping_Length_Distribution.txt",row.names=F,col.names=T, quote=F, sep = '\t')

meanmappinglengths <- cbind.data.frame(as.data.frame (mappinglengths) [, grep(paste(groupnames, collapse="|"), colnames(mappinglengths), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (mappinglengths) [, grep(x, colnames(mappinglengths))] )  ))
write.table(meanmappinglengths,"Mean_Mapping_Length_Distribution.txt",row.names=F,col.names=T, quote=F, sep = '\t')

ggplot(melt(meanmappinglengths[meanmappinglengths$Reads == "reads",],id.vars=c("Length","Reads", "Mapping")), aes_string(x = "Length", y = "value", fill = "Mapping")) + geom_area(position = 'stack') + facet_wrap(~variable, nrow=3, scales='free') + labs(title="Mapping Length Distribution", x = "Sequence Length [nt]", y= "Reads", fill = "Mapping") + scale_fill_manual(values = c("#4DAF4A", "#984EA3", "#A65628", "#FFFF33",  "#FF7F00", "#E41A1C")) + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 20)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mapping_Length_Distribution.png", width= 14, height = 8, dpi = 600)

ggplot(melt(meanmappinglengths[meanmappinglengths$Reads == "reads",],id.vars=c("Length","Reads", "Mapping")), aes_string(x = "Length", y = "value", fill = "Mapping")) + geom_area(position = 'stack') + facet_wrap(~variable, nrow=3) + labs(title="Mapping Length Distribution", x = "Sequence Length [nt]", y= "Reads", fill = "Mapping") + scale_fill_manual(values = c("#4DAF4A", "#984EA3", "#A65628", "#FFFF33",  "#FF7F00", "#E41A1C")) + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 20)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mapping_Length_Distribution_fixed_axis.png", width= 14, height = 8, dpi = 600)

ggplot(melt(meanmappinglengths[meanmappinglengths$Reads == "non.identical.sequences",],id.vars=c("Length","Reads", "Mapping")), aes_string(x = "Length", y = "value", fill = "Mapping")) + geom_area(position = 'stack') + facet_wrap(~variable, nrow=3, scales='free') + labs(title="Mapping Length Unique Sequences Distribution", x = "Sequence Length [nt]", y= "Reads", fill = "Mapping") + scale_fill_manual(values = c("#4DAF4A", "#984EA3", "#A65628", "#FFFF33",  "#FF7F00", "#E41A1C")) + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 20)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mapping_Length_Unique_Sequences_Distribution.png", width= 14, height = 8, dpi = 600)

ggplot(melt(meanmappinglengths[meanmappinglengths$Reads == "non.identical.sequences",],id.vars=c("Length","Reads", "Mapping")), aes_string(x = "Length", y = "value", fill = "Mapping")) + geom_area(position = 'stack') + facet_wrap(~variable, nrow=3) + labs(title="Mapping Length Unique Sequences Distribution", x = "Sequence Length [nt]", y= "Reads", fill = "Mapping") + scale_fill_manual(values = c("#4DAF4A", "#984EA3", "#A65628", "#FFFF33",  "#FF7F00", "#E41A1C")) + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 20)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mapping_Length_Unique_Sequences_Distribution_fixed_axis.png", width= 14, height = 8, dpi = 600)

for (i in colnames(meanmappinglengths[-(1:3)])) {
  print(ggplot(data=meanmappinglengths,
               aes_string(x= "Length", y= as.name(i), colour= "Reads")) +
          geom_line() + labs(title="Mapping Length Distribution", x = "Sequence Length [nt]", y= "Reads", colour = "Total/ Unique Sequences") + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 20)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10)) + facet_wrap(~Mapping, nrow=3, scales='free')) 
  ggsave(filename=paste(i,"_Mapping_Length_Unique_Distribution.png",sep=""), width= 10, height = 8, dpi = 600)
}

for (i in colnames(meanmappinglengths[-(1:3)])) {
  print(ggplot(data=meanmappinglengths,
               aes_string(x= "Length", y= as.name(i), colour= "Reads")) +
          geom_line() + labs(title="Mapping Length Distribution", x = "Sequence Length [nt]", y= "Reads", colour = "Total/ Unique Sequences") + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 20)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10)) + facet_wrap(~Mapping, nrow=3)) 
  ggsave(filename=paste(i,"_Mapping_Length_Unique_Distribution_fixed_axis.png",sep=""), width= 10, height = 8, dpi = 600)
}

mapping <- data.table(copy(mappinglengths[mappinglengths$Reads == "reads",]))
mapping[,(1:2):=NULL]
setkey (mapping, Mapping)
mapping <-  mapping[,lapply(.SD,sum),by = key(mapping)] 
mapping <- rbind (short_adapter, mapping)
write.table(mapping,"Mapping_Distribution.txt",row.names=F,col.names=T, quote=F, sep = '\t')
ggplot(melt(mapping, id.vars="Mapping"), aes_string(x = "variable", y = "value", fill = "Mapping")) + geom_bar(stat = "identity") + labs(title="Mapping Distribution", x = "Samples", y= "Reads", fill = "Mapping") + scale_fill_manual(values = c("#F781BF",  "#377EB8", "#4DAF4A", "#984EA3", "#A65628", "#FFFF33",  "#FF7F00", "#E41A1C")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mapping_Distribution.png", width= (ifelse(length(sample_names)*0.5 <= 5, 5, (length(sample_names)*0.5))), dpi = 600, limitsize = FALSE)

meanmapping <- cbind.data.frame(as.data.frame (mapping)[, grep(paste(groupnames, collapse="|"), colnames(mapping), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (mapping) [, grep(x, colnames(mapping))] )))
meanmapping <- setnames(meanmapping,1,"Mapping")
write.table(meanmapping,"Mean_Mapping_Distribution.txt",row.names=F,col.names=T, quote=F, sep = '\t')
ggplot(melt(meanmapping, id.vars="Mapping"), aes_string(x = "variable", y = "value", fill = "Mapping")) + geom_bar(stat = "identity") + labs(title="Mapping Distribution", x = "Samples", y= "Reads", fill = "Mapping") + scale_fill_manual(values = c("#F781BF",  "#377EB8", "#4DAF4A", "#984EA3", "#A65628", "#FFFF33",  "#FF7F00", "#E41A1C")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mean_Mapping_Distribution.png", width= (ifelse(length(groupnames)*0.5 <= 5, 5, (length(groupnames)*0.5))), dpi = 600)

ggplot(melt(mapping, id.vars="Mapping"), aes_string(x = "variable", y = "value", fill = "Mapping")) + geom_bar(stat = "identity", position ="fill") + labs(title="Mapping Distribution", x = "Samples", y= "Relative frequency", fill = "Mapping") + scale_fill_manual(values = c("#F781BF",  "#377EB8", "#4DAF4A", "#984EA3", "#A65628", "#FFFF33",  "#FF7F00", "#E41A1C")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Relative_Mapping_Distribution.png", width= (ifelse(length(sample_names)*0.5 <= 5, 5, (length(sample_names)*0.5))), dpi = 600)

ggplot(melt(meanmapping, id.vars="Mapping"), aes_string(x = "variable", y = "value", fill = "Mapping")) + geom_bar(stat = "identity", position = "fill") + labs(title="Mapping Distribution", x = "Samples", y= "Relative frequency", fill = "Mapping") + scale_fill_manual(values = c("#F781BF",  "#377EB8", "#4DAF4A", "#984EA3", "#A65628", "#FFFF33",  "#FF7F00", "#E41A1C")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mean_Relative_Mapping_Distribution.png", width= (ifelse(length(groupnames)*0.5 <= 5, 5, (length(groupnames)*0.5))), dpi = 600)

my_files_i <- list.files(pattern="*isomir_readcount.txt")
my_data_i <- lapply(my_files_i, function(x) {dat=read.table(x, row.names=NULL, header=F, quote= "",col.names=c("gene",x)); return (dat)})
DT_i <- lapply(my_data_i, data.table)
isomir_readcount <- Reduce(function(...) merge(..., all = T, by = c("gene")), DT_i)
isomir_readcount <- setnames(isomir_readcount,c("isomir", sample_names))
isomir_readcount[is.na(isomir_readcount)] <- 0
write.table(isomir_readcount,"isomir_Readcount.txt",row.names=F,col.names=T, quote=F, sep = '\t')
meanisomir_readcount <- cbind.data.frame(as.data.frame (isomir_readcount) [, grep(paste(groupnames, collapse="|"), colnames(isomir_readcount), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (isomir_readcount) [, grep(x, colnames(isomir_readcount))] )  ))
meanisomir_readcount <- setnames(meanisomir_readcount,1,"isomir")
write.table(meanisomir_readcount,"Mean_isomir_Readcount.txt",row.names=F,col.names=T, quote=F, sep = '\t')
#write.table(meanisomir_readcount,snakemake@output[[1]],row.names=F,col.names=T, quote=F, sep = '\t')

multiread_readcount <- isomir_readcount[like (isomir,"multi")]
write.table(multiread_readcount,"multiread_Readcount.txt",row.names=F,col.names=T, quote=F, sep = '\t')
meanisomir_readcount <- data.table(meanisomir_readcount)
meanmultiread_readcount <- meanisomir_readcount[like (isomir,"multi")]
write.table(meanmultiread_readcount,"Mean_multiread_Readcount.txt",row.names=F,col.names=T, quote=F, sep = '\t')

isomir_readcount[, c("canonical","mod5","nt5","seq5","template5","mod3","nt3","seq3","template3","multiread","missmatch","ntmm","posmm","seed","seedseq","m8","m8seq","end3","end3seq","length","ntlength") := tstrsplit(isomir, "_", fixed=TRUE)]
isomir_readcount[,1:=NULL]
setcolorder(isomir_readcount, c("canonical","mod5","nt5","seq5","template5","mod3","nt3","seq3","template3","multiread","missmatch","ntmm","posmm","seed","seedseq","m8","m8seq","end3","end3seq","length","ntlength", sample_names))
isomir_readcount$ntlength <- as.integer(isomir_readcount$ntlength)
isomir_readcount$nt5 <- as.integer(isomir_readcount$nt5)
isomir_readcount$nt3 <- as.integer(isomir_readcount$nt3)

parentmiRNAsum <- copy(isomir_readcount)
parentmiRNAsum[,(2:21):=NULL]
parentmiRNAsum <- parentmiRNAsum[,lapply(.SD,sum),by = canonical]
write.table(parentmiRNAsum,"Parent_miRNA_Sum_Readcount.txt",row.names=F, sep = '\t')

iso1mmsum <- copy(isomir_readcount)
iso1mmsum <- iso1mmsum[ntmm <2]
write.table(data.frame(isomir = do.call(paste, c(iso1mmsum[ , 1:21], list(sep = '_'))), iso1mmsum[,-(1:21)]),"isomir_1mm_Readcount.txt",row.names = F, sep = '\t') 
iso1mmsum[,(2:21):=NULL]
iso1mmsum <- iso1mmsum[,lapply(.SD,sum),by = canonical]
write.table(iso1mmsum,"Parent_miRNA_1mm_Sum_Readcount.txt",row.names=F, sep = '\t')

iso3ntsum <- copy(isomir_readcount)
iso3ntsum <- iso3ntsum[nt5 <4 & nt5 >-4 & nt3 <4 &nt3 > -4]
iso3ntsum[,(2:21):=NULL]
iso3ntsum <- iso3ntsum[,lapply(.SD,sum),by = canonical]
write.table(iso3ntsum,"Parent_miRNA_3nt_Sum_Readcount.txt",row.names=F, sep = '\t')

isounisum <- copy(isomir_readcount)
isounisum <- isounisum[multiread == "unique"]
isounisum[,(2:21):=NULL]
isounisum <- isounisum[,lapply(.SD,sum),by = canonical]
write.table(isounisum,"Parent_miRNA_Unique_Sum_Readcount.txt",row.names=F, sep = '\t')

iso3nt1mmunisum <- copy(isomir_readcount)
iso3nt1mmunisum <- iso3nt1mmunisum[nt5 <4 & nt5 >-4 & nt3 <4 &nt3 > -4 & ntmm < 2 & multiread == "unique"]
iso3nt1mmunisum[,(2:21):=NULL]
iso3nt1mmunisum <- iso3nt1mmunisum[,lapply(.SD,sum),by = canonical]
write.table(iso3nt1mmunisum,"Parent_miRNA_3nt_1mm_Unique_Sum_Readcount.txt",row.names=F, sep = '\t')

canonicalsum <- copy(isomir_readcount)
canonicalsum <- canonicalsum[mod5 == "Can5" & mod3 == "Can3" & missmatch == "0mm"]
canonicalsum[,(2:21):=NULL]
canonicalsum <- canonicalsum[,lapply(.SD,sum),by = canonical]
write.table(canonicalsum,"Canonical_Strict_Sum_Readcount.txt",row.names=F, sep = '\t')

maturesum <- copy(isomir_readcount)
maturesum <- maturesum[(mod5 == "Can5" | mod5  == "Trim5") & (mod3 == "Can3" | mod3 == "Trim3")]
maturesum[,(2:21):=NULL]
maturesum <- maturesum[,lapply(.SD,sum),by = canonical]
write.table(maturesum,"Mature_miRNA_Sum_Readcount.txt",row.names=F, sep = '\t')

maturesum <- copy(isomir_readcount)
maturesum <- maturesum[(mod5 == "Can5" | mod5  == "Trim5") & (mod3 == "Can3" | mod3 == "Trim3") & (ntmm < 2)]
maturesum[,(2:21):=NULL]
maturesum <- maturesum[,lapply(.SD,sum),by = canonical]
write.table(maturesum,"Mature_miRNA_1mm_Sum_Readcount.txt",row.names=F, sep = '\t')

precursorsum <- copy(isomir_readcount)
precursorsum <- precursorsum[(template5 == "template") & (template3 == "template")]
precursorsum[,(2:21):=NULL]
precursorsum <- precursorsum[,lapply(.SD,sum),by = canonical]
write.table(precursorsum,"Precursor_Sum_Readcount.txt",row.names=F, sep = '\t')

precursor1mmsum <- copy(isomir_readcount)
precursor1mmsum <- precursor1mmsum[(template5 == "template") & (template3 == "template") & ntmm <2 ]
precursor1mmsum[,(2:21):=NULL]
precursor1mmsum <- precursor1mmsum[,lapply(.SD,sum),by = canonical]
write.table(precursor1mmsum,"Precursor_1mm_Sum_Readcount.txt",row.names=F, sep = '\t')

parentmiRNAsum <- copy(isomir_readcount)
parentmiRNAsum[,(2:11):=NULL]
parentmiRNAsum[,(3:11):=NULL]
parentmiRNAsum[, canonical := "isomir_total"]
parentmiRNAsum <- parentmiRNAsum[,lapply(.SD,sum),by = c("canonical","ntmm")]

canonicalsum <- copy(isomir_readcount)
canonicalsum <- canonicalsum[mod5 == "Can5" & mod3 == "Can3" & missmatch == "0mm"]
canonicalsum[,(2:11):=NULL]
canonicalsum[,(3:11):=NULL]
canonicalsum[, canonical := "canonical_strict"]
canonicalsum <- canonicalsum[,lapply(.SD,sum),by = c("canonical","ntmm")]

maturesum <- copy(isomir_readcount)
maturesum <- maturesum[(mod5 == "Can5" | mod5  == "Trim5") & (mod3 == "Can3" | mod3 == "Trim3")]
maturesum[,(2:11):=NULL]
maturesum[,(3:11):=NULL]
maturesum[, canonical := "mature_miRNA"]
maturesum <- maturesum[,lapply(.SD,sum),by = c("canonical","ntmm")]

precursorsum <- copy(isomir_readcount)
precursorsum <- precursorsum[(template5 == "template") & (template3 == "template")]
precursorsum[,(2:11):=NULL]
precursorsum[,(3:11):=NULL]
precursorsum[, canonical := "precursor_miRNA"]
precursorsum <- precursorsum[,lapply(.SD,sum),by = c("canonical","ntmm")]

isomirmapping <- rbind.data.frame(canonicalsum,maturesum,precursorsum,parentmiRNAsum)
write.table(isomirmapping,"Mapping_Yield_Missmatch_Distribution.txt",row.names=F, sep = '\t')
meanisomirmapping <- cbind.data.frame(as.data.frame (isomirmapping) [, grep(paste(groupnames, collapse="|"), colnames(isomirmapping), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (isomirmapping) [, grep(x, colnames(isomirmapping))] )  ))
meanisomirmapping$canonical <- factor(meanisomirmapping$canonical, levels = c("canonical_strict", "mature_miRNA", "precursor_miRNA", "isomir_total"))
meanisomirmapping$ntmm <- factor(meanisomirmapping$ntmm, levels = c("3", "2", "1", "0"))
write.table(meanisomirmapping,"Mean_Mapping_Yield_Missmatch_Distribution.txt",row.names=F, sep = '\t')

ggplot(melt(meanisomirmapping,id.vars=c("canonical","ntmm")), aes_string(x = "canonical", y = "value", fill = "ntmm")) + geom_bar(stat="identity") + facet_wrap(~variable, scales = "free") + labs(title="Mapping Yield Distribution", y= "Reads", fill = "Number of Missmatches", x= "isomir Classes") + scale_fill_manual(values = c("#4DAF4A", "#FF7F00", "#377EB8", "#E41A1C")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mapping_Yield_Missmatch_Distribution.png", width= 14, height = 8, dpi = 600)

ggplot(melt(meanisomirmapping,id.vars=c("canonical","ntmm")), aes_string(x = "canonical", y = "value", fill = "ntmm")) + geom_bar(stat="identity") + facet_wrap(~variable) + labs(title="Mapping Yield Distribution", y= "Reads", fill = "Number of Missmatches", x= "isomir Classes") + scale_fill_manual(values = c("#4DAF4A", "#FF7F00", "#377EB8", "#E41A1C")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mapping_Yield_Missmatch_Distribution_fixed_axis.png", width= 14, height = 8, dpi = 600)

parentmiRNAsum <- copy(isomir_readcount)
parentmiRNAsum[,(2:9):=NULL]
parentmiRNAsum[,(3:13):=NULL]
parentmiRNAsum[, canonical := "isomir_total"]
parentmiRNAsum <- parentmiRNAsum[,lapply(.SD,sum),by = c("canonical","multiread")]

canonicalsum <- copy(isomir_readcount)
canonicalsum <- canonicalsum[mod5 == "Can5" & mod3 == "Can3" & missmatch == "0mm"]
canonicalsum[,(2:9):=NULL]
canonicalsum[,(3:13):=NULL]
canonicalsum[, canonical := "canonical_strict"]
canonicalsum <- canonicalsum[,lapply(.SD,sum),by = c("canonical","multiread")]

maturesum <- copy(isomir_readcount)
maturesum <- maturesum[(mod5 == "Can5" | mod5  == "Trim5") & (mod3 == "Can3" | mod3 == "Trim3")]
maturesum[,(2:9):=NULL]
maturesum[,(3:13):=NULL]
maturesum[, canonical := "mature_miRNA"]
maturesum <- maturesum[,lapply(.SD,sum),by = c("canonical","multiread")]

precursorsum <- copy(isomir_readcount)
precursorsum <- precursorsum[(template5 == "template") & (template3 == "template")]
precursorsum[,(2:9):=NULL]
precursorsum[,(3:13):=NULL]
precursorsum[, canonical := "precursor_miRNA"]
precursorsum <- precursorsum[,lapply(.SD,sum),by = c("canonical","multiread")]

isomirmapping <- rbind.data.frame(canonicalsum,maturesum,precursorsum,parentmiRNAsum)
write.table(isomirmapping,"Mapping_Yield_Unique_Distribution.txt",row.names=F, sep = '\t')
meanisomirmapping <- cbind.data.frame(as.data.frame (isomirmapping) [, grep(paste(groupnames, collapse="|"), colnames(isomirmapping), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (isomirmapping) [, grep(x, colnames(isomirmapping))] )  ))
meanisomirmapping$canonical <- factor(meanisomirmapping$canonical, levels = c("canonical_strict", "mature_miRNA", "precursor_miRNA", "isomir_total"))
write.table(meanisomirmapping,"Mean_Mapping_Yield_Unique_Distribution.txt",row.names=F, sep = '\t')

ggplot(melt(meanisomirmapping,id.vars=c("canonical","multiread")), aes_string(x = "canonical", y = "value", fill = "multiread")) + geom_bar(stat="identity") + facet_wrap(~variable, scales ="free") + labs(title="Mapping Yield Distribution", y= "Reads", fill = "Multiread", x= "isomir Classes") + scale_fill_manual(values =c("#377EB8", "#E41A1C")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mapping_Yield_Unique_Distribution.png", width= 14, height = 8, dpi = 600)

ggplot(melt(meanisomirmapping,id.vars=c("canonical","multiread")), aes_string(x = "canonical", y = "value", fill = "multiread")) + geom_bar(stat="identity") + facet_wrap(~variable) + labs(title="Mapping Yield Distribution", y= "Reads", fill = "Multiread", x= "isomir Classes") + scale_fill_manual(values =c("#377EB8", "#E41A1C")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mapping_Yield_Unique_Distribution_fixed_axis.png", width= 14, height = 8, dpi = 600)

lengthsum <- copy(isomir_readcount)  
lengthsum[,(1:11):=NULL]
lengthsum[,(2:9):=NULL]
setkey (lengthsum, ntmm, ntlength)
lengthsum <-  lengthsum[,lapply(.SD,sum),by = key(lengthsum)]
write.table(lengthsum[order(ntlength),],"Isomir_Length_Missmatch_Distribution_Overview.txt",row.names=F, sep = '\t')
meanlengthsum <- cbind.data.frame(as.data.frame (lengthsum) [, grep(paste(groupnames, collapse="|"), colnames(lengthsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (lengthsum) [, grep(x, colnames(lengthsum))] )  ))
write.table(meanlengthsum,"Mean_Isomir_Length_Missmatch_Distribution_Overview.txt",row.names=F, sep = '\t')
meanlengthsum$ntmm <- factor(meanlengthsum$ntmm)
meanlengthsum$ntmm <- factor(meanlengthsum$ntmm, levels = rev(levels(meanlengthsum$ntmm)))
meanlengthsum$ntlength <- as.integer(meanlengthsum$ntlength)
ggplot(melt(meanlengthsum,id.vars=c("ntmm","ntlength")), aes_string(x = "ntlength", y = "value", fill = "ntmm")) + geom_area(position = 'stack') + facet_wrap(~variable, nrow=3, scales = 'free') + labs(title="Isomir Length Missmatch Distribution", x = "Sequence Length [nt]", y= "Reads", fill = "Number of Missmatches") + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 20)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10)) + scale_fill_manual(values = c("#4DAF4A", "#FF7F00", "#377EB8", "#E41A1C"))
ggsave(filename="Isomir_Length_Missmatch_Distribution.png", width= 12, dpi = 600)
for (i in colnames(meanlengthsum[-(1:2)])) {
  plot <- print(ggplot(meanlengthsum, aes_string(x = "ntlength", y = as.name(i), fill = "ntmm")) + geom_area(position = 'stack') + labs(title= paste(i, " Isomir Length Missmatch Distribution",sep=""), x = "Sequence Length [nt]", y= "Reads", fill = "Number of Missmatches") + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 20)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10)) + scale_fill_manual(values = c("#4DAF4A", "#FF7F00", "#377EB8", "#E41A1C")))
  ggsave(filename=paste(i,"_Isomir_Length_Missmatch_Distribution.png",sep=""), dpi = 600)
}
lengthsum[,(1):=NULL]
lengthsum <- lengthsum[,lapply(.SD,sum),by = ntlength]
lengthsum <- lengthsum[ order(lengthsum$ntlength),]
write.table(lengthsum,"Isomir_Length_Distribution_Overview.txt",row.names=F, sep = '\t')
meanlengthsum <- cbind.data.frame(as.data.frame (lengthsum) [, grep(paste(groupnames, collapse="|"), colnames(lengthsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (lengthsum) [, grep(x, colnames(lengthsum))] )  ))
meanlengthsum <- setnames(meanlengthsum,1,"ntlength")
write.table(meanlengthsum,"Mean_Isomir_Length_Distribution_Overview.txt",row.names=F, sep = '\t')
ggplot(data=melt(meanlengthsum,id="ntlength"),
       aes(x=ntlength, y=value, colour=variable)) +
  geom_line(size= 1.5) + labs(title="Isomir Length Distribution", x = "Sequence Length [nt]", y= "Reads", colour = "Groups")
ggsave(filename="Isomir_Length_Distribution.png", dpi = 600)

mod5sum <- copy(isomir_readcount)
mod5sum[,(1):=NULL]
mod5sum[,(5:20):=NULL]
setkey (mod5sum, mod5, nt5, seq5, template5)
mod5sum <-  mod5sum[,lapply(.SD,sum),by = key(mod5sum)] 
write.table(mod5sum,"5'end_Modifications_Sum_Template_Overview.txt",row.names=F, sep = '\t')
meanmod5sum <- cbind.data.frame(as.data.frame (mod5sum) [, grep(paste(groupnames, collapse="|"), colnames(mod5sum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (mod5sum) [, grep(x, colnames(mod5sum))] )  ))
write.table(meanmod5sum,"Mean_5'end_Modifications_Sum_Template_Overview.txt",row.names=F, sep = '\t')
mod5sum <- data.table(mod5sum)
mod5sum <- mod5sum[nt5 <4 & nt5 >-4]
mod5sum[,(2):=NULL]
meanmod5sum <- cbind.data.frame(as.data.frame (mod5sum) [, grep(paste(groupnames, collapse="|"), colnames(mod5sum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (mod5sum) [, grep(x, colnames(mod5sum))] )  ))
ggplot(melt(meanmod5sum,id.vars=c("mod5","seq5","template5")),aes(seq5,value,fill=variable))+
  geom_bar(position="dodge",stat="identity")+
  facet_wrap(~template5 + mod5,nrow=4, scales = 'free') + labs(title="5'end Modification Nucleotide Distribution", x = "5' Sequence", y= "Reads", fill = "Groups") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename="5'end_Modification_Nucleotide_Distribution_3nt.png", dpi = 600, height = 8 ,width= 16, limitsize = F)

mod5sum[,(3):=NULL]
setkey (mod5sum, mod5, seq5)
mod5sum <-  mod5sum[,lapply(.SD,sum),by = key(mod5sum)] 
write.table(mod5sum,"5'end_Modifications_Sum_Overview.txt",row.names=F, sep = '\t') 
meanmod5sum <- cbind.data.frame(as.data.frame (mod5sum) [, grep(paste(groupnames, collapse="|"), colnames(mod5sum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (mod5sum) [, grep(x, colnames(mod5sum))] )  ))
write.table(meanmod5sum,"Mean_5'end_Modifications_Sum_Overview.txt",row.names=F, sep = '\t') 
mod5sum <- copy(isomir_readcount)
mod5sum[,(1:2):=NULL]
mod5sum[,(2:19):=NULL]
mod5sum <- mod5sum[,lapply(.SD,sum),by = nt5]
meanmod5sum <- cbind.data.frame(as.data.frame (mod5sum) [, grep(paste(groupnames, collapse="|"), colnames(mod5sum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (mod5sum) [, grep(x, colnames(mod5sum))] )  ))
meanmod5sum <- setnames(meanmod5sum,1,"nt5")

ggplot(melt(meanmod5sum,id.vars=c("nt5")),aes(nt5,value,color=variable)) +
  geom_line(size= 1.5) + labs(title="5'end Modification Length Distribution", x = "5' Modfication Length", y= "Reads", color = "Groups") + scale_x_reverse(labels = comma, breaks = c(-6:6))
ggsave(filename="5'end_Modification_Length_Distribution_Lines.png", dpi = 600)

meanmod5sum$nt5 <- factor(meanmod5sum$nt5, levels = c("3", "2", "1", "0", "-1", "-2", "-3", "-4", "-5", "-6"))
ggplot(melt(meanmod5sum,id.vars=c("nt5")),aes(nt5,value,fill=variable))+
  geom_bar(position="dodge",stat="identity") + labs(title="5'end Modification Length Distribution", x = "5' Modfication Length", y= "Reads", fill = "Groups")
ggsave(filename="5'end_Modification_Length_Distribution_Bars.png", dpi = 600)


mod3sum <- copy(isomir_readcount)
mod3sum[,(1:5):=NULL]
mod3sum[,(5:16):=NULL]
setkey (mod3sum, mod3, nt3, seq3, template3)
mod3sum <-  mod3sum[,lapply(.SD,sum),by = key(mod3sum)] 
write.table(mod3sum,"3'end_Modifications_Sum_Template_Overview.txt",row.names=F, sep = '\t') 
mod3sum <- data.table(mod3sum)
mod3sum <- mod3sum[nt3 <4 & nt3 >-4]
mod3sum[,(2):=NULL]
meanmod3sum <- cbind.data.frame(as.data.frame (mod3sum) [, grep(paste(groupnames, collapse="|"), colnames(mod3sum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (mod3sum) [, grep(x, colnames(mod3sum))] )  ))
write.table(meanmod3sum,"Mean_3'end_Modifications_Sum_Template_Overview.txt",row.names=F, sep = '\t') 
ggplot(melt(meanmod3sum,id.vars=c("mod3","seq3","template3")),aes(seq3,value,fill=variable))+
  geom_bar(position="dodge",stat="identity") +
  facet_wrap(~template3 + mod3,nrow=4, scales = 'free') + labs(title="3'end Modification Nucleotide Distribution", x = "3' Sequence", y= "Reads", fill = "Groups") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename="3'end_Modification_Nucleotide_Distribution.png", dpi = 600, height=8, width= 16, limitsize = F)

mod3sum[,(3):=NULL]
setkey (mod3sum, mod3, seq3)
mod3sum <-  mod3sum[,lapply(.SD,sum),by = key(mod3sum)] 
write.table(mod3sum,"3'end_Modifications_Sum_Overview.txt",row.names=F, sep = '\t') 
meanmod3sum <- cbind.data.frame(as.data.frame (mod3sum) [, grep(paste(groupnames, collapse="|"), colnames(mod3sum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (mod3sum) [, grep(x, colnames(mod3sum))] )  ))
write.table(meanmod3sum,"Mean_3'end_Modifications_Sum_Overview.txt",row.names=F, sep = '\t') 
mod3sum <- copy(isomir_readcount)
mod3sum[,(1:6):=NULL]
mod3sum[,(2:15):=NULL]
setkey (mod3sum, nt3)
mod3sum <- mod3sum[,lapply(.SD,sum),by = nt3]
meanmod3sum <- cbind.data.frame(as.data.frame (mod3sum) [, grep(paste(groupnames, collapse="|"), colnames(mod3sum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (mod3sum) [, grep(x, colnames(mod3sum))] )  ))
meanmod3sum <- setnames(meanmod3sum,1,"nt3")

ggplot(melt(meanmod3sum,id.vars=c("nt3")),aes(nt3,value,color=variable))+
  geom_line(size= 1.5) + labs(title="3'end Modification Length Distribution", x = "3' Modfication Length", y= "Reads", color = "Groups") + scale_x_continuous(labels = comma, breaks = c(-6:6)) 
ggsave(filename="3'end_Modification_Length_Distribution_Lines.png", dpi = 600)

ggplot(melt(meanmod3sum,id.vars=c("nt3")),aes(nt3,value,fill=variable))+
  geom_bar(position="dodge",stat="identity") + labs(title="3'end Modification Length Distribution", x = "3' Modfication Length", y= "Reads", fill = "Groups") + scale_x_continuous(labels = comma, breaks = c(-6:6)) 
ggsave(filename="3'end_Modification_Length_Distribution_Bars.png", dpi = 600)


templatesum <- copy(isomir_readcount)  
templatesum[,(1:4):=NULL]
templatesum[,(2:4):=NULL]
templatesum[,(3):=NULL]
templatesum[,(4:13):=NULL]
setkey (templatesum, template5, template3, missmatch)
templatesum <- templatesum[,lapply(.SD,sum),by = key(templatesum)]
write.table(templatesum,"Isomir_Template_Overview.txt",row.names=F, sep = '\t')
meantemplatesum <- cbind.data.frame(as.data.frame (templatesum) [, grep(paste(groupnames, collapse="|"), colnames(templatesum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (templatesum) [, grep(x, colnames(templatesum))] )  ))
write.table(meantemplatesum,"Mean_Isomir_Template_Overview.txt",row.names=F, sep = '\t')
templatesum$template <- paste (templatesum$template5, templatesum$template3)
templatesum[,(1:2):=NULL]
meantemplatesum <- cbind.data.frame(as.data.frame (templatesum) [, grep(paste(groupnames, collapse="|"), colnames(templatesum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (templatesum) [, grep(x, colnames(templatesum))] )  ))
facet_names <- c('0mm' = "No Internal Missmatches", 'Poly' = "Polymorphic Isomirs")
ggplot(melt(meantemplatesum,id.vars=c("missmatch","template")),aes(variable,value,fill=template))+
  geom_bar(stat="identity") + facet_wrap(~missmatch,nrow=2, labeller = as_labeller(facet_names)) + labs(title="Hairpin Template Distribution", x = "Groups",   y= "Reads", fill = "5'/3' End mapping to Hairpin") + scale_fill_manual(values = c("#F781BF",  "#377EB8", "#4DAF4A", "#984EA3", "#A65628", "#FFFF33",  "#FF7F00", "#E41A1C"))
ggsave(filename="Template_Distribution.png", dpi = 600, width= 9, limitsize = F)

missmatchsum <- copy(isomir_readcount)  
missmatchsum[,(1:11):=NULL]
missmatchsum[,(2:10):=NULL]
setkey (missmatchsum, ntmm)
missmatchsum <- missmatchsum[,lapply(.SD,sum),by = key(missmatchsum)]
write.table(missmatchsum,"Isomir_Missmatch_Overview.txt",row.names=F, sep = '\t')
meanmissmatchsum <- cbind.data.frame(as.data.frame (missmatchsum) [, grep(paste(groupnames, collapse="|"), colnames(missmatchsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (missmatchsum) [, grep(x, colnames(missmatchsum))] )  ))
meanmissmatchsum <- setnames(meanmissmatchsum,1,"ntmm")
write.table(meanmissmatchsum,"Mean_Isomir_Missmatch_Overview.txt",row.names=F, sep = '\t')

ggplot(melt(meanmissmatchsum,id.vars=c("ntmm")),aes(ntmm,value,fill=variable))+
  geom_bar(position="dodge",stat="identity") + labs(title="Polymorphic Isomir Missmatch Distribution", x = "Number of internal Missmatches", y= "Reads", fill = "Groups")
ggsave(filename="Polymorphic_Missmatch_Distribution.png", dpi = 600)

modsum <- copy(isomir_readcount)
modsum[,(1):=NULL]
modsum[,(2:4):=NULL]
modsum[,(3:16):=NULL]
setkey (modsum,  mod5, mod3, ntlength)
modsum <- modsum[,lapply(.SD,sum),by = key(modsum)]    
write.table(modsum,"Modifications_Length_Sum_Overview.txt",row.names=F, sep = '\t')
meanmodsum <- cbind.data.frame(as.data.frame (modsum) [, grep(paste(groupnames, collapse="|"), colnames(modsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (modsum) [, grep(x, colnames(modsum))] )  ))
write.table(meanmodsum,"Mean_Modifications_Length_Sum_Overview.txt",row.names=F, sep = '\t')
modsum$mod <- paste (modsum$mod5, modsum$mod3)
modsum[,(1:2):=NULL]
setcolorder(modsum, c("ntlength", "mod", sample_names))
meanmodsum <- cbind.data.frame(as.data.frame (modsum) [, grep(paste(groupnames, collapse="|"), colnames(modsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (modsum) [, grep(x, colnames(modsum))] )  ))
ggplot(melt(meanmodsum,id.vars=c("mod","ntlength")), aes_string(x = "ntlength", y = "value", fill = "mod")) + geom_area(position = 'stack') + facet_wrap(~variable, nrow=3) + labs(title="Isomir Length Modification Distribution", x = "Sequence Length [nt]", y= "Reads", fill = "5'/3' End Modification") + scale_fill_manual(values = c("#999999", "#F781BF", "#984EA3",  "#4DAF4A", "#377EB8", "#E41A1C", "#A65628", "#FFFF33",  "#FF7F00"))
ggsave(filename="Isomir_Length_Modification_Distribution.png", dpi = 600, width = 9)
for (i in colnames(meanmodsum[-(1:2)])) {
  plot <- print(ggplot(meanmodsum, aes_string(x = "ntlength", y = as.name(i), fill = "mod")) + geom_area(position = 'stack') + labs(title= paste(i, " Isomir Length Modification Distribution",sep=""), x = "Sequence Length [nt]", y= "Reads", fill = "5'/3' End Modification") + scale_fill_manual(values = c("#999999", "#F781BF", "#984EA3",  "#4DAF4A", "#377EB8", "#E41A1C", "#A65628", "#FFFF33",  "#FF7F00")))
  ggsave(filename=paste(i,"_Isomir_Length_Modification_Distribution.png",sep=""), dpi = 600)
}
modsum <- copy(isomir_readcount)
modsum[,(1):=NULL]
modsum[,(2:4):=NULL]
modsum[,(3:17):=NULL]
setkey (modsum,  mod5, mod3)
modsum <- modsum[,lapply(.SD,sum),by = key(modsum)]    
write.table(modsum,"Modifications_Sum_Overview.txt",row.names=F, sep = '\t')
meanmodsum <- cbind.data.frame(as.data.frame (modsum) [, grep(paste(groupnames, collapse="|"), colnames(modsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (modsum) [, grep(x, colnames(modsum))] )  ))
write.table(meanmodsum,"Mean_Modifications_Sum_Overview.txt",row.names=F, sep = '\t')
modsum <- cbind.data.frame(mod = paste (modsum$mod5, modsum$mod3), modsum)
modsum[,(2:3):=NULL]
meanmodsum <- cbind.data.frame(as.data.frame (modsum) [, grep(paste(groupnames, collapse="|"), colnames(modsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (modsum) [, grep(x, colnames(modsum))] )  ))
meanmodsum <- setnames(meanmodsum,1,"mod")
ggplot(melt(meanmodsum,id="mod"), aes_string(x = "variable", y = "value", fill = "mod")) + geom_bar(position="dodge",stat="identity") + labs(title= "Isomir Modification Distribution",sep="", x = "Groups", y= "Reads", fill = "5'/3' End Modification") + scale_fill_manual(values = c("#999999", "#F781BF", "#984EA3",  "#4DAF4A", "#377EB8", "#E41A1C", "#A65628", "#FFFF33",  "#FF7F00"))
ggsave(filename="Modification_Distribution.png", dpi = 600, width= 9, limitsize = F)

endsum <- copy(isomir_readcount)
endsum[,(1:18):=NULL]
endsum[,(2:3):=NULL]
setkey (endsum,  end3seq)
endsum <- endsum[,lapply(.SD,sum),by = key(endsum)]    
meanendsum <- cbind.data.frame(as.data.frame (endsum) [, grep(paste(groupnames, collapse="|"), colnames(endsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (endsum) [, grep(x, colnames(endsum))] )  ))
meanmodsum <- setnames(meanmodsum,1,"end3seq")
write.table(meanendsum,"Mean_End_Sequence_Sum_Overview.txt",row.names=F, sep = '\t')

endsum <- copy(isomir_readcount)
endsum[,(1:5):=NULL]
endsum[,(2):=NULL]
endsum[,(3:12):=NULL]
endsum[,(4:5):=NULL]
endsum$seq3 <- ifelse (endsum$mod3=="Add3",substr(endsum$end3seq,7-nchar(as.character(endsum$seq3)),7-nchar(as.character(endsum$seq3))),substr(endsum$seq3,nchar(as.character(endsum$seq3)),nchar(as.character(endsum$seq3))))
endsum$end3seq <- substr(endsum$end3seq,7,7)
setkey (endsum, mod3, seq3, end3seq)
endsum <- endsum[,lapply(.SD,sum),by = key(endsum)] 
endsum <- endsum[seq3 != "N"]
endsum <- endsum[end3seq != "N"]
meanendsum <- cbind.data.frame(as.data.frame (endsum) [, grep(paste(groupnames, collapse="|"), colnames(endsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (endsum) [, grep(x, colnames(endsum))] )  ))
write.table(meanendsum,"Mean_End_Sequence_Modification_Sum_Overview.txt",row.names=F, sep = '\t')
facet_names <- c('A' = "Canonical miRNA ended with A", 'C' = "Canonical miRNA ended with C", 'G' = "Canonical miRNA ended with G", 'T' = "Canonical miRNA ended with T", 'X' = "No Modifications")
for (i in colnames(meanendsum[-(1:3)])) {
  plot <- print(ggplot(meanendsum, aes_string(x = "end3seq", y = as.name(i), fill = "mod3")) + geom_bar(position="dodge",stat="identity") + facet_wrap(~seq3, nrow=3, labeller = as_labeller(facet_names), scales = "free") + labs(title=paste(i," Sequence Changes at 3' End",sep=""), x = "Isomir ends with", y= "Reads", fill = "3' End Modification"))
  ggsave(filename=paste(i,"_End_Sequence_Modification_Distribution.png",sep=""), dpi = 600, width = 9)
}
endsum[,(1:2):=NULL]
setkey (endsum, end3seq)
endsum <- endsum[,lapply(.SD,sum),by = key(endsum)]
meanendsum <- cbind.data.frame(as.data.frame (endsum) [, grep(paste(groupnames, collapse="|"), colnames(endsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (endsum) [, grep(x, colnames(endsum))] )  ))
meanendsum <- setnames(meanendsum,1,"end3seq")
ggplot(melt(meanendsum,id="end3seq"), aes_string(x = "variable", y = "value", fill = "end3seq")) + geom_bar(position="dodge",stat="identity") + labs(title= "Ending Nucleotide Distribution",sep="", x = "Groups", y= "Reads", fill = "Ending Nucleotide")
ggsave(filename="End_Nucleotide_Distribution.png", dpi = 600, width= 9, limitsize = F)

seedsum <- copy(isomir_readcount)
seedsum[,(2:14):=NULL]
seedsum[,(3):=NULL]
seedsum[,(4:7):=NULL]
setkey (seedsum,  canonical, seedseq, m8seq)
seedsum <- seedsum[,lapply(.SD,sum),by = key(seedsum)]    
meanseedsum <- cbind.data.frame(as.data.frame (seedsum) [, grep(paste(groupnames, collapse="|"), colnames(seedsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (seedsum) [, grep(x, colnames(seedsum))] )  ))
write.table(meanseedsum,"Mean_Canonical_Seed_M8_Sum_Overview.txt",row.names=F, sep = '\t')
seedsum[,(3):=NULL]
setkey (seedsum,  canonical, seedseq)
seedsum <- seedsum[,lapply(.SD,sum),by = key(seedsum)]    
meanseedsum <- cbind.data.frame(as.data.frame (seedsum) [, grep(paste(groupnames, collapse="|"), colnames(seedsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (seedsum) [, grep(x, colnames(seedsum))] )  ))
write.table(meanseedsum,"Mean_Canonical_Seed_Sum_Overview.txt",row.names=F, sep = '\t')
seedsum[,(1):=NULL]
setkey (seedsum, seedseq)
seedsum <- seedsum[,lapply(.SD,sum),by = key(seedsum)]    
meanseedsum <- cbind.data.frame(as.data.frame (seedsum) [, grep(paste(groupnames, collapse="|"), colnames(seedsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (seedsum) [, grep(x, colnames(seedsum))] )  ))
meanseedsum <- setnames(meanseedsum,1,"seedseq")
write.table(meanseedsum,"Mean_Seed_Sum_Overview.txt",row.names=F, sep = '\t')


isomir_readcount_missmatch_position_split <- data.table(isomir_readcount, key=c("canonical","mod5","nt5","seq5","template5","mod3","nt3","seq3","template3","multiread","missmatch","ntmm","seed","seedseq","m8","m8seq","end3","end3seq","length","ntlength", sample_names))
isomir_readcount_missmatch_position_split <- isomir_readcount_missmatch_position_split[ntmm != "0"]
isomir_readcount_missmatch_position_split <- isomir_readcount_missmatch_position_split[, list(posmm = unlist(strsplit(posmm, ","))), by=key(isomir_readcount_missmatch_position_split)]
isomir_readcount_missmatch_position_split[, c("posmm","reference") := tstrsplit(posmm, ":", fixed=TRUE)]
isomir_readcount_missmatch_position_split[, c("reference","mmseq") := tstrsplit(reference, ">", fixed=TRUE)]
setcolorder(isomir_readcount_missmatch_position_split, c("canonical","mod5","nt5","seq5","template5","mod3","nt3","seq3","template3","multiread","missmatch","ntmm","posmm","reference","mmseq","seed","seedseq","m8","m8seq","end3","end3seq","length","ntlength", sample_names))
isomir_readcount_missmatch_position_split <- isomir_readcount_missmatch_position_split[mmseq != "N"]
isomir_readcount_missmatch_position_split$ntmm <- as.integer(isomir_readcount_missmatch_position_split$ntmm)
isomir_readcount_missmatch_position_split$posmm <- as.integer(isomir_readcount_missmatch_position_split$posmm)
isomir_readcount_missmatch_position_split <- isomir_readcount_missmatch_position_split[, posmm := posmm+1]

positionsum <- copy(isomir_readcount_missmatch_position_split) 
positionsum[,(1:11):=NULL]
positionsum[,(5:11):=NULL]
setkey (positionsum, ntmm,posmm,reference,mmseq,ntlength)
positionsum <- positionsum[,lapply(.SD,sum),by = key(positionsum)]
write.table(positionsum,"Polymorphic_Position_Modification_Sum_Overview.txt",row.names=F, sep = '\t')
positionsum[,(1):=NULL]
setkey (positionsum,  posmm,reference,mmseq,ntlength)
positionsum <- positionsum[,lapply(.SD,sum),by = key(positionsum)]
write.table(positionsum,"Polymorphic_Position_Length_Sum_Overview.txt",row.names=F, sep = '\t')
pos5endsum <- copy(positionsum)
pos5endsum[,(3:4):=NULL]
setkey (pos5endsum,  posmm,reference)
pos5endsum <- pos5endsum[,lapply(.SD,sum),by = key(pos5endsum)]
write.table(pos5endsum,"Polymorphic_Position_5end_Replaced_Reference_Overview.txt",row.names=F, sep = '\t')
pos5endsum <- copy(positionsum)
pos5endsum[,(2):=NULL]
pos5endsum[,(3):=NULL]
setkey (pos5endsum,  posmm,mmseq)
pos5endsum <- pos5endsum[,lapply(.SD,sum),by = key(pos5endsum)]
write.table(pos5endsum,"Polymorphic_Position_5end_Replacing_Mismatch_Overview.txt",row.names=F, sep = '\t')
pos5endsum <- copy(isomir_readcount_missmatch_position_split)
pos5endsum[,(1):=NULL]
pos5endsum[,(2:11):=NULL]
pos5endsum[,(3:12):=NULL]
setkey (pos5endsum,  mod5, posmm)
pos5endsum <- pos5endsum[,lapply(.SD,sum),by = key(pos5endsum)]
write.table(pos5endsum,"Polymorphic_Position_5end_Sum_Overview.txt",row.names=F, sep = '\t')
meanpos5endsum <- cbind.data.frame(as.data.frame (pos5endsum) [, grep(paste(groupnames, collapse="|"), colnames(pos5endsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (pos5endsum) [, grep(x, colnames(pos5endsum))] )  ))
meanpos5endsum <- data.table(meanpos5endsum)
ggplot(melt(meanpos5endsum[posmm < 11],id=c("posmm","mod5")), aes_string(x = "posmm", y = "value", fill = "variable")) + geom_bar(position="dodge",stat="identity") + labs(title= "Mismatch 5'end Position Distribution",sep="", x = "Nukleotide Position", y= "Reads", fill = "Groups") + facet_wrap(~mod5,nrow=3, scales = 'free') + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mismatch 5'end Position Distribution.png", dpi = 600, width = 9, height = 6)


positionsum <- positionsum[, pos3mm := (ntlength-posmm+1)*(-1)]
setcolorder(positionsum, c("posmm", "pos3mm", "reference", "mmseq", "ntlength", sample_names))
pos3endsum <- copy(positionsum)
pos3endsum[,(1):=NULL]
pos3endsum[,(3:4):=NULL]
setkey (pos3endsum,  pos3mm, reference)
pos3endsum <- pos3endsum[,lapply(.SD,sum),by = key(pos3endsum)]
write.table(pos3endsum,"Polymorphic_Position_3end_Replaced_Reference_Overview.txt",row.names=F, sep = '\t')
pos3endsum <- copy(positionsum)
pos3endsum[,(1):=NULL]
pos3endsum[,(2):=NULL]
pos3endsum[,(3):=NULL]
setkey (pos3endsum,  pos3mm, mmseq)
pos3endsum <- pos3endsum[,lapply(.SD,sum),by = key(pos3endsum)]
write.table(pos3endsum,"Polymorphic_Position_3end_Replacing_MM_Overview.txt",row.names=F, sep = '\t')
pos3endsum <- copy(isomir_readcount_missmatch_position_split)
pos3endsum <- pos3endsum[, pos3mm := (ntlength-posmm+1)*(-1)]
pos3endsum[,(1:5):=NULL]
pos3endsum[,(2:18):=NULL]
setcolorder(pos3endsum, c("mod3", "pos3mm", sample_names))
setkey (pos3endsum,  mod3, pos3mm)
pos3endsum <- pos3endsum[,lapply(.SD,sum),by = key(pos3endsum)]
write.table(pos3endsum,"Polymorphic_Position_3end_Sum_Overview.txt",row.names=F, sep = '\t')
meanpos3endsum <- cbind.data.frame(as.data.frame (pos3endsum) [, grep(paste(groupnames, collapse="|"), colnames(pos3endsum), invert=T)],sapply(groupnames, function(x) rowMeans(as.data.frame (pos3endsum) [, grep(x, colnames(pos3endsum))] )  ))
meanpos3endsum <- data.table(meanpos3endsum)
ggplot(melt(meanpos3endsum[pos3mm > -11],id=c("pos3mm","mod3")), aes_string(x = "pos3mm", y = "value", fill = "variable")) + geom_bar(position="dodge",stat="identity") + labs(title= "Mismatch 3'end Position Distribution",sep="", x = "Nukleotide Position", y= "Reads", fill = "Groups") + facet_wrap(~mod3,nrow=3, scales = 'free') + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10)) + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10))
ggsave(filename="Mismatch 3'end Position Distribution.png", dpi = 600, width = 9, height = 6)

meanpos3endsum <- setnames(meanpos3endsum,1:2,c("mod", "posmm"))
meanpos5endsum <- setnames(meanpos5endsum,1:2,c("mod", "posmm"))
meanposendsum <- rbind(meanpos5endsum, meanpos3endsum)
meanposendsum$mod <- sub("[0-9]+","",meanposendsum$mod)
meanposendsum <- meanposendsum[posmm > -11 & posmm < 11]
newrow <- data.frame(c ("Add", "Can", "Trim"), "<- 5|3 ->", matrix(c(rep.int(0,length(meanposendsum)-2)),nrow=3,ncol=length(meanposendsum)-2))
colnames(newrow) <- colnames(meanposendsum)
meanposendsum <- rbind( meanposendsum, newrow)
meanposendsum <- melt(meanposendsum,id=c("posmm","mod"))
meanposendsum$posmm <- factor(meanposendsum$posmm, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "<- 5|3 ->", "-10", "-9", "-8", "-7", "-6", "-5", "-4", "-3", "-2", "-1"))
ggplot(meanposendsum, aes_string(x = "posmm", y = "value", fill = "variable")) + geom_bar(position="dodge",stat="identity") + labs(title= "Mismatch 5'/3'end Position Distribution",sep="", x = "Nukleotide Position", y= "Reads", fill = "Groups") + facet_wrap(~mod,nrow=3, scales ="free_y") + scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10)) + geom_vline(xintercept = 11, size= 1)
ggsave(filename="Mismatch Position Distribution.png", dpi = 600, width = 9, height = 6)

