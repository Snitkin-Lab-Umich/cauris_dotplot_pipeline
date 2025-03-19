
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4){
	print('Please provide exactly four arguments, in this order: [input_coord_file] [subject_contig_data] [query_contig_data] [output_file]')
} else {
	input_file = args[1]
	subject_contig_data = args[2]
	query_contig_data = args[3]
	output_file = args[4]
}

# read in coordinate data
plotdata = read.csv(args[1], sep = '\t',header = F,skip=4)
colnames(plotdata) = c('subjectStart','subjectEnd','queryStart','queryEnd','subjectSeqLen','querySeqLen','percentIdent','subjectScafLen','queryScafLen','subjectCov','queryCov','subjectTag','queryTag')

# read in contig data and assign chromosome names, from largest to smallest
make_contig_data = function(contig_data_frame){
	contig_lengths = contig_data_frame$contig_length
	name_key = sapply(1:length(contig_lengths),function(x){paste('chromosome_',x,sep='')})
	names(name_key) = as.character(sort(contig_lengths,decreasing=TRUE))
	return(name_key)
}

subjectContig = read.csv(subject_contig_data, sep = '\t',header = T)
subjectContigData = make_contig_data(subjectContig)
queryContig = read.csv(query_contig_data, sep = '\t',header = T)
queryContigData = make_contig_data(queryContig)

# use this information to rename the contigs in the coord file
plotdata$subjectChr = factor(subjectContigData[as.character(plotdata$subjectScafLen)])
plotdata$queryChr = factor(queryContigData[as.character(plotdata$queryScafLen)])

# create an additional data frame with only one entry for each unique combination
make_limitdata = function(subjectContigData,queryContigData){
	limitdata = data.frame('subjectChr' = NA, 'subjectChrSize' = NA,'queryChr' = NA, 'queryChrSize' = NA)
	r = 1
	for (i in 1:length(subjectContigData)){
		for (j in 1:length(queryContigData)){
			limitdata[r,] = c(subjectContigData[i],as.integer(names(subjectContigData)[i]),queryContigData[j],as.integer(names(queryContigData)[j]))
			r = r + 1
		}
	}
	limitdata$subjectChr = factor(limitdata$subjectChr)
	limitdata$queryChr = factor(limitdata$queryChr)
    limitdata$subjectChrSize = as.numeric(limitdata$subjectChrSize)
    limitdata$queryChrSize = as.numeric(limitdata$queryChrSize)
	return(limitdata)
}
limitdata = make_limitdata(subjectContigData,queryContigData)

# get axis names
xlabel1 = strsplit(subject_contig_data,'_contig_data')[[1]][1]
xlabel2 = strsplit(xlabel1,'/')[[1]]
xlabel2 = xlabel2[length(xlabel2)]
ylabel1 = strsplit(query_contig_data,'_contig_data')[[1]][1]
ylabel2 = strsplit(ylabel1,'/')[[1]]
ylabel2 = ylabel2[length(ylabel2)]


pdf(output_file)
ggplot() + facet_grid(queryChr~subjectChr,scales = 'free',drop = F) + xlim(0,NA) + ylim(0,NA)+ coord_cartesian(clip = "off")  + ggplot2::theme_bw() +
	geom_segment(data=plotdata, aes(x=subjectStart, xend=subjectEnd, y=queryStart, yend=queryEnd)) + 
	geom_point(data=limitdata,aes(x=subjectChrSize,y=queryChrSize),color='white',size=0.01) + theme(strip.text.x = element_text(size = 6),strip.text.y = element_text(size = 6),axis.text.x = element_blank()) +
	xlab(xlabel2) + ylab(ylabel2)
dev.off()

