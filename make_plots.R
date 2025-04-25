
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5){
	print('Please provide exactly five arguments, in this order: [input_coord_file] [subject_contig_data] [query_contig_data] [highlight_data] [output_file]')
} else {
	input_file = args[1]
	subject_contig_data = args[2]
	query_contig_data = args[3]
	highlight_data = args[4]
	output_file = args[5]
}

# read in coordinate data
plotdata = read.csv(input_file, sep = '\t',header = F,skip=4)
colnames(plotdata) = c('subjectStart','subjectEnd','queryStart','queryEnd','subjectSeqLen','querySeqLen','percentIdent','subjectScafLen','queryScafLen','subjectCov','queryCov','subjectTag','queryTag')


# read in contig data and assign chromosome names, from largest to smallest
make_contig_data = function(contig_data_frame){
	contig_lengths = sort(contig_data_frame$contig_length,decreasing=TRUE)
	name_key = sapply(1:length(contig_lengths),function(x){paste('C',x,'_',contig_lengths[x],sep='')})
	names(name_key) = as.character(contig_lengths)
	return(name_key)
}

subjectContig = read.csv(subject_contig_data, sep = '\t',header = T)
subjectContigData = make_contig_data(subjectContig)
queryContig = read.csv(query_contig_data, sep = '\t',header = T)
queryContigData = make_contig_data(queryContig)

# This section can subset the plots to specific contigs 
# This should eventually be added as an option, but for now it needs to be manual
#plotdata = plotdata[(plotdata$subjectScafLen > 450000 & plotdata$queryScafLen > 450000) | (plotdata$subjectTag == 'contig_74'),]
plotdata = plotdata[(plotdata$subjectTag == 'contig_42') & (plotdata$queryTag == 'scaffold_6'),]
subjectContigData = subjectContigData[names(subjectContigData) %in% as.character(plotdata$subjectScafLen)]
queryContigData = queryContigData[names(queryContigData) %in% as.character(plotdata$queryScafLen)]

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


# make a copy of limitdata, except for red lines to highlight the transposon
# linedata = limitdata
# linedata$xstart = 0
# linedata$xend = 0
# linedata$ystart = 0
# linedata$yend = 0
# linedata2 = linedata
# if (grepl('GCA_000597905.1_ASM59790v1',subject_contig_data)){
# 	linedata[linedata$subjectChr == 'C3_25284',]$xstart = 900
# 	linedata[linedata$subjectChr == 'C3_25284',]$xend = 10905
# }
# if (grepl('MI_KPC_689_flye_medaka_polypolish',subject_contig_data)){
# 	linedata[linedata$subjectChr == 'C2_390225',]$xstart = 39525
# 	linedata[linedata$subjectChr == 'C2_390225',]$xend = 49530
# }
# if (grepl('MI_KPC_59_flye_medaka_polypolish',subject_contig_data)){
# 	linedata[linedata$subjectChr == 'C4_43621',]$xstart = 36394
# 	linedata[linedata$subjectChr == 'C4_43621',]$xend = 26389
# }
# if (grepl('MI_KPC_112_flye_medaka_polypolish',subject_contig_data)){
# 	linedata[linedata$subjectChr == 'C4_43620',]$xstart = 36931
# 	linedata[linedata$subjectChr == 'C4_43620',]$xend = 26926
# }
# if (grepl('ML2TQ8',subject_contig_data)){
# 	linedata[linedata$subjectChr == 'C39_46260',]$xstart = 13193
# 	linedata[linedata$subjectChr == 'C39_46260',]$xend = 7806
# }
# if (grepl('Chi_Caur_5',query_contig_data)){
# 	linedata2[linedata2$queryChr == 'C6_945903',]$ystart = 941142
# 	linedata2[linedata2$queryChr == 'C6_945903',]$yend = 944236
# }
# #print(linedata)

# get axis names
xlabel1 = strsplit(subject_contig_data,'_contig_data')[[1]][1]
xlabel2 = strsplit(xlabel1,'/')[[1]]
xlabel2 = xlabel2[length(xlabel2)]
ylabel1 = strsplit(query_contig_data,'_contig_data')[[1]][1]
ylabel2 = strsplit(ylabel1,'/')[[1]]
ylabel2 = ylabel2[length(ylabel2)]


# generate plot
#pdf(output_file)
# ggplot() + facet_grid(queryChr~subjectChr,scales = 'free',drop = F) + xlim(0,NA) + ylim(0,NA)+ coord_cartesian(clip = "off")  + ggplot2::theme_bw() +
# 	geom_segment(data=plotdata, aes(x=subjectStart, xend=subjectEnd, y=queryStart, yend=queryEnd)) + 
# 	geom_point(data=limitdata,aes(x=subjectChrSize,y=queryChrSize),color='white',size=0.01) + theme(strip.text.x = element_text(size = 6),strip.text.y = element_text(size = 6),axis.text.x = element_blank()) +
# 	xlab(xlabel2) + ylab(ylabel2)
gplot1 = ggplot() + facet_grid(queryChr~subjectChr,scales = 'free',drop = F) + xlim(0,NA) + ylim(0,NA)+ coord_cartesian(clip = "off")  + ggplot2::theme_bw() +
	geom_segment(data=plotdata, aes(x=subjectStart, xend=subjectEnd, y=queryStart, yend=queryEnd)) + 
	geom_point(data=limitdata,aes(x=subjectChrSize,y=queryChrSize),color='white',size=0.01) + 
	theme(strip.text.x = element_text(size = 6),strip.text.y = element_text(size = 6),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1)) +
	xlab(xlabel2) + ylab(ylabel2)
#dev.off()



# generate a key that converts between the original and the new contig names
# the new contig name (such as C6_972738) is the value, and the original contig name (such as contig_6) is used as the vector name
subjectContigKey = plotdata$subjectChr
names(subjectContigKey) = plotdata$subjectTag
subjectContigKey = subjectContigKey[!duplicated(subjectContigKey)]
print(subjectContigKey)
queryContigKey = plotdata$queryChr
names(queryContigKey) = plotdata$queryTag
queryContigKey = queryContigKey[!duplicated(queryContigKey)]

# make a default data frame that won't draw any lines
linedata = limitdata
linedata$xstart = 0
linedata$xend = 0
linedata$ystart = 0
linedata$yend = 0

if (highlight_data != 'NA'){
	# read in the highlight data
	highlightdata = read.csv(highlight_data, sep = '\t',header = T)
	# for each row/column that needs to be highlighted, add it to the plot
	for (i in 1:nrow(highlightdata)){
		linedata2 = linedata
		htype = highlightdata$type[i]
		hname = highlightdata$name[i]
		hcontig = highlightdata$contig[i]
		hstart = highlightdata$start[i]
		hend = highlightdata$end[i]
		print(c(htype,hname,subject_contig_data,hcontig,names(subjectContigKey)))
		if ((htype=='subject') & (grepl(hname,subject_contig_data)) & (hcontig %in% names(subjectContigKey))){
			hcontig2 = subjectContigKey[hcontig]
			linedata2[linedata2$subjectChr == hcontig2,]$xstart = hstart
			linedata2[linedata2$subjectChr == hcontig2,]$xend = hend
			print(linedata2)
			gplot1 = gplot1 + geom_segment(data=linedata2,aes(x=xstart,xend=xend,y=ystart,yend=yend,color='red',linewidth=5),show.legend = F)
		}
		if ((htype=='query') & (grepl(hname,query_contig_data)) & (hcontig %in% names(queryContigKey))){
			hcontig2 = queryContigKey[hcontig]
			linedata2[linedata2$queryChr == hcontig2,]$ystart = hstart
			linedata2[linedata2$queryChr == hcontig2,]$yend = hend
			gplot1 = gplot1 + geom_segment(data=linedata2,aes(x=xstart,xend=xend,y=ystart,yend=yend,color='red',linewidth=5),show.legend = F)
		}
	}
}

pdf(output_file)
gplot1
dev.off()



