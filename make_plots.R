
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5){
	('Please provide exactly five arguments, in this order: [input_coord_file] [subject_contig_data] [query_contig_data] [highlight_data] [output_file]')
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
	# order the df from largest to smallest contig
	ocdf = contig_data_frame[order(contig_data_frame$contig_length,decreasing=TRUE),]
	# make a name_key, in this format: C{contig_number}_{contig_length}
	ocdf$name_key = sapply(1:nrow(ocdf),function(x){paste('C',x,'_',ocdf$contig_length[x],sep='')})
	return(ocdf)
	#contig_lengths = sort(contig_data_frame$contig_length,decreasing=TRUE)
	#name_key = sapply(1:length(contig_lengths),function(x){paste('C',x,'_',contig_lengths[x],sep='')})
	#names(name_key) = as.character(contig_lengths)
	#return(name_key)
}

subjectContig = read.csv(subject_contig_data, sep = '\t',header = T)
# order the contig data and add an extra column
subjectContigDF = make_contig_data(subjectContig)
# use this to make a key to convert contig names from mummer's names to mine
subjectContigData = subjectContigDF$name_key
names(subjectContigData) = subjectContigDF$contig_name
#subjectContigData = make_contig_data(subjectContig)

queryContig = read.csv(query_contig_data, sep = '\t',header = T)
queryContigDF = make_contig_data(queryContig)
# use this to make a key to convert contig names from mummer's names to mine
queryContigData = queryContigDF$name_key
names(queryContigData) = queryContigDF$contig_name
#queryContigData = make_contig_data(queryContig)

#print(subjectContigDF)
#print(queryContigDF)

# This section can subset the plots to specific contigs 
# This should eventually be added as an option, but for now it needs to be manual
# plotdata = plotdata[(plotdata$subjectScafLen > 450000 & plotdata$queryScafLen > 450000) | (plotdata$subjectTag == 'contig_74'),]
# plotdata = plotdata[(plotdata$subjectTag == 'contig_42') & (plotdata$queryTag == 'scaffold_6'),]
# subjectContigData = subjectContigData[names(subjectContigData) %in% as.character(plotdata$subjectScafLen)]
# queryContigData = queryContigData[names(queryContigData) %in% as.character(plotdata$queryScafLen)]

# use this information to rename the contigs in the coord file
# rename with contig lengths
#plotdata$subjectChr = factor(subjectContigData[as.character(plotdata$subjectScafLen)])
#plotdata$queryChr = factor(queryContigData[as.character(plotdata$queryScafLen)])
# rename with contig tags
plotdata$subjectChr = factor(subjectContigData[plotdata$subjectTag])
plotdata$queryChr = factor(queryContigData[plotdata$queryTag])

#print(head(plotdata))
#print(tail(plotdata))

make_limitdata = function(subjectContigDF,queryContigDF){
	limitdata = data.frame('subjectChr' = NA, 'subjectChrSize' = NA,'queryChr' = NA, 'queryChrSize' = NA)
	r = 1
	for (i in 1:nrow(subjectContigDF)){
		for (j in 1:nrow(queryContigDF)){
			limitdata[r,] = c(subjectContigDF$name_key[i],subjectContigDF$contig_length[i],queryContigDF$name_key[j],queryContigDF$contig_length[j])
			r = r + 1
		}
	}
	limitdata$subjectChr = factor(limitdata$subjectChr)
	limitdata$queryChr = factor(limitdata$queryChr)
    limitdata$subjectChrSize = as.numeric(limitdata$subjectChrSize)
    limitdata$queryChrSize = as.numeric(limitdata$queryChrSize)
	return(limitdata)
}
limitdata = make_limitdata(subjectContigDF,queryContigDF)

#print(limitdata)


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
	geom_point(data=limitdata,aes(x=subjectChrSize,y=queryChrSize),color='white',size=0.01) +
	geom_segment(data=plotdata, aes(x=subjectStart, xend=subjectEnd, y=queryStart, yend=queryEnd)) +  
	theme(strip.text.x = element_text(size = 6),strip.text.y = element_text(size = 6),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1)) +
	xlab(xlabel2) + ylab(ylabel2)
#dev.off()



# generate a key that converts between the original and the new contig names
# the new contig name (such as C6_972738) is the value, and the original contig name (such as contig_6) is used as the vector name
# subjectContigKey = plotdata$subjectChr
# names(subjectContigKey) = plotdata$subjectTag
# subjectContigKey = subjectContigKey[!duplicated(subjectContigKey)]
# queryContigKey = plotdata$queryChr
# names(queryContigKey) = plotdata$queryTag
# queryContigKey = queryContigKey[!duplicated(queryContigKey)]

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
		print(c(htype,hname,subject_contig_data,hcontig,names(subjectContigData)))
		if ((htype=='subject') & (grepl(hname,subject_contig_data)) & (hcontig %in% names(subjectContigData))){
			print('RUNNING')
			hcontig2 = subjectContigData[hcontig]
			print(hcontig2)
			linedata2[linedata2$subjectChr == hcontig2,]$xstart = hstart
			linedata2[linedata2$subjectChr == hcontig2,]$xend = hend
			print(linedata2)
			gplot1 = gplot1 + geom_segment(data=linedata2,aes(x=xstart,xend=xend,y=ystart,yend=yend,color='red',linewidth=5),show.legend = F)
		}
		if ((htype=='query') & (grepl(hname,query_contig_data)) & (hcontig %in% names(queryContigData))){
			hcontig2 = queryContigData[hcontig]
			linedata2[linedata2$queryChr == hcontig2,]$ystart = hstart
			linedata2[linedata2$queryChr == hcontig2,]$yend = hend
			gplot1 = gplot1 + geom_segment(data=linedata2,aes(x=xstart,xend=xend,y=ystart,yend=yend,color='red',linewidth=5),show.legend = F)
		}
	}
}

pdf(output_file)
gplot1
dev.off()



