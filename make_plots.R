
# this will install ggplot if you don't have it already
if (!('ggplot2' %in% installed.packages()[,'Package'])){
	install.packages("ggplot2", repos = "https://cloud.r-project.org")
}

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
	# read in the data frame
	df = read.csv(contig_data_frame, sep = '\t',header = T)
	# order the df from largest to smallest contig
	ocdf = df[order(df$contig_length,decreasing=TRUE),]
	# make a name_key as a column in the df, in this format: C{contig_number}_{contig_length}
	ocdf$name_key = sapply(1:nrow(ocdf),function(x){paste('C',x,'_',ocdf$contig_length[x],sep='')})
	# generate a vector that converts mummer's contig names to these new names
	cdata_vec = ocdf$name_key
	names(cdata_vec) = ocdf$contig_name
	# return both
	return(list('DF' = ocdf,'contig_key' = cdata_vec))
}


subjectContigOut = make_contig_data(subject_contig_data)
subjectContigDF = subjectContigOut$DF
subjectContigData = subjectContigOut$contig_key

queryContigOut = make_contig_data(query_contig_data)
queryContigDF = queryContigOut$DF
queryContigData = queryContigOut$contig_key


# use this information to rename the contigs in the coord file 
# the factor levels must be the chromosome names in ascending order for the subject and descending order for the query
plotdata$subjectChr = factor(subjectContigData[plotdata$subjectTag], levels = subjectContigDF$name_key)
plotdata$queryChr = factor(queryContigData[plotdata$queryTag], levels = rev(queryContigDF$name_key))


make_limitdata = function(subjectContigDF,queryContigDF){
	limitdata = data.frame('subjectChr' = NA, 'subjectChrSize' = NA,'queryChr' = NA, 'queryChrSize' = NA)
	r = 1
	for (i in 1:nrow(subjectContigDF)){
		for (j in 1:nrow(queryContigDF)){
			limitdata[r,] = c(subjectContigDF$name_key[i],subjectContigDF$contig_length[i],queryContigDF$name_key[j],queryContigDF$contig_length[j])
			r = r + 1
		}
	}
	# use the same levels as above
	limitdata$subjectChr = factor(limitdata$subjectChr, levels = subjectContigDF$name_key)
	limitdata$queryChr = factor(limitdata$queryChr, levels = rev(queryContigDF$name_key))
    limitdata$subjectChrSize = as.numeric(limitdata$subjectChrSize)
    limitdata$queryChrSize = as.numeric(limitdata$queryChrSize)
	return(limitdata)
}
limitdata = make_limitdata(subjectContigDF,queryContigDF)


# get axis names
xlabel1 = strsplit(subject_contig_data,'_contig_data')[[1]][1]
xlabel2 = strsplit(xlabel1,'/')[[1]]
xlabel2 = xlabel2[length(xlabel2)]
ylabel1 = strsplit(query_contig_data,'_contig_data')[[1]][1]
ylabel2 = strsplit(ylabel1,'/')[[1]]
ylabel2 = ylabel2[length(ylabel2)]


# generate plot
gplot1 = ggplot() + facet_grid(queryChr~subjectChr,scales = 'free',drop = F) + xlim(0,NA) + ylim(0,NA)+ coord_cartesian(clip = "off")  + ggplot2::theme_bw() +
	geom_point(data=limitdata,aes(x=subjectChrSize,y=queryChrSize),color='white',size=0.01) +
	geom_segment(data=plotdata, aes(x=subjectStart, xend=subjectEnd, y=queryStart, yend=queryEnd)) +  
	theme(strip.text.x = element_text(size = 6),strip.text.y = element_text(size = 6),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1)) +
	xlab(xlabel2) + ylab(ylabel2)


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
		if ((htype=='subject') & (grepl(hname,subject_contig_data)) & (hcontig %in% names(subjectContigData))){
			hcontig2 = subjectContigData[hcontig]
			linedata2[linedata2$subjectChr == hcontig2,]$xstart = hstart
			linedata2[linedata2$subjectChr == hcontig2,]$xend = hend
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
nothing = dev.off()



