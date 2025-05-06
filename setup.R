
# this will install ggplot if you don't have it already
if (!('ggplot2' %in% installed.packages()[,'Package'])){
	install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
