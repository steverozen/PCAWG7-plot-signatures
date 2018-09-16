#!/usr/bin/Rscript





## input: indel.catalog.csv 
args = commandArgs(trailingOnly=TRUE)
inputfile = args[1]

library(RColorBrewer)
dinuc.class.col = brewer.pal(10, "Paired")

## plot function
plot_dinus <- function(counts, display.name, mut.type=NULL, classes=TRUE, grid.line=TRUE){
	num.classes = length(counts)
	cols = rep(dinuc.class.col, c(9,6,9,6,9,6,6,9,9,9))
	## get ylim
	ymax = ifelse(max(counts) * 1.3 > 1, 1, max(counts) * 1.3)
	# barplot
	bp = barplot(counts, main='', xaxt='n', ylim=c(0,ymax), xlab=NA, yaxt='n', ylab=NA, 
	             lwd=if(num.classes<200)3 else 1, space=1.35, border=NA, col=cols, xpd=NA)
	if (grid.line == TRUE){
	  # draw box and grid lines
	  segments(bp[1]-0.5, seq(0,ymax,ymax/4), bp[num.classes]+0.5, seq(0,ymax,ymax/4), col='grey90')
	  rect(xleft=bp[1]-1, 0, xright=bp[num.classes]+1, ymax, col=NA, border='grey60', lwd = 1.5)
	  # draw y axis ticks
	  segments(bp[1]-1, seq(0,ymax,ymax/4), bp[1]+1, seq(0,ymax,ymax/4), col='grey60')
	  segments(bp[num.classes]-1, seq(0,ymax,ymax/4), bp[num.classes]+1, seq(0,ymax,ymax/4), col='grey60')
	  # draw y axis
	  y.axis.values = seq(0,ymax,ymax/4)
	  y.axis.labels = format(round(y.axis.values*100, 1), nsmall = 1)
	  text(0.35, y.axis.values, labels=paste(y.axis.labels, '%', sep=''), 
	        las=1, adj=1, xpd=NA, family='Courier', cex=0.75)
	} else {
	  rect(xleft=bp[1]-1, 0, xright=bp[num.classes]+1, ymax, col=NA, border='grey60', lwd = 0.5)
	  # draw y axis ticks
	  segments(bp[1]-1, seq(0,ymax,ymax/2), bp[1]+1, seq(0,ymax,ymax/2), col='grey60', lwd = 0.5)
	  segments(bp[num.classes]-1, seq(0,ymax,ymax/2), bp[num.classes]+1, seq(0,ymax,ymax/2), col='grey60', lwd = 0.5)
	  # draw y axis
	  y.axis.values = seq(0,ymax,ymax/2)
	  y.axis.labels = format(round(y.axis.values*100, 1), nsmall = 1)
	  # text(0.35, y.axis.values, labels=paste(y.axis.labels, '%', sep=''), 
	  #      las=1, adj=1, xpd=NA, family='Courier', cex=0.75)
	}
	
	maj.class.names = paste(c('AC', 'AT', 'CC', 'CG', 'CT', 'GC', 'TA', 'TC', 'TG', 'TT'), 'NN', sep='>')
	#maj.class.names = c('AC', 'AT', 'CC', 'CG', 'CT', 'GC', 'TA', 'TC', 'TG', 'TT')
	if (classes == TRUE){
	  # draw lines above each class:
	  x.left = bp[c(0,9,15,24,30,39,45,51,60,69)+1]-0.5
	  x.right = bp[c(9,15,24,30,39,45,51,60,69,78)]+0.5
	  rect(xleft=x.left, ymax*1.01, xright=x.right, ymax*1.08, col=dinuc.class.col, border=NA, xpd=NA)
	  # mutation class labels at the top of the figure:
	  text((x.left+x.right)/2, ymax*1.123, labels=maj.class.names, cex=1, xpd=NA)
	} else if (classes == 'top'){
    # draw lines above each class:
    x.left = bp[c(0,9,15,24,30,39,45,51,60,69)+1]-0.5
    x.right = bp[c(9,15,24,30,39,45,51,60,69,78)]+0.5
    rect(xleft=x.left, ymax*1.01, xright=x.right, ymax*1.1, col=dinuc.class.col, border=NA, xpd=NA)
    ## mutation class labels at the top of the figure:
    #text((x.left+x.right)/2, ymax*1.2, labels=maj.class.names, cex=0.5, xpd=NA)
	} else if (classes == 'bottom'){
	  # draw lines above each class:
	  x.left = bp[c(0,9,15,24,30,39,45,51,60,69)+1]-0.5
	  x.right = bp[c(9,15,24,30,39,45,51,60,69,78)]+0.5
	  rect(xleft=x.left, ymax*(-0.1), xright=x.right, ymax*(-0.01), col=dinuc.class.col, border=NA, xpd=NA)
	  # # mutation class labels at the top of the figure:
	  # text((x.left+x.right)/2, ymax*(-0.2), labels=maj.class.names, cex=0.5, xpd=NA)
	} 
	
  # x axis label
  if (!is.null(mut.type)){
    text(bp, -ymax/100, labels=mut.type, cex=0.75, srt=90, adj=1, family='Courier', col='grey60', xpd=NA)
  }
  # sig name
	display.name = gsub('Signature\\.Dinucs\\.0', 'DBS', display.name)
	display.name = gsub('Signature\\.Dinucs\\.', 'DBS', display.name)
	if (grid.line == TRUE){
    text(1.5, ymax*7/8, labels=display.name, adj=0, cex=1.5, font=2)
  } else{
    text(1.5, ymax*7/8, labels=display.name, adj=0, cex=1, font=2)
  }
}

plot_dinucs_file <- function(inputfile, outputfile = NULL){
	data <- read.csv(inputfile)
	names.data <- read.csv(inputfile, stringsAsFactors = F, header = F)
	names.data <- t(names.data[1,])
	names(data) <- names.data
	
	na.sigs <- names(which(is.na(colSums(data[,2:ncol(data)]))))
	data <- data[,!names(data)%in%na.sigs]
	outputfile1 <- gsub('csv', 'pdf', inputfile)
	outputfile2 <- gsub('csv', 'combined.pdf', inputfile)
	pdf(outputfile1, width=12, height=2.75, onefile=T, useDingbats = F)
	par(mar=c(1.8,0.5,2.2,0), oma=c(0,0,0,0))
	for (i in 2:ncol(data)){
	  mut.type = gsub('.*>','',data[,1])
	  plot_dinus(data[,i], colnames(data)[i], mut.type=mut.type)
	}
	invisible(dev.off())
	
	## generate seperate png files
	outputfolder <- gsub('\\.csv', '', inputfile)
	for (i in 2:ncol(data)){
	  png(paste(outputfolder, "/", colnames(data)[i], ".png", sep=""), width=12, height=2.75, units = 'in', res = 300)
	  par(mar=c(1.8,0.5,2.2,0), oma=c(0,0,0,0))
	  mut.type = gsub('.*>','',data[,1])
	  plot_dinus(data[,i], colnames(data)[i], mut.type=mut.type)
	  dev.off()
	}
	
	cairo_pdf(outputfile2, width=11.6929, height=2.2977, onefile=T)
	par(mfrow=c(3,4), mar=c(0.3,0.1,0.3,0.1), oma=c(2,2,2,0))
  for (i in 2:5){
	  plot_dinus(data[,i], colnames(data)[i], classes='top', grid.line=F)
  }
	for (i in 6:8){
	  plot_dinus(data[,i], colnames(data)[i], classes=F, grid.line=F)
	}
	for (i in 9:12){
	  plot_dinus(data[,i], colnames(data)[i], classes='bottom', grid.line=F)
	}
	invisible(dev.off())
}

## plot the indel catalogs
plot_dinucs_file(inputfile)

### for Jaegil's sigs combined plot
if (FALSE){
  inputfile = '~/Desktop/ICGC-upload/Mutation_Signatures/approach_2/approach_2_DBS.signature.20180420.csv'
  data <- read.csv(inputfile)
  outputfile2 <- gsub('csv', 'combined.pdf', inputfile)
  cairo_pdf(outputfile2, width=11.6929, height=2.9277, onefile=T)
  par(mfrow=c(4,4), mar=c(0.3,0.1,0.3,0.1), oma=c(2,2,2,0))
  for (i in 2:5){
    plot_dinus(data[,i], colnames(data)[i], classes='top', grid.line=F)
  }
  for (i in 6:12){
    plot_dinus(data[,i], colnames(data)[i], classes=F, grid.line=F)
  }
  for (i in 13:16){
    plot_dinus(data[,i], colnames(data)[i], classes='bottom', grid.line=F)
  }
  invisible(dev.off())
}