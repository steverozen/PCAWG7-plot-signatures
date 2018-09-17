###########################################################################
###########################################################################
# Code to plot unstranded SBS (single base substitution) spectra and signatures
#
# 2018 09 17
# 
# v0.2
# 
# An alpha version
# 
# Copyright (C) 2018 Steven G. Rozen and Mi Ni Huang
# 
# The code is released under GPL-3
# https://www.gnu.org/licenses/gpl-3.0.en.html
# you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option)
# any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# Contact: steverozen@gmail.com
###########################################################################
###########################################################################


class.col  <- c('#03bcee',
                '#010101',
                '#e32926',
                '#999999',
                '#a1ce63',
                '#ebc6c4')

## plot function
plot_snvs <- function(counts, display.name, mut.type=NULL, classes=TRUE, grid.line=TRUE){
	num.classes = length(counts)
	cols = rep(class.col, each = 16)
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
	  #      las=1, adj=1, xpd=NA, family='Courier', cex=0.5)
	}
	
	maj.class.names = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
	if (classes == TRUE){
	  # draw lines above each class:
	  x.left = bp[seq(0,80, 16)+1] -0.5
	  x.right = bp[seq(16, 96, 16)] +0.5
	  rect(xleft=x.left, ymax*1.01, xright=x.right, ymax*1.08, col=class.col, border=NA, xpd=NA)
	  # mutation class labels at the top of the figure:
	  text((x.left+x.right)/2, ymax*1.123, labels=maj.class.names, cex=1, xpd=NA)
	} else if (classes == 'top'){
    # draw lines above each class:
	  x.left = bp[seq(0,80, 16)+1] -0.5
	  x.right = bp[seq(16, 96, 16)] +0.5
	  rect(xleft=x.left, ymax*1.01, xright=x.right, ymax*1.1, col=class.col, border=NA, xpd=NA)
    ## mutation class labels at the top of the figure:
    #text((x.left+x.right)/2, ymax*1.2, labels=maj.class.names, cex=1, xpd=NA)
	} else if (classes == 'bottom'){
	  # draw lines above each class:
	  x.left = bp[seq(0,80, 16)+1] -0.5
	  x.right = bp[seq(16, 96, 16)] +0.5
	  rect(xleft=x.left, ymax*(-0.1), xright=x.right, ymax*(-0.01), col=class.col, border=NA, xpd=NA)
	  # # mutation class labels at the top of the figure:
	  # text((x.left+x.right)/2, ymax*(-0.2), labels=maj.class.names, cex=1, xpd=NA)
	} 
	
  # x axis label
  if (!is.null(mut.type)){
    text(bp, -ymax/100, labels=mut.type, cex=0.75, srt=90, adj=1, family='Courier', col='grey60', xpd=NA)
  }
  # sig name
  # display.name = gsub('Signature\\.Subs\\.0', 'SBS', display.name)
  # display.name = gsub('Signature\\.Subs\\.', 'SBS', display.name)
  display.name = gsub('Signature\\.', 'SBS', display.name)
  # display.name = paste(display.name, "v2")
  if (grid.line == TRUE){
    text(1.5, ymax*7/8, labels=display.name, adj=0, cex=1.5, font=2)
  } else{
    text(1.5, ymax*7/8, labels=display.name, adj=0, cex=1, font=2)
  }
}

# Plot 96-channel mutational signatures using COSMIC conventions. 
#
# CAUTION: the input file format, inlucing the order of rows must
# be identical to that use on the PCAWG7 site on Synapse, *and 
# ths is NOT checked*
plot_snvs_file <- function(inputfile, # Must be in csv format
                                      # Filename must end in .csv
                                      
                           also.all.on.one.page=F
                           # Set this to T to also print out all signatures in a
                           # single pdf
                           )
{
  require(Cairo)
  data <- read.csv(inputfile) # Note - no checking that inputfile is correct!
  names.data <- read.csv(inputfile, header = F, stringsAsFactors = F)
	names.data <- t(names.data[1,])
	names(data) <- names.data

	outputfile1 <- gsub('csv', 'pdf', inputfile)
	outputfile2 <- gsub('csv', 'combined.pdf', inputfile)
	
	pdf(outputfile1, width=12, height=2.75, onefile=T, useDingbats = F)
	par(mar=c(1.8,0.5,2.2,0), oma=c(0,0,0,0))
	
	## for local sigs
	for (i in 3:ncol(data)){
	  mut.type = gsub('>.*','',data[,2])
	  plot_snvs(data[,i], colnames(data)[i], mut.type=mut.type)
	}
	invisible(dev.off())
	
	## generate seperate png files
	outputfolder <- gsub('\\.csv', '', inputfile)
	
	if (!file.exists(outputfolder)){
	  dir.create(outputfolder)
	}
	
	for (i in 3:ncol(data)){
	  png(paste(outputfolder, "/", colnames(data)[i], ".png", sep=""), width=12, height=2.75, units = 'in', res = 300)
	  par(mar=c(1.8,0.5,2.2,0), oma=c(0,0,0,0))
	  mut.type = gsub('>.*','',data[,1])
	  plot_snvs(data[,i], colnames(data)[i], mut.type=mut.type)
	  dev.off()
	}
	
	if (also.all.on.one.page) {
	  ## remove artifact sigs
	  data <- data[,-c(52:67)]
	  ## remove SBS27 and SBS43
	  data <- data[,-c(34,50)]
	  cairo_pdf(outputfile2, width=11.6929, height=8.2677, onefile=T)
	  par(mfrow=c(13,4), mar=c(0.3,0.1,0.3,0.1), oma=c(2,2,2,0))
	  for (i in 3:6){
	    plot_snvs(data[,i], colnames(data)[i], classes='top', grid.line=F)
	  }
	  for (i in 7:47){
	    plot_snvs(data[,i], colnames(data)[i], classes=F, grid.line=F)
	  }
	  for (i in 48:51){
	    plot_snvs(data[,i], colnames(data)[i], classes='bottom', grid.line=F)
	  }
	  invisible(dev.off())
	}
}
