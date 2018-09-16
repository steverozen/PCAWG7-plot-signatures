#!/usr/bin/Rscript

###########################################################################
###########################################################################
# Code to plot insertion/deletion (indel, ID) spectra and signatures
#
# 2018 09 14
# 
# v0.1
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

## input: indel.catalog.csv 
args = commandArgs(trailingOnly=TRUE)
inputfile = args[1]

# indel.class.col  <- c('#FFE5CC',
#                       '#FFB266',
#                       '#CCFFCC',
#                       '#66FF66',
#                       '#FFCCCC',
#                       '#FF9999',
#                       '#FF6666',
#                       '#FF3333',
#                       '#CCE5FF',
#                       '#99CCFF',
#                       '#66B2FF',
#                       '#3399FF',
#                       '#E5CCFF',
#                       '#CC99FF',
#                       '#B266FF',
#                       '#9933FF')

indel.class.col  <- c('#fdbe6f',
                      '#ff8001',
                      '#b0dd8b',
                      '#36a12e',
                      '#fdcab5',
                      '#fc8a6a',
                      '#f14432',
                      '#bc141a',
                      '#d0e1f2',
                      '#94c4df',
                      '#4a98c9',
                      '#1764ab',
                      '#e2e2ef',
                      '#b6b6d8',
                      '#8683bd',
                      '#61409b')


## plot function
plot_indels <- function(counts, display.name, mut.type=NULL, classes=TRUE, grid.line=TRUE){
	num.classes = length(counts)
	cols = rep(indel.class.col, 
	           c(6,6,6,6,
	             6,6,6,6,
	             6,6,6,6,
	             1,2,3,5))
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
	
	#  mutation class labels and lines above each class
	maj.class.names = c("1bp deletion", "1bp insertion", ">1bp deletions at repeats\n(Deletion length)",
	                    ">1bp insertions at repeats\n(Insertion length)", "Deletions with microhomology\n(Deletion length)")
	#maj.class.names.2 = c("1bp \ndeletion", "1bp\ninsertion", ">1bp deletions\nat repeats\n(Deletion length)",
	#                    ">1bp insertions\nat repeats\n(Insertion length)", "Deletions with\nmicrohomology\n(Deletion length)")
	#maj.class.names.3 = c("1bp\ndeletion", "1bp\ninsertion", "(Deletion length)\>1bp deletions\nat repeats",
	#                     "(Insertion length)\>1bp insertions\nat repeats", "(Deletion length)\nDeletions with\nmicrohomology")
	
	x.left  <-  bp[c(seq(0, 66, 6), 72, 73, 75, 78) +1] -0.5
	x.right  <- bp[c(seq(6, 72, 6), 73, 75, 78, 83)] +0.5
	class.pos  <- c( (x.left[seq(1, 4, 2)] + x.right[seq(2, 5, 2)])/2,
	                 (x.left[c(5, 9)] + x.right[c(8, 12)])/2, 
	                 (x.left[13] +x.right[length(x.left)])/2)
	
	## headings in rects
	category.lab  <-  c(rep(c('C', 'T'), 2), rep(c('2', '3', '4', '5+'), 3)) 
	category.col  <-  c(rep(c('black', 'white'), 2), rep(c('black', 'black', 'black', 'white'), 3)) 
	
	if (classes == TRUE){
	  # draw lines above each class:
	  rect(xleft=x.left, ymax*1.01, xright=x.right, ymax*1.09, col=indel.class.col, border=NA, xpd=NA)
	  text((x.left + x.right)/2, ymax*1.05, labels=category.lab, cex=0.8, col=category.col, xpd=NA)
	  # mutation class labels at the top of the figure:
	  text(class.pos, ymax*1.2, labels=maj.class.names, cex=1, xpd=NA)
	} else if (classes == 'top'){
    # draw lines above each class:
	  rect(xleft=x.left, ymax*1.01, xright=x.right, ymax*1.15, col=indel.class.col, border=NA, xpd=NA)
	  # text((x.left + x.right)/2, ymax*1.08, labels=category.lab, cex=0.5, xpd=NA)
    ## mutation class labels at the top of the figure:
    # text(class.pos, ymax*1.4, labels=maj.class.names.2, cex=0.6, xpd=NA)
	} else if (classes == 'bottom'){
	  ## draw lines above each class:
	  rect(xleft=x.left, ymax*(-0.15), xright=x.right, ymax*(-0.01), col=indel.class.col, border=NA, xpd=NA)
	  # text((x.left + x.right)/2, ymax*(-0.08), labels=category.lab, cex=0.5, xpd=NA)
	  ## mutation class labels at the top of the figure:
	  # text(class.pos, ymax*(-0.4), labels=maj.class.names.3, cex=0.6, xpd=NA)
	} 
	
  # x axis label
  if (is.null(mut.type)){
    mut.type  <-c (rep( c('1', '2', '3', '4', '5', '6+'), 2), rep( c('0', '1', '2', '3', '4', '5+'), 2), 
                   rep( c('1', '2', '3', '4', '5', '6+'), 4), rep( c('0', '1', '2', '3', '4', '5+'), 4), 
                   '1', '1', '2', '1', '2','3', '1', '2', '3', '4', '5+')
    ## V1
    rect(xleft=x.left, -ymax*0.09, xright=x.right, -ymax*0.01, col=indel.class.col, border=NA, xpd=NA)
    ## V2
    # rect(xleft=x.left-0.5, -ymax*0.09, xright=x.right+0.5, -ymax*0.01, col=c('white', 'grey90'), border=NA, xpd=NA)
    
    text(bp, -ymax/100, labels=mut.type, cex=0.75, adj=c(0.5,3.5), family='Courier', col='black', xpd=NA)
    
    # bottom.pos  <- c((x.left[1] + x.right[4])/2,
    #                  class.pos[3:length(class.pos)])
    
    # bottom.lab  <- c('Homopolyer length', 'Number of repeat units after deletion',
    #                  'Number of repeat units before insertion', 'Microhomology length')
    
    bottom.pos  <- c((x.left[1] + x.right[2])/2, (x.left[3] + x.right[4])/2,
                     class.pos[3:length(class.pos)])
    
    bottom.lab  <- c('Homopolymer length', 'Homopolymer length', 'Number of repeat units',
                     'Number of repeat units', 'Microhomology length')
                     
    text(bottom.pos, -ymax/100, labels=bottom.lab, cex=0.8, adj=c(0.5,4.25), xpd=NA)
  }
  # sig name
	display.name = gsub('Signature\\.Indels\\.0', 'ID', display.name)
	display.name = gsub('Signature\\.Indels\\.', 'ID', display.name)
	if (grid.line == TRUE){
    text(1.5, ymax*7/8, labels=display.name, adj=0, cex=1.5, font=2)
  } else{
    text(1.5, ymax*7/8, labels=display.name, adj=0, cex=1, font=2)
  }
}

plot_idel_file <- function(inputfile, outputfile = NULL){
	data <- read.csv(inputfile)
	names.data <- read.csv(inputfile, stringsAsFactors = F, header = F)
	names.data <- t(names.data[1,])
	names(data) <- names.data
	
	outputfile1 <- gsub('csv', 'pdf', inputfile)
	outputfile2 <- gsub('csv', 'combined.pdf', inputfile)
	pdf(outputfile1, width=12, height=3, onefile=T, useDingbats = F)
	par(mar=c(2.5,0.5,3,0), oma=c(0,0,0,0))
	for (i in 2:ncol(data)){
	  plot_indels(data[,i], colnames(data)[i])
	}
	invisible(dev.off())
	
	## generate seperate png files
	outputfolder <- gsub('\\.csv', '', inputfile)
	for (i in 2:ncol(data)){
	  png(paste(outputfolder, "/", colnames(data)[i], ".png", sep=""), width=12, height=3, units = 'in', res = 300)
	  par(mar=c(2.5,0.5,3,0), oma=c(0,0,0,0))
	  plot_indels(data[,i], colnames(data)[i])
	  dev.off()
	}
	
	cairo_pdf(outputfile2, width=11.6929, height=3.8677, onefile=T)
	par(mfrow=c(5,4), mar=c(0.3,0.1,0.3,0.1), oma=c(3.5,2,3.5,0))
  for (i in 2:5){
    plot_indels(data[,i], colnames(data)[i], mut.type=NA, classes='top', grid.line=F)
  }
	for (i in 6:14){
	  plot_indels(data[,i], colnames(data)[i], mut.type=NA, classes=F, grid.line=F)
	}
	for (i in 15:18){
	  plot_indels(data[,i], colnames(data)[i], mut.type=NA, classes='bottom', grid.line=F)
	}
	invisible(dev.off())
}

## plot the indel catalogs
plot_idel_file(inputfile)

### for Jaegil's sigs combined plot
if (FALSE){
  inputfile = '~/Desktop/ICGC-upload/Mutation_Signatures/approach_2/SignatureAnalyzer_ID.signature.csv'
  data <- read.csv(inputfile)
  outputfile2 <- gsub('csv', 'combined.pdf', inputfile)
  cairo_pdf(outputfile2, width=11.6929, height=4.8677, onefile=T)
  par(mfrow=c(8,4), mar=c(0.3,0.1,0.3,0.1), oma=c(2,2,2,0))
  for (i in 2:5){
    plot_indels(data[,i], colnames(data)[i], mut.type=NA, classes='top', grid.line=F)
  }
  for (i in 6:26){
    plot_indels(data[,i], colnames(data)[i], mut.type=NA, classes=F, grid.line=F)
  }
  for (i in 27:30){
    plot_indels(data[,i], colnames(data)[i], mut.type=NA, classes='bottom', grid.line=F)
  }
  invisible(dev.off())
}