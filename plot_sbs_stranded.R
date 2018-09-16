#!/usr/bin/Rscript

###########################################################################
###########################################################################
# Code to plot *stranded* SBS (single base substitution) spectra and signatures
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

maj.class.names = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

class.col  <- c('#03bcee',
                '#010101',
                '#e32926',
                '#999999',
                '#a1ce63',
                '#ebc6c4')

bg.class.col  <- c('#DCF8FF',
                   '#E9E9E9',
                   '#FFC7C7',
                   '#F7F7F7',
                   '#E5F9DF',
                   '#F9E7E7')

strand.col <- c('#394398', 
                '#e83020')

## plot function
plot_snvs_stranded <- function(counts, display.name, mut.type=NULL, classes=TRUE){
	num.classes = length(counts)
	cols = rep(strand.col, num.classes/2)
	## get ylim
	ymax = ifelse(max(counts) * 1.3 > 1, 1, max(counts) * 1.3)
	# barplot: side by side
	mat = matrix(counts, nrow = 2, ncol = num.classes/2)
	bp = barplot(mat, beside = T, main='', xaxt='n', ylim=c(0,ymax), xlab=NA, yaxt='n', ylab=NA, 
	             lwd=if(num.classes<200)3 else 1, border=NA, col=cols, xpd=NA)

	# draw lines above each class:
	x.left = bp[seq(0, 160, 32)+1] -0.5
	x.right = bp[seq(32, 192, 32)] +0.5
	rect(xleft=x.left, ymax*1.01, xright=x.right, ymax*1.08, col=class.col, border=NA, xpd=NA)
	# mutation class labels at the top of the figure:
	text((x.left+x.right)/2, ymax*1.15, labels=maj.class.names, cex=1.5, xpd=NA)
	
	# background color
	rect(xleft=x.left-0.5, 0, xright=x.right+0.5, ymax, col=bg.class.col, border='grey90', lwd = 1.5)
	
	# plot again ...
	barplot(mat, beside = T, main='', xaxt='n', ylim=c(0,ymax), xlab=NA, yaxt='n', ylab=NA, 
	        lwd=if(num.classes<200)3 else 1, border=NA, col=cols, xpd=NA, add=T)
	
	# draw grid lines
	segments(bp[1]-1, seq(0,ymax,ymax/4), bp[num.classes]+1, seq(0,ymax,ymax/4), col='grey35', lwd = 0.25)
  # draw y axis
	y.axis.values = seq(0,ymax,ymax/4)
	y.axis.labels = format(round(y.axis.values*100, 1), nsmall = 1)
	text(0.35, y.axis.values, labels=paste(y.axis.labels, '%', sep=''), 
	     las=1, adj=1, xpd=NA, cex=1)
	
  # x axis label
  if (!is.null(mut.type)){
    context.pos <- (bp[seq(1,191,2)] + bp[seq(2,192,2)])/2
    text(context.pos, -ymax/100, labels=mut.type, cex=0.85, srt=90, adj=1, family='Courier', col='grey60', xpd=NA)
  }
  # sig name
  text(1.5, ymax*7/8, labels=display.name, adj=0, cex=1.5, font=2)
}

plot_side_strand <- function(SB.counts){
  num.classes = length(SB.counts)
  cols = rep(strand.col, num.classes/2)
  ## get xlim
  xmax = ifelse(max(SB.counts) * 1.3 > 1, 1, max(SB.counts) * 1.3)
  xmax = ifelse(xmax < 0.5, 0.5, xmax)
  # barplot: side by side
  mat = matrix(SB.counts, nrow = 2, ncol = num.classes/2)
  mat = mat[2:1,6:1]
  bp = barplot(mat, beside = T, horiz = T, main='', xaxt='n', xlim=c(0,xmax), xlab=NA, yaxt='n', ylab=NA, 
               lwd=if(num.classes<200)3 else 1, border=NA, col=cols[2:1], xpd=NA)
  # draw x axis
  x.axis.max = ceiling (xmax * 10)/10
  x.axis.labels = seq(0,x.axis.max,x.axis.max/2)
  segments(0, 0.5, x.axis.max, 0.5, col='black', lwd = 1.5, xpd=NA)
  text(x.axis.labels, -1, labels=paste0(x.axis.labels*100, '%'), las=1, xpd=NA, cex=1)
  # draw y axis
  segments(0, 0.5, 0, bp[num.classes]+1, col='black', lwd = 1.5, xpd=NA)
  text(-x.axis.max*0.01, (bp[seq(1,11,2)]+bp[seq(2,12,2)])/2, labels=rev(maj.class.names), 
       adj=1, las=2, xpd=NA, cex=1.5)
  # add legend
  y.pos <- ifelse(which(mat==max(mat)) <= 6, bp[num.classes], bp[6])
  legend(x.axis.max*0.4, y.pos, fill = strand.col, xpd=NA, bty='n',
         legend = c('Transcribed strand', 'Untranscribed strand'), cex=1.5)
  
}



### test ###
inputfile <- '~/Desktop/PCAWG7-calls/Ludmil-SBS-strand-bias/signatures_strand_bias.csv'
data <- read.csv(inputfile, header = T)
names(data) <- gsub('Signature\\.Subs\\.0', 'SBS', names(data))
names(data) <- gsub('Signature\\.Subs\\.', 'SBS', names(data))
names(data) <- gsub('\\.', '-', names(data))
## switch SBS44 and SBS60
data[,c(52,68)] <- data[,c(68,52)]
write.csv(data, '~/Desktop/PCAWG7-calls/Ludmil-SBS-strand-bias/sigProfiler_TSB_signatures.csv', quote=F, row.names = F)

### sort data in the plotting order
data <- data[order(data[,'Strand']),]
data <- data[order(data[,'Subtype']),]
data <- data[order(data[,'Type']),]

## get mut type context
mut.type <- data[order(data[,'Strand']),][1:96, 'Subtype']

for (i in 4:ncol(data)){
  subdata <- data[,c(1:3,i)]
  sig <- names(data)[i]
  counts <- data[,sig]
  ### sort strand bias counts: T/U
  subdata <- subdata[order(subdata[,'Strand']),]
  subdata <- subdata[order(subdata[,'Type']),]
  SB.counts <- unlist(lapply(c(1:12), function(i) sum(subdata[(16*(i-1)+1):(16*i),sig])))
  outname <- paste('~/Desktop/PCAWG7-calls/Ludmil-SBS-strand-bias/SB_sigs/', sig, '-SB.png', sep='')
  png(outname, width=12, height=2, units = 'in', res = 300)
  par(mar=c(1.8,0.5,2.4,0.5), oma=c(0,0,0,4.5))
  layout(matrix(c(1,1,1,2), 1, 4, byrow = TRUE))
  plot_snvs_stranded(counts, sig, mut.type=mut.type)
  plot_side_strand(SB.counts)
  dev.off()
}


## one pdf
outname <- paste('~/Desktop/PCAWG7-calls/Ludmil-SBS-strand-bias/sigProfiler_TSB_signatures.pdf')
cairo_pdf(outname, width=12, height=2, onefile = T)
for (i in 4:ncol(data)){
  subdata <- data[,c(1:3,i)]
  sig <- names(data)[i]
  counts <- data[,sig]
  ### sort strand bias counts: T/U
  subdata <- subdata[order(subdata[,'Strand']),]
  subdata <- subdata[order(subdata[,'Type']),]
  SB.counts <- unlist(lapply(c(1:12), function(i) sum(subdata[(16*(i-1)+1):(16*i),sig])))
  par(mar=c(1.8,0.5,2.4,0.5), oma=c(0,0,0,4.5))
  layout(matrix(c(1,1,1,2), 1, 4, byrow = TRUE))
  plot_snvs_stranded(counts, sig, mut.type=mut.type)
  plot_side_strand(SB.counts)
}
dev.off()



### test for exome ###
inputfile <- '~/Desktop/PCAWG7-calls/Ludmil-SBS-strand-bias/exome_strand_bias.csv'
data <- read.csv(inputfile, header = T)
names(data) <- gsub('\\.', '-', names(data))

### sort data in the plotting order
data <- data[order(data[,'Strand']),]
data <- data[order(data[,'Subtype']),]
data <- data[order(data[,'Type']),]

## get mut type context
mut.type <- data[order(data[,'Strand']),][1:96, 'Subtype']

for (i in 4:ncol(data)){
  subdata <- data[,c(1:3,i)]
  sig <- names(data)[i]
  counts <- data[,sig]
  ### sort strand bias counts: T/U
  subdata <- subdata[order(subdata[,'Strand']),]
  subdata <- subdata[order(subdata[,'Type']),]
  SB.counts <- unlist(lapply(c(1:12), function(i) sum(subdata[(16*(i-1)+1):(16*i),sig])))
  outname <- paste('~/Desktop/PCAWG7-calls/Ludmil-SBS-strand-bias/Exome_SB_sigs/', sig, '-SB.png', sep='')
  png(outname, width=12, height=2, units = 'in', res = 300)
  par(mar=c(1.8,0.5,2.4,0.5), oma=c(0,0,0,4.5))
  layout(matrix(c(1,1,1,2), 1, 4, byrow = TRUE))
  plot_snvs_stranded(counts, sig, mut.type=mut.type)
  plot_side_strand(SB.counts)
  dev.off()
}


## one pdf
outname <- paste('~/Desktop/PCAWG7-calls/Ludmil-SBS-strand-bias/sigProfiler_exome_TSB_signatures.pdf')
cairo_pdf(outname, width=12, height=2, onefile = T)
for (i in 4:ncol(data)){
  subdata <- data[,c(1:3,i)]
  sig <- names(data)[i]
  counts <- data[,sig]
  ### sort strand bias counts: T/U
  subdata <- subdata[order(subdata[,'Strand']),]
  subdata <- subdata[order(subdata[,'Type']),]
  SB.counts <- unlist(lapply(c(1:12), function(i) sum(subdata[(16*(i-1)+1):(16*i),sig])))
  par(mar=c(1.8,0.5,2.4,0.5), oma=c(0,0,0,4.5))
  layout(matrix(c(1,1,1,2), 1, 4, byrow = TRUE))
  plot_snvs_stranded(counts, sig, mut.type=mut.type)
  plot_side_strand(SB.counts)
}
dev.off()
