#!/usr/bin/Rscript

###########################################################################
###########################################################################
# Code to plot unstranded SBS (single base substitution) spectra and signatures
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

source("plot_sbs.R")

## plot the unstranded SBS catalogs
plot_snvs_file(inputfile)

### for Jaegil's sigs combined plot
if (FALSE){
  inputfile = '~/Desktop/ICGC-upload/Mutation_Signatures/approach_2/SignatureAnalyzer_COMPOSITE_SBS_W96.signature.csv'
  data <- read.csv(inputfile)
  outputfile2 <- gsub('csv', 'combined.pdf', inputfile)
  cairo_pdf(outputfile2, width=11.6929, height=9.2677, onefile=T)
  par(mfrow=c(15,4), mar=c(0.3,0.1,0.3,0.1), oma=c(2,2,2,0))
  for (i in 3:6){
    plot_snvs(data[,i], colnames(data)[i], classes='top', grid.line=F)
  }
  for (i in 7:58){
    plot_snvs(data[,i], colnames(data)[i], classes=F, grid.line=F)
  }
  for (i in 59:62){
    plot_snvs(data[,i], colnames(data)[i], classes='bottom', grid.line=F)
  }
  invisible(dev.off())
}


