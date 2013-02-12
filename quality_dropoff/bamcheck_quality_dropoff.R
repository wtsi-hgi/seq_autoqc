###############################################################################
# bamcheck_quality_dropoff.R
###############################################################################
# 
# Copyright (c) 2012, 2013 Genome Research Ltd.
# 
# Authors: Mari Niemi <mn2@sanger.ac.uk> and Joshua Randall <jr17@sanger.ac.uk>
# 
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 3 of the License, or (at your option) any later 
# version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details.
# 
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see <http://www.gnu.org/licenses/>.
# 
###############################################################################

###############################################################################
# Usage
###############################################################################
usage='Suggested R command line: \
R --vanilla --slave --args bamcheck="1234_5#6.bam.bamcheck" outdir="./" < bamcheck_mari_test.R'



###############################################################################
# Description
###############################################################################
#
# This R script processes the First Fragment Quality (FFQ) and Last Fragment
# Quality (LFQ) sections of bamcheck [1] output. 
#
# [1] https://github.com/VertebrateResequencing/vr-codebase/blob/master/scripts/bamcheck.c
###############################################################################

###############################################################################
# Required libraries
###############################################################################
require(plyr)

###############################################################################
# Default Parameters
###############################################################################
outdir = "./"
iqr.threshold = 12

###############################################################################
# Process command-line arguments into variables, overriding defaults.
###############################################################################
for (e in commandArgs(trailingOnly=TRUE)) {
  ta = strsplit(e,"=",fixed=TRUE)
  if(!is.null(ta[[1]][2])) {
    assign(ta[[1]][1],ta[[1]][2])
  } else {
    assign(ta[[1]][1],TRUE)
  }
}

###############################################################################
# Usage function prints a message and dies with specified exit status
###############################################################################
usage <- function(message, status=0) {
  cat(usage)
  cat(message,"\n")
  cat("Exiting with status: ",status,"\n")
  quit(save="no", status=status, runLast=FALSE)
}

###############################################################################
# Check that necessary arguments have been provided
###############################################################################
if(!exists("bamcheck")) {
  usage()
  die("You must specify an input bamcheck file!", -1)
}
if(!exists("outfile")) {
  die("You must specify an output bamcheck file!", -1)
}
iqr.threshold = as.numeric(iqr.threshold)

###############################################################################
# Protect filename from shell special characters in pipe command
# (does not support spaces in filenames for now)
###############################################################################
protect.shell.pattern='[^\\w\\-\\@\\.\\#\\/]'
bamcheck.file <- gsub(pattern=protect.shell.pattern, replacement='', x=bamcheck, perl=T)

output.file.base <- gsub(pattern='//', replacement='/', fixed=T, x=paste(outdir,basename(bamcheck.file),sep="/"))

message(paste("bamcheck file:",bamcheck,"outdir:",outdir,"iqr.threshold:",iqr.threshold))

###############################################################################
# Read bamcheck data using pipe egrep to extract First Fragment Quality
# and Last Fragment Quality (FFQ and LFQ) data only
###############################################################################
data <- read.table(pipe(paste("egrep '^(FFQ|LFQ)'", bamcheck.file)), as.is=T, sep="\t")

calc.row.interquartile.range <- function(qualdiff) {
                   rowdata <- rep.int(x=seq(from=0,to=length(qualdiff[1,])-1), times=qualdiff[1,]);
                   lowquartile <- quantile(rowdata, probs=c(0.25));
                   highquartile <- quantile(rowdata, probs=c(0.75));
                   quartilediff <- highquartile-lowquartile
                   return(data.frame(quartilediff))
                 }

FFQ.quartilerange <- adply(.data=data[data$V1=="FFQ", -1],
                 .margins=c(1),
                 .fun=calc.row.interquartile.range)
names(FFQ.quartilerange)[1] <- "read.cycle"
FFQ.quartilerange$fwdrev = "fwd"

LFQ.quartilerange <- adply(.data=data[data$V1=="LFQ", -1],
                 .margins=c(1),
                 .fun=calc.row.interquartile.range)
names(LFQ.quartilerange)[1] <- "read.cycle"
LFQ.quartilerange$fwdrev = "rev"

rev.quartilerange <- rbind(FFQ.quartilerange, LFQ.quartilerange)



##################################################
# The calculated quartile range is in column 2 
##################################################
rev.quartilerange <- rev.quartilerange[,2]

embed.qual.range <- data.frame(embed(rev.quartilerange,3))
names(embed.qual.range)<-c("next.q", "q", "prev.q")

embed.qual.range.scored <- as.matrix(embed.qual.range.scored <- ifelse(embed.qual.range>=iqr.threshold, 1, 0))

##########################################################################
# Adding a column summing the number of cycles (0-3/3 contiguous
# cycles) with a quartile difference of more than iqr.threshold. 
#
# If you want to find out whether there's been a problem that caused massive
# differences in quality range (likely to be a hardwae failure) maybe
# count the length of the consecutive 3-point scoring cycles and see
# if it's longer than a threshold value. We'll call this the 
# contiguous cycle dropoff count. 
#
# Something like >=10 contiguous.cycle.dropoff.count could indicate failure 
# and >=5 contigious.cycle.dropoff.count could be a warning.
# 
##########################################################################

embed.scores.summed <- cbind(embed.qual.range.scored, rowSums(embed.qual.range.scored))

colnames(embed.scores.summed) <- c("next.q", "q", "prev.q", "sum.q")

embed.score.sum <- rle(embed.scores.summed[,4])

contiguous.cycle.dropoff.count <- max(embed.score.sum$length[embed.score.sum$values == 3], 0)

###############################################################################
# Output contiguous.cycle.dropoff.count as bamcheck-style Summary Number (SN) entry
###############################################################################
outdata <- data.frame(section="SN", variable=c("contiguous cycle dropoff count:"), value=c(contiguous.cycle.dropoff.count))
write.table(file=outfile, x=outdata, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE)
