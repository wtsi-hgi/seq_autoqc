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
R --vanilla --slave --args bamcheck="1234_5#6.bam.bamcheck" \
outdir="./" < bamcheck_quality_dropoff.R'

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
require(reshape2)
require(Hmisc)

###############################################################################
# Default Parameters
###############################################################################
outdir = "./"
high.iqr.threshold = 12
runmed.k = 15
ignore.edge.cycles = 3

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
high.iqr.threshold = as.numeric(high.iqr.threshold)

##########################################################################
# Calculate quantiles, IQR, and mean for a read cycle
##########################################################################
readcycle.stats <- function(read.df, var) {
  quantiles <- wtd.quantile(x=read.df[,var], weights=read.df$count, probs=c(0.25, 0.5, 0.75));
  q1 <- quantiles[1]
  median <- quantiles[2]
  q3 <- quantiles[3]
  iqr <- q3 - q1
  mean <- wtd.mean(x=read.df[,var], weights=read.df$count)
  retdf <- data.frame(q1, median, q3, iqr, mean)
  names(retdf) <- paste(sep=".", var, c("q1","median","q3","iqr","mean"))
  return(retdf)
}

##########################################################################
# Calculate quality stats for each read cycle
##########################################################################
quality.per.read.cycle <- function(df, dir) {
  qprc <- ddply(.data=subset(df, subset=(direction==dir)), .variables=c("read.cycle"), .fun=readcycle.stats, "quality")
  qprc$quality.median.runmed <- as.vector(runmed(qprc$quality.median, runmed.k))
  qprc$quality.mean.runmed <- as.vector(runmed(qprc$quality.mean, runmed.k))

  ########################################################################
  # Drop the first and last few read cycles (if ignore.edge.cycles > 0)
  ########################################################################
  qprc <- subset(qprc, subset=((read.cycle > ignore.edge.cycles) & (read.cycle < (max(qprc$read.cycle)-ignore.edge.cycles+1))))

  return(qprc)
}

##########################################################################
# Calculate the max number of contiguous true values
##########################################################################
max.contiguous <- function(bools) {
  embed.bools.3 <- data.frame(embed(bools,3))
  names(embed.bools.3) <- c("next", "current", "prev")

  embed.bools.3.summed <- cbind(embed.bools.3, rowSums(embed.bools.3))

  colnames(embed.bools.3.summed) <- c("next", "current", "prev", "sum")

  embed.bools.3.summed.rle <- rle(embed.bools.3.summed$sum)

  max.contiguous.read.cycles <- max(embed.bools.3.summed.rle$length[embed.bools.3.summed.rle$values == 3], 0) 
  rle.index <- which(embed.bools.3.summed.rle$lengths == max.contiguous.read.cycles & embed.bools.3.summed.rle$values == 3)-1
  start.read.cycle <- sum(embed.bools.3.summed.rle$lengths[1:(rle.index)])+1
  end.read.cycle <- start.read.cycle + max.contiguous.read.cycles - 1
		 
  return( data.frame(start.read.cycle, end.read.cycle, max.contiguous.read.cycles) )
}

##########################################################################
# Calculate the max number of contiguous declining values
##########################################################################
max.declining <- function(values) {
  embed.values.2 <- data.frame(embed(values,2))
  names(embed.values.2)<-c("next", "current")
  ret.df <- max.contiguous(embed.values.2[,1] <= embed.values.2[,2])
  ret.df$high.value <- values[ret.df$start.read.cycle]
  ret.df$low.value <- values[ret.df$end.read.cycle]
  return(ret.df)
}

##########################################################################
# Contiguous stats
##########################################################################
contiguous.quality.stats <- function(quality.stats) {
  return(data.frame(
	   high.iqr = max.contiguous(quality.stats$quality.iqr >= high.iqr.threshold),
           mean.runmed.decline = max.declining(quality.stats$quality.mean.runmed)
	  ))
}

##########################################################################
# Plot quality stats and save to files under outplotbase
##########################################################################
plot.quality.stats <- function(quality.stats.df, direction) {
  qs.plot <- ggplot(data=melt(quality.stats.df, id.vars=c("read.cycle","quality.q1","quality.q3","quality.iqr"), measure.vars=c("quality.median","quality.mean","quality.mean.runmed")), mapping=aes(x=read.cycle)) + geom_path(mapping=aes(y=value, colour=variable)) + geom_ribbon(mapping=aes(ymin=quality.q1, ymax=quality.q3), alpha=0.25) + scale_colour_brewer(palette="Dark2")
  ggsave(plot=qs.plot, filename=paste(sep=".", outplotbase, direction, "pdf"))
}

###############################################################################
# Protect filename from shell special characters in pipe command
# (does not support spaces in filenames for now)
###############################################################################
protect.shell.pattern='[^\\w\\-\\@\\.\\#\\/]'
bamcheck.file <- gsub(pattern=protect.shell.pattern, replacement='', x=bamcheck, perl=T)
message(paste("bamcheck file:", bamcheck, "high.iqr.threshold:", high.iqr.threshold, "outfile:", outfile))

###############################################################################
# Read bamcheck data using pipe egrep to extract First Fragment Quality
# and Last Fragment Quality (FFQ and LFQ) data only
###############################################################################
data <- read.table(pipe(paste("egrep '^(FFQ|LFQ)'", bamcheck.file)), as.is=T, sep="\t")
names(data)[1] <- "direction"
names(data)[2] <- "read.cycle"
names(data)[3:length(data[1,])] <- (seq_along(names(data)[3:length(data[1,])])-1)
data$direction <- factor(x=data$direction, levels=c("FFQ","LFQ"), labels=c("forward","reverse"))

###############################################################################
# "Melt" the data into long (instead of wide) format (filtering out count==0)
###############################################################################
data.melt <- subset(melt(data, id.vars=c("direction", "read.cycle"), variable.name="quality", value.name="count"), subset=(count>0))
data.melt$quality <- as.numeric(as.character(data.melt$quality))

###############################################################################
# Calculate basic stats (median, mean, quartiles, iqr) for each read cycle
###############################################################################
fwd.quality.stats <- quality.per.read.cycle(data.melt, "forward")
rev.quality.stats <- quality.per.read.cycle(data.melt, "reverse")

###############################################################################
# Calculate contiguous stats
###############################################################################
fwd.contiguous.quality.stats <- contiguous.quality.stats(fwd.quality.stats)
rev.contiguous.quality.stats <- contiguous.quality.stats(rev.quality.stats)
contiguous.quality.stats <- data.frame(quality.dropoff.fwd=fwd.contiguous.quality.stats, quality.dropoff.rev=rev.contiguous.quality.stats, quality.dropoff.high.iqr.threshold=high.iqr.threshold, quality.dropoff.runmed.k=runmed.k, quality.dropoff.ignore.edge.cycles=ignore.edge.cycles)

###############################################################################
# Output contiguous.high.iqr as bamcheck-style Summary Number (SN) entry
###############################################################################
outdata <- data.frame(section = "SN", melt(contiguous.quality.stats))
write.table(file=outfile, x=outdata, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE)

###############################################################################
# Optionally output plots
###############################################################################
if(exists("outplotbase")) {
  require(ggplot2)
  plot.quality.stats(fwd.quality.stats, "forward")
  plot.quality.stats(rev.quality.stats, "reverse")
}

