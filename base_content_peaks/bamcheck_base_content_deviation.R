###############################################################################
# bamcheck_base_content_peaks.R
###############################################################################
# 
# Copyright (c) 2012, 2013 Genome Research Ltd.
# 
# Author: Joshua Randall <joshua.randall@sanger.ac.uk>
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
R --vanilla --slave --args bamcheck=1234_5#6.bam.bamcheck outfile=1234_5#6.bam.bamcheck runmed.k=25 baseline.method="runmed" < bamcheck_indel_peaks.R'


###############################################################################
# Description
###############################################################################
#
# This R script processes the ACGT content per cycle (GCC) section of 
# bamcheck [1] output. It fits a baseline separately to each of A, C, G, T 
# base content, substracts the baseline from the percentage at each read cycle, 
# and then calculates the percentage above and below the baseline. It outputs 
# the overall percent above and below (across all read cycles) as key/value pairs 
# that can be used as thresholds for downstream QC. 
#
# Several methods are available for fitting the baseline. The default is 
# baseline.method="median" which uses a horizontal line at Y=median as the 
# baseline. This assumes that the base content at each base will not 
# systematically change across all read cycles. You can also use 
# baseline.method="mean" to fit a horizontal line at Y=mean if you prefer.
#
# An alternative if the content shape has systematic bias (and that is not 
# a problem in and of itself) whould be to use baseline.method="runmed" which 
# uses a running median. This is calculated over a sliding window of size 
# runmed.k=25 (by default, but this can be specified as a command-line 
# parameter). This baseline estimator does not make any assumptions about 
# the baseline shape.
#
# Normal runs should be expected to fit the horiztonal line very well, having a 
# very small percentage above and below baseline, and depending on the extent to 
# which small peaks could influence the dataset, a threshold of <1% might be 
# appropriate for passing lanes. 
# 
# [1] https://github.com/VertebrateResequencing/vr-codebase/blob/master/scripts/bamcheck.c
###############################################################################

###############################################################################
# Required libraries
###############################################################################
require(plyr)
require(reshape2)

###############################################################################
# Default Parameters
###############################################################################
runmed.k=25
baseline.method="mean"

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
  die("You must specify an input bamcheck file!", -1)
}
if(!exists("outfile")) {
  die("You must specify an output bamcheck file!", -1)
}

###############################################################################
# Protect filename from shell special characters in pipe command
# (does not support spaces in filenames for now)
###############################################################################
protect.shell.pattern='[^\\w\\-\\@\\.\\#/]'
bamcheck.file <- gsub(pattern=protect.shell.pattern, replacement='', x=bamcheck, perl=T)

###############################################################################
# Read bamcheck data using pipe grep to extract ACGT content (GCC) data only
###############################################################################
gcc.data <- read.table(pipe(paste("grep '^GCC'", bamcheck.file)), colClasses=c("NULL","integer","numeric","numeric","numeric","numeric"), col.names=c("GCC","read.cycle","base.content.a","base.content.c","base.content.g","base.content.t"))

# older versions of samtools had a bug in which GCC read cycle started with 0 instead of 1
# if that is the case here, fix it so that it starts with 1
if(gcc.data$read.cycle[1]==0) {
  gcc.data$read.cycle <- gcc.data$read.cycle + 1
}

###############################################################################
# calc.baselines calculates baselines for given values
###############################################################################
calc.baselines <- function(df, runmed.k) {
  # "fit" a baseline of y=mean
  mean.baseline <- rep.int(x=mean(df$value), times=length(df$value))

  # "fit" a baseline of y=median
  median.baseline <- rep.int(x=median(df$value), times=length(df$value))

  # fit a baseline as a running median across a sliding window of runmed.k reads
  runmed.baseline <- runmed(df$value, runmed.k)

  return(data.frame(df, mean.baseline=mean.baseline, median.baseline=median.baseline, runmed.baseline=runmed.baseline))
}

###############################################################################
# bc.baseline.delta subtracts a baseline and provides values above and below 
###############################################################################
bc.baseline.delta <- function(df, baseline.method="mean", runmed.k=25) {
  df <- calc.baselines(df, runmed.k)
  if (baseline.method=="mean") {
    df$baseline <- df$mean.baseline
  } else if(baseline.method=="median") {
    df$baseline <- df$median.baseline
  } else if(baseline.method=="runmed") {
    df$baseline <- df$runmed.baseline
  } else {
    stop(paste(baseline.method,"is not a valid baseline.method"))
  }
  df <- within(df, {
	  values.above.baseline = ifelse((value - baseline) > 0, value-baseline, 0)
	  values.below.baseline = ifelse((value - baseline) < 0, abs(value-baseline), 0)
	})
  return(df)
}

###############################################################################
# bc.peaks calculates sums and percentages above and below baseline
###############################################################################
bc.peaks <- function(df) {
  ret.df <- data.frame(
    mean.above.baseline = mean(df$values.above.baseline),
    mean.below.baseline = mean(df$values.below.baseline),
    max.above.baseline = max(df$values.above.baseline),
    max.below.baseline = max(df$values.below.baseline)
  )
  ret.df <- within(ret.df, {
    total.mean.baseline.deviation = mean.above.baseline + mean.below.baseline
    max.baseline.deviation = max(max.above.baseline, max.below.baseline)
  })
  return(ret.df)
}

###############################################################################
# Plot bc.peaks diagnostics
###############################################################################
plot.bc.baseline <- function(gcc.data.baseline) {
  require(ggplot2)
  bc.plot <- ggplot(data=melt(gcc.data.baseline, id.vars=c("read.cycle","base","baseline","values.below.baseline","values.above.baseline")), mapping=aes(x=read.cycle)) + facet_grid(facets=.~base, scales="fixed") + geom_path(mapping=aes(y=value, colour=variable)) + scale_y_continuous(limits=c(0,100)) + ylab("Base Content %") + xlab("Read Cycle") + ggtitle(bamcheck)
  ggsave(plot=bc.plot, filename=paste(sep=".", outplotbase, "pdf"), width=12, height=5)
}


###############################################################################
# "Melt" the data into long (instead of wide) format (filtering out count==0)
###############################################################################
gcc.data.melt <- melt(gcc.data, id.vars="read.cycle", variable.name="base")


###############################################################################
# Calculate A, C, G, T base content percentages over and under baseline
###############################################################################
gcc.data.baseline <- ddply(.data=gcc.data.melt, .variables=c("base"), .fun=bc.baseline.delta, baseline.method=baseline.method, runmed.k=runmed.k)
gcc.data.peaks <- ddply(.data=gcc.data.baseline, .variables=c("base"), .fun=bc.peaks)
gcc.data.peaks.melt <- melt(gcc.data.peaks)

###############################################################################
# Output bc.percents and percentages as bamcheck-style Summary Number (SN) entries
###############################################################################
outdata <- data.frame(section="SN", variable=paste(sep=".", gcc.data.peaks.melt$base, gcc.data.peaks.melt$variable), value=gcc.data.peaks.melt$value)
write.table(file=outfile, x=outdata, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE)


###############################################################################
# Optionally output plots
###############################################################################
if(exists("outplotbase")) {
  plot.bc.baseline(gcc.data.baseline)
}

