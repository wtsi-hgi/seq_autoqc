###############################################################################
# bamcheck_indel_peaks.R
###############################################################################
# 
# Copyright (c) 2012 Genome Research Ltd.
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
R --vanilla --slave --args bamcheck=1234_5#6.bam.bamcheck output=1234_5#6.bam.bamcheck k=25 baseline.method="runmed" < bamcheck_indel_peaks.R'


###############################################################################
# Description
###############################################################################
#
# This R script processes the Indel Count (IC) section of bamcheck [1] output. 
# It fits a baseline separately to insertion count vs read cycle and deletion 
# count vs read cycle, substracts the baseline from the counts, and then counts 
# the number of insertions and deletions that are above and below the baseline 
# and calculates the percentage of insertions that are above and below the 
# baseline and the percentage of deletions that are above and below the 
# baseline. It outputs both the counts above and below baseline and the percent 
# above and below as key/value pairs that can be used as thresholds for 
# downstream QC. 
#
# Several methods are available for fitting the baseline. The default is 
# baseline.method="runmed" which uses a running median. This is calculated over 
# a sliding window of size k=25 (by default, but this can be specified as a 
# command-line parameter). This baseline estimator is robust to peaks in  
# indel count vs read cycle of the size we have observed (normally <10 reads 
# wide) and does not make any assumptions about the baseline shape.
#
# A good alternative is baseline.method="runmed_quadratic" which uses the same 
# procedure as above but then fits a second degree polynomial to the medians, 
# forcing the baseline to a quadratic shape. In practice, this appears to fit 
# very well to indel count vs read cycle observations, but I'm not sure what 
# the reason is behind this, so I think it is safer to use the robust running 
# median as the default. 
#
# When a high percentage of indels (say >20%) are above or below the baseline, 
# this indicates a likely problem with one or more read cycles and a manual 
# inspection of the indel count vs read cycle plot can be conducted to identify 
# which read cycle(s) are effected (although this could also be automated). 
#
# Normal runs should be expected to have a small percentage of indels above and 
# below baseline, and depending on the extent to which small peaks could influence 
# the dataset, a threshold of <3% might be appropriate for passing lanes. 
# 
# In between those two extremes is a more ambiguous area in which further 
# investigation (i.e. by visual inspection of the indel count vs read cycle plot) 
# may be warranted. 
#
# [1] https://github.com/VertebrateResequencing/vr-codebase/blob/master/scripts/bamcheck.c
###############################################################################


###############################################################################
# Default Parameters
###############################################################################
k=25
baseline.method="runmed"

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

###############################################################################
# Protect filename from shell special characters in pipe command
# (does not support spaces in filenames for now)
###############################################################################
protect.shell.pattern='[^\\w\\-\\@\\.\\#/]'
bamcheck.file <- gsub(pattern=protect.shell.pattern, replacement='', x=bamcheck, perl=T)

###############################################################################
# Read bamcheck data using pipe grep to extract Indel Count (IC) data only
###############################################################################
ic.data <- read.table(pipe(paste("grep '^IC'", bamcheck.file)), colClasses=c("NULL","integer","integer","integer","integer","integer"), col.names=c("IC","read.cycle","fwd.insertion.count","rev.insertion.count","fwd.deletion.count","rev.deletion.count"))


###############################################################################
# indel.peaks subtracts a baseline and counts indels above and below baseline
###############################################################################
indel.peaks <- function(read.cycle, count, baseline.method="runmed", k=25) {
  if(baseline.method=="runmed") {
    # baseline is a running median across a sliding window of k reads
    baseline <- runmed(count, k)
  } 
  else if(baseline.method=="runmed_quadratic") {
    # quadratic fitted to running median across a sliding window of k reads
    count_runmed <- runmed(count, k)
    baseline <- lm(count_runmed~poly(read.cycle,degree=2))$fitted.values
  } 
  else if(baseline.method=="quadratic") {
    # quadratic fit to in/del counts (not recommended as it will overestimate baseline)
    warn("quadratic baseline.method is not recommended as it will likely overestimate baseline in the presence of peaks, resulting in an undercount of indels above baseline")
    baseline <- lm(count~poly(read.cycle,degree=2))$fitted.values
  } 
  else if(baseline.method=="runmin_quadratic") {
    # quadratic fit to running minimum of in/del counts (not recommended as it will underestimate baseline)
    warn("runmin_quadratic baseline.method is not recommended as it will underestimate baseline, resulting in an overcount of indels above baseline and an undercount of indels below baseline")
    require("caTools")
    count_runmin <- runmin(count, k)
    baseline <- lm(count_runmin~poly(read.cycle,degree=2))$fitted.values
  } 
  else {
    stop(paste(baseline.method,"is not a valid baseline.method"))
  } 
  count_minus_baseline <- count - baseline
  count_above_baseline <- sum(count_minus_baseline[count_minus_baseline>0])
  count_below_baseline <- abs(sum(count_minus_baseline[count_minus_baseline<0]))
  count_total <- sum(count)
  percent_above_baseline <- count_above_baseline / count_total
  percent_below_baseline <- count_below_baseline / count_total
  return(list(above.count=count_above_baseline, below.count=count_below_baseline, total.count=count_total, percent.above=percent_above_baseline, percent.below=percent_below_baseline, read.cycle=read.cycle, count=count, baseline=baseline, count.minus.baseline=count_minus_baseline))
}

###############################################################################
# Plot indel.peaks diagnostics
###############################################################################
# require(ggplot2)
# ggplot(data=data.frame(indel.peaks(peaky$readcycle, peaky$delcount, baseline.method="runmed")), mapping=aes(x=read.cycle)) + geom_line(mapping=aes(y=count), colour="black") + geom_line(mapping=aes(y=baseline), colour="red") + geom_line(mapping=aes(y=count.minus.baseline), colour="blue") 


###############################################################################
# Calculate insertion and deletion counts and percentages 
###############################################################################
fwd.insertion.peaks <- indel.peaks(read.cycle=ic.data$read.cycle, count=ic.data$fwd.insertion.count, baseline.method=baseline.method, k=k)
rev.insertion.peaks <- indel.peaks(read.cycle=ic.data$read.cycle, count=ic.data$rev.insertion.count, baseline.method=baseline.method, k=k)
fwd.deletion.peaks <- indel.peaks(read.cycle=ic.data$read.cycle, count=ic.data$fwd.deletion.count, baseline.method=baseline.method, k=k)
rev.deletion.peaks <- indel.peaks(read.cycle=ic.data$read.cycle, count=ic.data$rev.deletion.count, baseline.method=baseline.method, k=k)


###############################################################################
# Output counts and percentages as bamcheck-style Summary Number (SN) entries
###############################################################################
outdata <- data.frame(section="SN", variable=c("fwd percent insertions above baseline:", "fwd percent insertions below baseline:", "fwd percent deletions above baseline:", "fwd percent deletions below baseline:", "rev percent insertions above baseline:", "rev percent insertions below baseline:", "rev percent deletions above baseline:", "rev percent deletions below baseline:"), value=c(100.0*fwd.insertion.peaks$percent.above, 100.0*fwd.insertion.peaks$percent.below, 100.0*fwd.deletion.peaks$percent.above, 100.0*fwd.deletion.peaks$percent.below,100.0*rev.insertion.peaks$percent.above, 100.0*rev.insertion.peaks$percent.below, 100.0*rev.deletion.peaks$percent.above, 100.0*rev.deletion.peaks$percent.below))
write.table(file=outfile, x=outdata, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE)
