###############################################################################
# base-content-deviation.r 
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
# baseline_method="median" which uses a horizontal line at Y=median as the 
# baseline. This assumes that the base content at each base will not 
# systematically change across all read cycles. You can also use 
# baseline_method="mean" to fit a horizontal line at Y=mean if you prefer.
#
# An alternative if the content shape has systematic bias (and that is not 
# a problem in and of itself) whould be to use baseline_method="runmed" which 
# uses a running median. This is calculated over a sliding window of size 
# runmed_k=25 (by default, but this can be specified as a command-line 
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
# calc_baselines calculates baselines for given values
###############################################################################
calc_baselines <- function(df, runmed_k) {
  # "fit" a baseline of y=mean
  mean.baseline <- rep.int(x=mean(df$value), times=length(df$value))

  # "fit" a baseline of y=median
  median.baseline <- rep.int(x=median(df$value), times=length(df$value))

  # fit a baseline as a running median across a sliding window of runmed_k reads
  runmed.baseline <- runmed(df$value, runmed_k)

  return(data.frame(df, mean.baseline=mean.baseline, median.baseline=median.baseline, runmed.baseline=runmed.baseline))
}

###############################################################################
# subtract_baseline subtracts a baseline and provides values above and below 
###############################################################################
subtract_baseline <- function(df, baseline_method="mean", runmed_k=25) {
  df <- calc_baselines(df, runmed_k)
  if (baseline_method=="mean") {
    df$baseline <- df$mean.baseline
  } else if (baseline_method=="median") {
    df$baseline <- df$median.baseline
  } else if (baseline_method=="runmed") {
    df$baseline <- df$runmed.baseline
  } else {
    stop(paste(baseline_method,"is not a valid baseline_method."))
  }
  df <- within(df, {
	  values.above.baseline = ifelse((value - baseline) > 0, value-baseline, 0)
	  values.below.baseline = ifelse((value - baseline) < 0, abs(value-baseline), 0)
	})
  return(df)
}

###############################################################################
# calculate_deviation calculates sums and percentages above and below baseline
###############################################################################
calculate_deviation <- function(df) {
  ret_df <- data.frame(
    mean.above.baseline = mean(df$values.above.baseline),
    mean.below.baseline = mean(df$values.below.baseline),
    max.above.baseline = max(df$values.above.baseline),
    max.below.baseline = max(df$values.below.baseline)
  )
  ret_df <- within(ret_df, {
    total.mean.baseline.deviation = mean.above.baseline + mean.below.baseline
    max.baseline.deviation = max(max.above.baseline, max.below.baseline)
  })
  return(ret_df)
}

###############################################################################
# Plot baseline diagnostics
###############################################################################
plot_baseline_diagnostics <- function(gcc_baseline) {
  gcc_plot <- ggplot(data=melt(gcc_baseline, id.vars=c("read.cycle","base","baseline","values.below.baseline","values.above.baseline")), mapping=aes(x=read.cycle)) + facet_grid(facets=.~base, scales="fixed") + geom_path(mapping=aes(y=value, colour=variable)) + scale_y_continuous(limits=c(0,100)) + ylab("Base Content %") + xlab("Read Cycle") + ggtitle(bamcheck)
  ggsave(plot=gcc_plot, filename=paste(sep=".", outplotbase, "pdf"), width=12, height=5)
}


###############################################################################
# base_content_deviation calculates and returns summary numbers on 
# bamcheck GCC (base content) section.
###############################################################################
base_content_deviation <- function(gcc_data, baseline_method="mean", runmed_k=25, outplotbase) {

  #############################################################################
  # "Melt" the data into long (instead of wide) format (filtering out count==0)
  #############################################################################
  gcc_data.melt <- melt(gcc_data, id.vars="read.cycle", variable.name="base")

  #############################################################################
  # Calculate A, C, G, T base content percentages over and under baseline
  #############################################################################
  gcc_baseline <- ddply(.data=gcc_data.melt, .variables=c("base"), .fun=subtract_baseline, baseline_method=baseline_method, runmed_k=runmed_k)
  gcc_data_peaks <- ddply(.data=gcc_baseline, .variables=c("base"), .fun=calculate_deviation)
  gcc_data_peaks_melt <- melt(gcc_data_peaks)

  ###############################################################################
  # Optionally output plots
  ###############################################################################
  if (exists("outplotbase")) {
    plot_baseline_diagnostics(gcc_baseline)
  }

  #############################################################################
  # Output bc.percents and percentages as bamcheck-style Summary Number (SN) 
  #############################################################################
  outdata <- data.frame(section="SN", variable=paste(sep=".", gcc_data_peaks_melt$base, gcc_data_peaks_melt$variable), value=gcc_data_peaks_melt$value)
#  write.table(file=outfile, x=outdata, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE)
  return(outdata)
}

