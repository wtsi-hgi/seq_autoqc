###############################################################################
# quality-dropoff.r
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
# Description
###############################################################################
#
# This R script processes the First Fragment Quality (FFQ) and Last Fragment
# Quality (LFQ) sections of bamcheck [1] output. 
#
# [1] https://github.com/VertebrateResequencing/vr-codebase/blob/master/scripts/bamcheck.c
###############################################################################


##########################################################################
# Calculate new Summary Number entries based on quality dropoff
##########################################################################
quality_dropoff <- function(bamcheck, runmed_k = 25, ignore_edge_cycles = 3, 
		            high_iqr_threshold = 12, outplotbase = "") {

  data <- remove_rightmost_missing_columns(
       	    rbind(data.frame(direction="forward", bamcheck$data$FFQ), 
                  data.frame(direction="reverse", bamcheck$data$LFQ)))

  ###############################################################################
  # "Melt" the data into long (instead of wide) format (filtering out count==0)
  ###############################################################################
  data_melt <- subset(melt(data, id.vars=c("direction", "read.cycle"), 
  	                   variable.name="quality", value.name="count"), 
		      subset=(count>0))
  levels(data_melt$quality) <- str_replace(string = levels(data_melt$quality), 
  			       pattern = '^Q', replacement = '')
  data_melt$quality <- as.numeric(as.character(data_melt$quality))

  ###############################################################################
  # Calculate basic stats (median, mean, quartiles, iqr) for each read cycle
  ###############################################################################
  fwd_quality_stats <- calculate_quality_per_read_cycle(data_melt, "forward", 
  		       				        runmed_k, ignore_edge_cycles)
  rev_quality_stats <- calculate_quality_per_read_cycle(data_melt, "reverse", 
  		       					runmed_k, ignore_edge_cycles)

  ###############################################################################
  # Calculate contiguous stats
  ###############################################################################
  fwd_contiguous_quality_stats <- calculate_contiguous_quality_stats(
  			            fwd_quality_stats, high_iqr_threshold)
  rev_contiguous_quality_stats <- calculate_contiguous_quality_stats(
  			            rev_quality_stats, high_iqr_threshold)
  contiguous_quality_stats <- data.frame(
  			        quality.dropoff.fwd = fwd_contiguous_quality_stats, 
				quality.dropoff.rev = rev_contiguous_quality_stats, 
				quality.dropoff.high.iqr.threshold = high_iqr_threshold, 
				quality.dropoff.runmed.k = runmed_k, 
				quality.dropoff.ignore.edge.cycles = ignore_edge_cycles)

  ###############################################################################
  # Optionally output plots
  ###############################################################################
  if (outplotbase != "") {
    plot_quality_stats(fwd_quality_stats, "forward")
      plot_quality_stats(rev_quality_stats, "reverse")
  }

  outdata <- data.frame(melt(contiguous_quality_stats))
  outdata$variable <- paste(outdata$variable, ":", sep="")

  return(outdata)
}


##########################################################################
# Calculate quantiles, IQR, and mean for a read cycle
##########################################################################
calculate_readcycle_stats <- function(read.df, var) {
  quantiles <- Hmisc::wtd.quantile(x=read.df[,var], weights=read.df$count, probs=c(0.25, 0.5, 0.75));
  q1 <- quantiles[1]
  median <- quantiles[2]
  q3 <- quantiles[3]
  iqr <- q3 - q1
  mean <- Hmisc::wtd.mean(x=read.df[,var], weights=read.df$count)
  retdf <- data.frame(q1, median, q3, iqr, mean)
  names(retdf) <- paste(sep=".", var, c("q1","median","q3","iqr","mean"))
  return(retdf)
}

##########################################################################
# Calculate quality stats for each read cycle
##########################################################################
calculate_quality_per_read_cycle <- function(df, dir, runmed_k, 
				             ignore_edge_cycles) {
  qprc <- ddply(.data=subset(df, subset=(direction==dir)), .variables=c("read.cycle"), .fun=calculate_readcycle_stats, "quality")
  qprc$quality.median.runmed <- as.vector(runmed(x = qprc$quality.median, k = runmed_k))
  qprc$quality.mean.runmed <- as.vector(runmed(x = qprc$quality.mean, k = runmed_k))

  ########################################################################
  # Drop the first and last few read cycles (if ignore_edge_cycles > 0)
  ########################################################################
  qprc <- subset(qprc, subset=((read.cycle > ignore_edge_cycles) & (read.cycle < (max(qprc$read.cycle)-ignore_edge_cycles+1))))

  return(qprc)
}

##########################################################################
# Calculate the max number of contiguous true values
##########################################################################
max_contiguous <- function(bools) {
  embed_bools_3 <- data.frame(embed(bools,3))
  names(embed_bools_3) <- c("next", "current", "prev")

  embed_bools_3_summed <- cbind(embed_bools_3, rowSums(embed_bools_3))

  colnames(embed_bools_3_summed) <- c("next", "current", "prev", "sum")

  embed_bools_3_summed_rle <- rle(embed_bools_3_summed$sum)

  max_contiguous_read_cycles <- max(embed_bools_3_summed_rle$length[embed_bools_3_summed_rle$values == 3], 0) 
  rle_index <- which(embed_bools_3_summed_rle$lengths == max_contiguous_read_cycles & embed_bools_3_summed_rle$values == 3)-1
  if (length(rle_index) > 0) {
    start_read_cycle <- sum(embed_bools_3_summed_rle$lengths[1:(rle_index)])+1
    end_read_cycle <- start_read_cycle + max_contiguous_read_cycles - 1
  } else {
    start_read_cycle <- 0
    end_read_cycle <- 0
  }
		 
  return( data.frame(start.read.cycle = start_read_cycle, 
  	  	     end.read.cycle = end_read_cycle, 
		     max.contiguous.read.cycles = max_contiguous_read_cycles) )
}

##########################################################################
# Calculate the max number of contiguous declining values
##########################################################################
max_declining <- function(values) {
  embed_values_2 <- data.frame(embed(values,2))
  names(embed_values_2)<-c("next", "current")
  ret_df <- max_contiguous(embed_values_2[,1] <= embed_values_2[,2])
  ret_df$high.value <- values[ret_df$start.read.cycle]
  ret_df$low.value <- values[ret_df$end.read.cycle]
  return(ret_df)
}

##########################################################################
# Contiguous stats
##########################################################################
calculate_contiguous_quality_stats <- function(quality_stats, high_iqr_threshold) {
  return(data.frame(
	   high.iqr = max_contiguous(quality_stats$quality.iqr >= high_iqr_threshold),
           mean.runmed.decline = max_declining(quality_stats$quality.mean.runmed)
	  ))
}

##########################################################################
# Plot quality stats and save to files under outplotbase
##########################################################################
plot_quality_stats <- function(quality_stats_df, direction) {
  qs_plot <- ggplot(data=melt(quality_stats_df, id.vars=c("read.cycle","quality.q1","quality.q3","quality.iqr"), measure.vars=c("quality.median","quality.mean","quality.mean.runmed")), mapping=aes(x=read.cycle)) + geom_path(mapping=aes(y=value, colour=variable)) + geom_ribbon(mapping=aes(ymin=quality.q1, ymax=quality.q3), alpha=0.25) + scale_colour_brewer(palette="Dark2")
  ggsave(plot=qs.plot, filename=paste(sep=".", outplotbase, direction, "pdf"))
}

