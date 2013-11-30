###############################################################################
# bamcheck_augment_summary.R
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

library(optparse)

option_list <- list(
  make_option(c("--indel-baseline-method"), action = "store", 
  	      type = "character", default = "runmed", 
	      help = paste("Baseline method for indel (one of 'runmed',",
	                   "'runmed_quadratic', 'quadratic', or 'runmin_quadratic')")),
  make_option(c("--indel-runmed-k"), action = "store",
  	      type = "integer", default = 25, 
	      help = "Number of read cycles to include in running median for indel peak detection"),
  make_option(c("--base-content-baseline-method"), action = "store",
  	      type = "character", default = "mean", 
	      help = paste("Baseline method for base content deviation",
	                   "(one of 'mean', 'median', or 'runmed')")),
  make_option(c("--base-content-runmed-k"), action = "store", 
              type = "integer", default = 25, 
	      help = "Number of read cycles to include in running median for base content deviation"),
  make_option(c("--quality-dropoff-runmed-k"), action = "store", 
              type = "integer", default = 25, 
     	      help = "Number of read cycles to include in running median for quality dropoff"),
  make_option(c("--quality-dropoff-ignore-edge-cycles"), action = "store", 
              type = "integer", default = 3, 
	      help = "Number of cycles at beginning and end of read to exclude from quality dropoff"),
  make_option(c("--quality-dropoff-high-iqr-threshold"), action = "store", 
  	      type = "numeric", default = 10, 
	      help = "Interquartile range above which to consider the IQR high for quality dropoff"),
  make_option(c("--plot-base-path"), action = "store", 
  	      type = "character", default ="",
  	      help = "Base path to output plots (default: no plots)"))

usage <- paste(sep="\n",
               "Usage: %prog [options] <bamcheck_input> <bamcheck_output>",
	       "Use opiton [-h|--help] to list additional options.")

options_args <- parse_args(OptionParser(option_list = option_list, usage = usage), 
	                   positional_arguments = TRUE)

if(length(options_args$args) != 2) {
  stop(paste("Must specify exactly two arguments", usage))
}

bamcheck_input <- options_args$args[1]
bamcheck_output <- options_args$args[2]


library(bamcheckr)

bamcheck <- read_bamcheck(bamcheck_input)

bamcheck$data$SN <- rbind(bamcheck$data$SN, 
		    	  indel_peaks(bamcheck = bamcheck, 
			              baseline_method = options_args$options[["indel-baseline-method"]],
				      runmed_k = options_args$options[["indel-runmed-k"]],
				      outplotbase = options_args$options[["plot-base-path"]]),
			  quality_dropoff(bamcheck = bamcheck, 
			                  runmed_k = options_args$options[["quality-dropoff-runmed-k"]],
					  ignore_edge_cycles = options_args$options[["quality-dropoff-ignore-edge-cycles"]],
					  high_iqr_threshold = options_args$options[["quality-dropoff-high-iqr-threshold"]],
					  outplotbase = options_args$options[["plot-base-path"]]), 
			  base_content_deviation(bamcheck = bamcheck,
			                         baseline_method = options_args$options[["base-content-baseline-method"]],
						 runmed_k = options_args$options[["base-content-runmed-k"]],
						 outplotbase = options_args$options[["plot-base-path"]]))

# Add provenance notes to header
bamcheck$comments$HEADER <- append(bamcheck$comments$HEADER, "# It was modified by bamcheck_augment_summary.R")
bamcheck$comments$HEADER <- append(bamcheck$comments$HEADER, paste("# R command line was: ", paste0(commandArgs(), collapse=" ")))

# Write augmented bamcheck file
write_bamcheck(bamcheck = bamcheck, bamcheck_file = bamcheck_output)
