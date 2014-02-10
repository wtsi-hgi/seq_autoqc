###############################################################################
# bamcheck_augment_summary.R
###############################################################################
# 
# Copyright (c) 2012, 2014 Genome Research Ltd.
# 
# Author: Joshua Randall <joshua.randall@sanger.ac.uk>
#         Sam Nicholls <sn8@sanger.ac.uk>
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

library(methods)
library(optparse)

options(scipen=9)
indel_baseline_methods <- c("runmed", "runmed_quadratic", "quadratic", "runmin_quadratic")
base_content_baseline_methods <- c("mean", "median", "runmed")

option_list <- list(
  make_option(c("--indel-baseline-method"), action = "store", 
  	      type = "character", default = "runmed", metavar = "method",
	      help = paste(sep = "",
	                   "Baseline method for indels (one of",
			   paste0(collapse = ", ", sprintf("'%s'", indel_baseline_methods)),
			   ")")),
  make_option(c("--indel-runmed-k"), action = "store",
  	      type = "numeric", default = 25, metavar = "k",
	      help = paste(sep=" ",
	                   "Number of read cycles to include in running median for",
	             	   "indel peak detection")),
  make_option(c("--base-content-runmed-k"), action = "store", 
              type = "numeric", default = 25, metavar = "k",
	      help = paste(sep=" ",
	                   "Number of read cycles to include in running median",
			   "for base content deviation")),
  make_option(c("--base-content-baseline-method"), action = "store",
  	      type = "character", default = "mean", metavar = "method",
	      help = paste(sep="",
	                   "Baseline method for base content deviation (one of",
			   paste0(collapse = ", ", sprintf("'%s'", base_content_baseline_methods)),
			   ")")),
  make_option(c("--quality-dropoff-runmed-k"), action = "store", 
              type = "numeric", default = 25, metavar = "k",
     	      help = paste(sep=" ",
	                   "Number of read cycles to include in running median",
			   "for quality dropoff")),
  make_option(c("--quality-dropoff-ignore-edge-cycles"), action = "store", 
              type = "numeric", default = 3, metavar = "cycle-count",
	      help = paste(sep=" ",
	                   "Number of cycles at beginning and end of read to",
			   "exclude from quality dropoff")),
  make_option(c("--quality-dropoff-high-iqr-threshold"), action = "store", 
  	      type = "numeric", default = 10, metavar = "iqr", 
	      help = paste(sep=" ",
	      	           "Interquartile range above which to consider the",
			   "IQR high for quality dropoff")),
  make_option(c("--plot-base-path"), action = "store", 
  	      type = "character", default ="", metavar = "path",
  	      help = "Base path to output plots (default: no plots)"))

usage <- "Usage: Rscript %prog [-h|--help] [options] <bamcheck_input> <bamcheck_output>"

op <- OptionParser(option_list = option_list, usage = usage)
options_args <- parse_args(op, positional_arguments = TRUE)

if(length(options_args$args) != 2) {
  error_msg <- "Unknown error"
  possible_opts = options_args$args[substr(options_args$args, 1, 1) == "-"]
  if (length(possible_opts) > 0) {
    error_msg <- paste0(sep = "\n", 
                        "Unrecognized option: ", 
	  	        possible_opts)
  } else {
    error_msg <- paste(sep = "\n",
                     sprintf("Must specify exactly two arguments (had %d)", length(options_args$args)),
    		     "Arguments were:", 
		     paste0("  ", options_args$args))
  }
  print_help(op)
  stop(call. = FALSE, error_msg)
}

stop_if_invalid_opt_nonmissing <- function(opts, optkey) {
  if (is.na(opts[[optkey]])) {
    stop(call. = FALSE, sprintf("%s is invalid or missing", paste(sep = "", "--", optkey)))
  }
}

stop_if_invalid_opt_posnum <- function(opts, optkey) {
  stop_if_invalid_opt_nonmissing(opts, optkey)
  if (! is.numeric(opts[[optkey]])) {
    stop(call. = FALSE, sprintf("%s=%s is not numeric", paste(sep = "", "--", optkey), opts[[optkey]]))
  }
  if (! opts[[optkey]] > 0) {
    stop(call. = FALSE, sprintf("%s=%s is not > 0", paste(sep = "", "--", optkey), opts[[optkey]]))
  }
}

stop_if_invalid_opt_posint <- function(opts, optkey) {
  stop_if_invalid_opt_posnum(opts, optkey)
  if (! as.integer(opts[[optkey]]) == opts[[optkey]]) {
    stop(call. = FALSE, sprintf("%s=%s is not an integer", paste(sep = "", "--", optkey), opts[[optkey]]))
  }
}

stop_if_invalid_opt_method <- function(opts, optkey, valid_methods) {
  stop_if_invalid_opt_nonmissing(opts, optkey)
  if (! opts[[optkey]] %in% valid_methods) {
    stop(call. = FALSE, sprintf("%s is not a valid %s (must be one of %s)", opts[[optkey]], paste(sep = "", "--", optkey), paste0(collapse = ", ", valid_methods)))
  }
}

# validate options
stop_if_invalid_opt_posint(options_args$options, "indel-runmed-k")
stop_if_invalid_opt_method(options_args$options, "indel-baseline-method", indel_baseline_methods)
stop_if_invalid_opt_posint(options_args$options, "base-content-runmed-k")
stop_if_invalid_opt_method(options_args$options, "base-content-baseline-method", base_content_baseline_methods)
stop_if_invalid_opt_posint(options_args$options, "quality-dropoff-runmed-k")
stop_if_invalid_opt_posint(options_args$options, "quality-dropoff-ignore-edge-cycles")

# set input and output files
bamcheck_input <- options_args$args[1]
bamcheck_output <- options_args$args[2]

# require bamcheckr package
library(bamcheckr)

# read from bamcheck file
bamcheck <- read_bamcheck(bamcheck_input)

# augment bamcheck summary numbers 
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

# add provenance notes to header
bamcheck$comments$HEADER <- append(bamcheck$comments$HEADER, "# It was modified by bamcheck_augment_summary.R")
bamcheck$comments$HEADER <- append(bamcheck$comments$HEADER, paste("# R command line was: ", paste0(commandArgs(), collapse=" ")))

# write augmented bamcheck file
write_bamcheck(bamcheck = bamcheck, bamcheck_file = bamcheck_output)

