###############################################################################
# Usage
###############################################################################
#usage='Suggested R command line: \
#R --vanilla --slave --args bamcheck="1234_5#6.bam.bamcheck" \
#outdir="./" < bamcheck_quality_dropoff.R'

###############################################################################
# Default Parameters
###############################################################################
#outdir = "./"
#high.iqr.threshold = 12
#runmed.k = 15
#ignore.edge.cycles = 3

###############################################################################
# Process command-line arguments into variables, overriding defaults.
###############################################################################
#for (e in commandArgs(trailingOnly=TRUE)) {
#  ta = strsplit(e,"=",fixed=TRUE)
#  if(!is.null(ta[[1]][2])) {
#    assign(ta[[1]][1],ta[[1]][2])
#  } else {
#    assign(ta[[1]][1],TRUE)
#  }
#}

###############################################################################
# Check that necessary arguments have been provided
###############################################################################
#if(!exists("bamcheck")) {
#  usage("You must specify an input bamcheck file!", -1)
#}
#if(!exists("outfile")) {
#  usage("You must specify an output bamcheck file!", -1)
#}
#high.iqr.threshold = as.numeric(high.iqr.threshold)

###############################################################################
# Protect filename from shell special characters in pipe command
# (does not support spaces in filenames for now)
###############################################################################
#protect.shell.pattern='[^\\w\\-\\@\\.\\#\\/]'
#bamcheck.file <- gsub(pattern=protect.shell.pattern, replacement='', x=bamcheck, perl=T)
#message(paste("bamcheck file:", bamcheck, "high.iqr.threshold:", high.iqr.threshold, "outfile:", outfile))

quality_columns <- paste("Q", seq(from = 0, to = 93), sep = "")

bamcheck_section_columns <- list(
  SN  = c("key", "value"),
  COV = c("range.spec", "range.min", "count"),
  FFQ = c("read.cycle", quality_columns),
  GCC = c("read.cycle", "A.percent", "C.percent", "G.percent", "T.percent"),
  GCD = c("GC.percent", "unique.sequence.percentiles", "depth.percentile.10th", "depth.percentile.25th", "depth.percentile.50th", "depth.percentile.75th", "depth.percentile.90th"),
  GCF = c("GC.content", "count"),
  GCL = c("GC.content", "count"),
  IC  = c("read.cycle", "fwd.insertion.count", "rev.insertion.count", "fwd.deletion.count", "rev.deletion.count"),
  ID  = c("indel.length", "insertion.count", "deletion.count"),
  IS  = c("insert.size", "total.pair.count", "inward.oriented.pair.count", "outward.oriented.pair.count", "other.pair.count"),
  LFQ = c("read.cycle", quality_columns),
  MPC = c("read.cycle", "mismatch.count", quality_columns),
  RL  = c("read.length", "count")
)

###############################################################################
# Read bamcheck file and return a list, containing comments and data
# comments is a list of sections containing a vector of comment lines
# data is a list of sections containing a data frame of data 
###############################################################################
read_bamcheck <- function(bamcheck_file) {
  bamcheck_lines <- read_bamcheck_lines(bamcheck_file)

  comments_list <- extract_comments_by_section(bamcheck_lines$comment_lines)

  data_list <- extract_data_by_section(bamcheck_lines$data_lines)

  section_order <- extract_section_order(bamcheck_lines$comment_lines)

  return(list(comments = comments_list, data = data_list, section.order = section_order))
}


###############################################################################
# Read bamcheck file and return a list, containing comment_lines and data_lines
# both are vectors of lines from that section
###############################################################################
read_bamcheck_lines <- function(bamcheck_file) {
  bamcheck_lines <- readLines("6658_7#3.bamcheck")

  is_comment_line <- grepl(x=bamcheck_lines, pattern='^[#]')

  return(list(comment_lines = bamcheck_lines[is_comment_line], data_lines = bamcheck_lines[!is_comment_line]))
}


###############################################################################
# Takes a vector of comment lines and returns an ordered vector of section
# labels
###############################################################################
extract_section_order <- function(bamcheck_comment_lines) {

  comment_line_section_labels <- label_comment_line_sections(bamcheck_comment_lines)
  section_line_labels <- comment_line_section_labels[!is.na(comment_line_section_labels)]

  return(section_line_labels)
}

###############################################################################
# Takes a vector of comment lines and returns a vector indicating the section 
# label on that line (or NA if there is no label on the comment line)
###############################################################################
label_comment_line_sections <- function(bamcheck_comment_lines) {
  comment_section_extract_regex <- '^[#].*grep[[:space:]]+\\^([A-Z]+)[[:space:]]'

  # find section headers 
  comment_line_section_labels <- str_match(string = bamcheck_comment_lines, pattern = comment_section_extract_regex)[,2]
  
  return(comment_line_section_labels)
}


###############################################################################
# Takes a vector of comment lines and returns a list with entries named by 
# section, containing the comment lines for that section
###############################################################################
extract_comments_by_section <- function(bamcheck_comment_lines) {

  comment_line_section_labels <- label_comment_line_sections(bamcheck_comment_lines)
  section_line_indicies <- which(!is.na(comment_line_section_labels))

  # bamcheck files may start with an overall header describing processing (up until the first section header)
  bamcheck_header <- bamcheck_comment_lines[1:(section_line_indicies[1]-1)]

  # pull out section header comments
  section_comment_indicies <- rbind(
    subset(adply(.data = embed(section_line_indicies, 2), .margins = c(1), 
                 .fun = function(x) {return(data.frame(start=x[2], end=x[1]-1))}), select = c("start", "end")),
    data.frame(start=section_line_indicies[length(section_line_indicies)], end=length(bamcheck_comment_lines)))
  section_comment_indicies$section <- extract_section_order(bamcheck_comment_lines)

  comments <- dlply(.data = section_comment_indicies, .variable = c("section"), 
  		    .fun = function(df) {return(bamcheck_comment_lines[df$start:df$end])})

  comments[["HEADER"]] <- bamcheck_header

  return(comments)
}


###############################################################################
# Takes a vector of data lines and returns a list with entries named by 
# section, containing a data frame for the data in that section
###############################################################################
extract_data_by_section <- function(bamcheck_data_lines) {
  # extract vector of all data sections that occur in the data
  data_section_labels <- str_extract(string = bamcheck_data_lines, pattern = "^([A-Z]+)")
  data_sections <- levels(as.factor(data_section_labels))
  names(data_sections) <- data_sections

  # create a list comprised of data frames for each section
  data <- llply(.data = data_sections, .fun = function(section) {
    return(subset(colsplit(string = bamcheck_data_lines[data_section_labels==section], 
             pattern = "\t", 
             names = c("section", bamcheck_section_columns[[section]])), select=c(-section)))
  })

  # check for and workaround bug in older samtools/bamcheck (read.cycle should start at 1, not 0)
  if(data$GCC$read.cycle[1] == 0) {
    data$GCC$read.cycle <- data$GCC$read.cycle + 1
  }

  return(data)
}
