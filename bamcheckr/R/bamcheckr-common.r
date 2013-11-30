
###############################################################################
# Usage
###############################################################################
usage.text='Suggested R command line: \
R --vanilla --slave --args bamcheck=1234_5#6.bam.bamcheck outfile=1234_5#6.bam.bamcheck runmed.k=25 baseline.method="runmed" < bamcheck_indel_peaks.R'


###############################################################################
# Usage function prints a message and dies with specified exit status
###############################################################################
usage <- function(message, status=0) {
  cat(usage.text)
  cat(message,"\n")
  cat("Exiting with status: ",status,"\n")
  quit(save="no", status=status, runLast=FALSE)
}

###############################################################################
# remove all the rightmost columns that consist of entirely missing data
###############################################################################
remove_rightmost_missing_columns <- function(df) {
  for(i in length(df[1,]):1) {
    if (!all(is.na(df[,i]))) { 
      remove_cols_from <- i + 1
      break;
    }
  }
  return(df[,1:(remove_cols_from-1)])
}


