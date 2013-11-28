
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

