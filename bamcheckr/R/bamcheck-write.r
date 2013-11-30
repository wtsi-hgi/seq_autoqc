###############################################################################
# Write bamcheck file from given bamcheck object
###############################################################################
write_bamcheck <- function(bamcheck, bamcheck_file) {
  writeLines(text = interleave_comments_data_by_section(bamcheck, c("HEADER", bamcheck$section.order)), con = bamcheck_file, sep = "\n")
}

###############################################################################
# Returns a character vector, ordered by sections and including comments 
# followed by data for each section
###############################################################################
interleave_comments_data_by_section <- function(bamcheck, sections) {
  return(
    unlist(alply(.data = sections, 
          .margins = c(1), 
          .fun = function(section) {
             if (!is.null(bamcheck$data[[section]])) {
               data_lines <- aaply(.data = remove_rightmost_missing_columns(bamcheck$data[[section]]), 
	                    .margins = c(1), 
			    .fun = function(row) {
			      return(
			        paste0(c(section, row), collapse="\t")
			      )}, .expand = FALSE)
             } else {
	       data_lines <- c()
	     }
             return(c(bamcheck$comments[[section]], data_lines))
          }
    )))
}


