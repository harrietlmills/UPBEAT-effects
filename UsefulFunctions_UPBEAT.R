
format_nicely <- function(col_toformat, dp_report=3){
  
  which_lessdp = which(abs(col_toformat)<(10^(-dp_report)))
  
  new_col = formatC(col_toformat, digits =dp_report,format="f", flag="#")
  
  new_col[which_lessdp] = formatC(col_toformat[which_lessdp], digits=dp_report,format="e", flag="#")
  
  return(new_col)
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
