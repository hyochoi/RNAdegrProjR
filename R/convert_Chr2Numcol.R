#' Convert character columns to numeric columns in a data frame
#'
#' @param df a data frame
#' @param colnums the column numbers to be converted
#' @references https://www.geeksforgeeks.org/how-to-convert-dataframe-column-from-character-to-numeric-in-r/#
#' @export

convert_Chr2Numcol = function(df, colnums) {
  vec <- c(colnums)
  df[ , vec] <- apply(df[ , vec,drop=F], 2, function(x) as.numeric(as.character(x)))

  return(df)
}
