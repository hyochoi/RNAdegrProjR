#' Convert a continuous variable to eCDF values from ecdf function
#'
#' @param data a continuous variable for x-axis in the eCDF plot; a 1 x the number of genes or samples matrix; rownames are geneSymbol or NewSampleId.
#' @param name a column name for eCDF Values (name_ecdf)
#' @param margin 1 and 2 return eCDF values for gene information and sample information, respectively.
#' @references https://statisticsglobe.com/extract-ecdf-values-from-function-r
#' @references https://stats.stackexchange.com/questions/30858/how-to-calculate-cumulative-distribution-in-r
#' @export

convert_Cont2eCDF = function(data, name, margin){

  # Create ecdf function
  fun_ecdf <- ecdf(data)

  # Apply ecdf function
  my_ecdf <- fun_ecdf(data)

  # Combine x & eCDF values
  data_ecdf1 <- data.frame(data, my_ecdf)
  data_ecdf2 <- cbind(as.matrix(rownames(data_ecdf1)), data_ecdf1)
  data_ecdf3 <- data_ecdf2[,-c(2)]

  if (margin==1) {
    colnames(data_ecdf3) <- c("geneSymbol", sprintf("%s_ecdf",name))

  } else if (margin==2) {
    colnames(data_ecdf3) <- c("NewSampleId", sprintf("%s_ecdf",name))

  } else {
    stop(margin," is not an option for margin.")
  }

  return(data_ecdf3)
}
