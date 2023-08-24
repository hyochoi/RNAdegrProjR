#' Convert a continuous variable to eCDF values from ecdf function (For gene information)
#'
#' @param data a continuous variable for x-axis in the eCDF plot; a 1 x the number of genes matrix; rownames are geneSymbol
#' @param name a column name for eCDF Values (name_ecdf)
#' @references https://statisticsglobe.com/extract-ecdf-values-from-function-r
#' @references https://stats.stackexchange.com/questions/30858/how-to-calculate-cumulative-distribution-in-r
#' @export

convert_Cont2eCDF.gene = function(data, name){
  # Create ecdf function
  fun_ecdf <- ecdf(data)

  # Apply ecdf function
  my_ecdf <- fun_ecdf(data)

  # Combine x & eCDF values
  data_ecdf1 <- data.frame(data, my_ecdf)
  geneSymbol <- as.matrix(rownames(data_ecdf1))

  # Plot eCDF values
  # plot(data_ecdf1$data, data_ecdf1$my_ecdf)

  data_ecdf2 <- cbind(geneSymbol, data_ecdf1)
  data_ecdf3 <- data_ecdf2[,-c(2)]
  colnames(data_ecdf3) <- c("geneSymbol", sprintf("%s_ecdf",name))

  return(data_ecdf3)
}
