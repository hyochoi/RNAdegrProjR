## -------------------------------------------------
## Miscellaneous/Helper for Data Processing
## -------------------------------------------------

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


#' Convert a matrix to a data frame using pivot_longer
#'
#' @param mat a matrix with observations in rows, variables in columns, and values in cells
#' @param rcenames row, column, and cell names of a matrix
#' @references https://dcl-wrangle.stanford.edu/pivot-basic.html
#' @import tidyr
#' @export

convert_pivot.longer = function(mat, rcenames) {

  row <- rownames(mat)
  mat1 <- data.frame(row, mat)
  mat2 <- tidyr::pivot_longer(mat1,
                              cols = !starts_with("row"),
                              names_to = "col",
                              values_to = "cell")
  mat3 <- mat2[order(mat2$col, decreasing = FALSE), ]
  colnames(mat3) <- rcenames

  return(mat3)
}


#' Repeat rows of a matrix
#'
#' @param mat a matrix
#' @param n the number of repeat rows
#' @references https://stackoverflow.com/questions/11121385/repeat-rows-of-a-data-frame
#' @export

repeach = function(mat, n){
  mat1 <- mat[rep(seq_len(nrow(mat)), each=n), ]
  mat2 = as.matrix(mat1)

  return(mat2)
}


#' Repeat copies of array (equivalent of repmat in Matlab)
#'
#' @param X a matrix
#' @param m m copies of X are arranged in row
#' @param n n copies of X are arranged in column
#' @references Hua Liu, Jinhong You & Jiguo Cao (2023). A Dynamic Interaction Semiparametric Function-on-Scalar Model, Journal of the American Statistical Association, 118:541, 360-373, DOI: 10.1080/01621459.2021.1933496
#' @export

repmat = function(X, m, n){
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx*n)), mx*m, nx*n, byrow=T)
}

                      
#' Adjust y-axis for a plot
#'
#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/yaxis.hy.R
#' @export

yaxis.hy <- function(mat){
  #  mat : d by n matrix
  tempmax <- max(mat) ;
  tempmin <- min(mat) ;
  templen <- tempmax-tempmin ;
  return(c(tempmin-0.002*templen, tempmax+0.002*templen)) ;
}
