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
