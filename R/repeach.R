#' Repeat rows of a matrix
#'
#' @param mat a matrix
#' @param n the number of repeat rows
#' @references https://stackoverflow.com/questions/11121385/repeat-rows-of-a-data-frame
#' @export

repeach = function(mat,n){
  mat1 <- mat[rep(seq_len(nrow(mat)), each=n), ]
  mat2 = as.matrix(mat1)

  return(mat2)
}
