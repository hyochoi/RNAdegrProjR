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
