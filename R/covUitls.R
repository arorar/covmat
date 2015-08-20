#' This is a utility function to compare two covariance matrices
#' 
#' @details
#' This method takes in two different covariance/correlation matrices computed
#' using two different methods and visullay compares them using the ellipse plot.
#' It produces a matrix with ellipses drawn in the upper triangle. The ellipse 
#' is drawn to be a contour of a standard bivariate normal with correlation given
#' by the correlation of the two assets. One ellipse is drawn for each covariance
#' matrix.
#' 
#' @param  cov1 covariance matrix using the first method
#' @param  cov2 covariance matrix using the second method
#' @param  labels strings indicating the type of methods used for comparison
#' @param  corr flag indicating if the supplied matrices are of type covariance
#'          or correlation
#' @author Rohit Arora
#' @export
#' 
compareCov <- function(cov1, cov2, labels, corr=FALSE) {
  
    if (length(labels) != 2) stop("There must be two labels")
    
    if (!all(dim(cov1) == dim(cov2))) stop("Matrix dimensions are unequal")
    if(nrow(cov1) != ncol(cov1)) stop("Matrix is not square")
    if(nrow(cov1) == 2) stop("Need more than 2 dims")
  
    complist <- list(list(corr=corr, cov=cov1), 
                     list(corr=corr, cov=cov2))
    names(complist) <- labels
    ellipsesPlot.covfm(complist)
}