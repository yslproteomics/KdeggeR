#' Remove Outlier peptides
#'
#' Removes rows whose median is has a probability greater than ptreshold in the empirical frequency distribution of all rows' median. 
#' Not recommended for less than 5 rows.
#'
#' @param x a numeric matrix or data.frame, with peptides as rows.
#' @param min.size the minimum number of rows in x before outlier removal is performed.
#' @param pthreshold the probability threshold above which values are removed.
#'
#' @return an object of the same class and ncol as x.
#'
#' @export
removeOutlierPeptides <- function(x, min.size=3, pthreshold=0.0001){
  
  x <- x[which(apply(x, 1, FUN=function(z){ !all(is.na(x)) })), , drop = FALSE]
  
  if(nrow(x) < min.size)	return(x)
  
  md <- apply(x,1,na.rm=T,FUN=median)
  f <- ecdf(md)
  x <- x[which(f(md) > pthreshold), , drop = FALSE]
}