#' compareParameters
#'
#' Runs a series of protein-level kloss calculations and compares the results. See ?compareParametersResults to plot the results.
#'
#' @param o a pSILAC object, with all peptide-level rates calculated (see ?calcAllRates)
#' @param replicates A vector of length equal to the number of samples, or a character designating the column in the design data.frame, indicating which samples are replicates of the same condition. This will be used to calculate distance/correlation between replicates.
#' @param ncores The number of cores to use. If NULL, will default to detected cores minus 1.
#' @param m A list of the comparisons to make.
#'
#' @return A list of the results, each with the following slots:
#' 'call': the function and parameters called.
#' 'kloss': the protein kloss table
#' 'params': the list of input parameters
#' 'cor': the pearson correlations between replicates.
#' 'dist': the euclidean distance between replicates.
#' 'MdAE': the median absolute error between replicates.
#'
#' @export
compareParameters <- function(o, replicates=NULL, ncores=NULL, m=NULL){
  if(class(o) != "pSILAC")	stop("o should be a pSILAC object.")
  if(is.null(o$RIA.kloss) | is.null(o$hol.kloss))	stop("Peptide-level k_loss should first be calculated.")
  if(is.null(ncores)){
    library(parallel)
    ncores <- detectCores() - 1
  }
  if(is.null(replicates))	replicates <- rep(1,length(unique(o$design$sample)))
  if(length(replicates)==1){
    if(replicates %in% colnames(o$design))	replicates <- o$design[which(!duplicated(o$design$sample)),replicates]
  }
  if(length(replicates) != length(unique(o$design$sample)))	stop("'replicates' should either indicate a column of the design data.frame, or should be a vector of length equal to the number of samples.")
  mnames <- c("method","ag.metric","ag.weights","unique.weights","in.all","removeOutliers","ncores","tryRobust")
  if(is.null(m)){
    m <- list(list("mixed","median","both",T,2,5,3,T),
	      list("mixed","mean","both",T,2,5,3,T),
	      list("RIA","remodel",NA,T,2,5,3,F),
	      list("hol","remodel",NA,T,2,5,3,T),
	      list("NLI","remodel",NA,T,2,5,3,T),
	      list("NLI","remodel",NA,T,2,5,3,F),
	      list("NLI","remodel",NA,T,2,5,3,T),
	      list("combined","remodel","variance",T,2,5,3,T),
	      list("combined","remodel","nbpoints",T,2,5,3,T),
	      list("RIA","median",NA,F,2,5,3,F),
	      list("RIA","median",NA,T,2,5,3,F),
	      list("RIA","mean","variance",T,2,5,3,F),
	      list("RIA","mean","nbpoints",T,2,5,3,F),
	      list("RIA","mean","both",T,F,5,3,F),
	      list("RIA","mean","both",T,2,5,3,F),
	      list("hol","median",NA,T,2,5,3,T),
	      list("hol","mean","variance",T,2,5,3,T),
	      list("hol","mean","nbpoints",T,2,5,3,T),
	      list("combined","median",NA,T,2,5,3,T),
	      list("combined","mean","variance",T,F,F,3,T),
	      list("combined","mean","both",T,F,F,3,T),
	      list("combined","mean","nbpoints",T,F,F,3,T),
	      list("combined","mean","variance",T,2,5,3,T),
	      list("combined","mean","both",T,2,5,3,T)
	 )
  }
  v <- cbind(1:length(replicates),1:length(replicates))
  ll <- list()
  llcor <- list()
  lldist <- list()
  llmdae <- list()
  for(i in 1:length(m)){
    params <- m[[i]]
    names(params) <- mnames
    if(params[[1]]=="mixed"){
      message(paste("calcMixedProteinsKloss(o, ",paste(paste(mnames[c(-1,-3,-8)],params[c(-1,-3,-8)],sep="="),collapse=", "),")",sep=""))
      calcMixedProteinsKloss(o, m[[i]][[2]], m[[i]][[4]], m[[i]][[5]], m[[i]][[6]], m[[i]][[7]])
    }else{
      message(paste("calcProteinsKloss(o, ",paste(paste(mnames,params,sep="="),collapse=", "),")",sep=""))
      o <- calcProteinsKloss(o, m[[i]][[1]], m[[i]][[2]], m[[i]][[3]], m[[i]][[4]], m[[i]][[5]], m[[i]][[6]], m[[i]][[7]], m[[i]][[8]])
    }
    cc <- .getReplicatesValues(cor(o$protein.kloss, use="pairwise"), replicates)
    ce <- .getReplicatesValues(as.matrix(dist(t(o$protein.kloss))), replicates)
    mdae <- .getReplicatesValues(MdAE(o$protein.kloss), replicates)
    llcor[[paste(m[[i]],collapse=".")]] <- cc
    lldist[[paste(m[[i]],collapse=".")]] <- ce
    llmdae[[paste(m[[i]],collapse=".")]] <- mdae
    ll[[i]] <- list( call=o$info$protk.call, kloss=o$protein.kloss, params=params, cor=cc, dist=ce, MdAE=mdae)
  }
  layout(matrix(1:3,nrow=3))
  boxplot(llcor,ylab="Correlation",main="Comparison across replicates",las=3)
  boxplot(lldist,ylab="Euclidean distance",las=3)
  boxplot(llmdae,ylab="Median absolute error",las=3)
  return(ll)
}

#' compareParametersResults
#'
#' Plots the results of the 'compareParameters' function, taking only the proteins which are given a value in all analyses.
#'
#' @param ll a list, as produced by the compareParameters function.
#' @param replicates A vector of length equal to the number of samples, indicating which samples are replicates of the same condition. 
#' If this is contained in the design data.frame of a pSILAC object, you may retrieve it through: o$design[which(!duplicated(o$design$sample)),"replicateColumnNames"]
#'
#' @return Nothing.
#'
#' @export
compareParametersResults <- function(ll, replicates){
  g <- row.names(ll[[1]][["kloss"]])
  llcor <- list()
  lldist <- list()
  llmdae <- list()
  for(i in 2:length(ll))	g <- g[which(g %in% row.names(ll[[i]][["kloss"]]))]
  lnames <- gsub("TRUE","T",gsub("FALSE","F",as.character(lapply(ll, FUN=function(x){ paste(x$params,collapse=".") })),fixed=T),fixed=T)
  llcor <- lapply(ll, FUN=function(x){ .getReplicatesValues(cor(x$kloss[g,], use="pairwise"), replicates) })
  lldist <- lapply(ll, FUN=function(x){ .getReplicatesValues(as.matrix(dist(t(x$kloss[g,]))), replicates) })
  llmdae <- lapply(ll, FUN=function(x){ .getReplicatesValues(MdAE(x$kloss[g,]), replicates) })
  names(llcor) <- lnames
  names(lldist) <- lnames
  names(llmdae) <- lnames
  layout(matrix(1:3,nrow=3))
  boxplot(llcor,ylab="Correlation",main="Comparison across replicates",las=3,outline=F)
  for(i in 1:length(llcor)) points(rep(i,length(llcor[[i]])), llcor[[i]], pch=20)
  boxplot(lldist,ylab="Euclidean distance",las=3,outline=F)
  for(i in 1:length(llcor)) points(rep(i,length(lldist[[i]])), lldist[[i]], pch=20)  
  boxplot(llmdae,ylab="Median absolute error",las=3,outline=F)
  for(i in 1:length(llcor)) points(rep(i,length(llmdae[[i]])), llmdae[[i]], pch=20)  
}

.getReplicatesValues <- function(x, replicates){
  x[cbind(1:ncol(x), 1:ncol(x))] <- NA
  return(as.numeric(unlist(lapply(unique(replicates),FUN=function(r){
    w <- which(replicates==r)
    x <- x[w,w]
    as.numeric(x[lower.tri(x)])
  }))))
}

#' Calculates median absolute error (MdAE)
#'
#' @param x Either a numeric vector, or a numeric matrix/data.frame with samples as columns, and features as rows.
#' @param y Should be NULL if x is a matrix/data.frame, or (if x is a vector) a numeric vector of length equal to that of x.
#'
#' @return Either the median absolute error (if x and y are vectors), or a symmetric matrix (of size=ncol(x)) of pair-wise MdAE if x is a matrix/data.frame.
#'
#' @export
MdAE <- function(x,y=NULL){
  if(!is.null(y)){
    if(is.matrix(x) | is.data.frame(x) | length(x) != length(y))	stop("If y is defined, x and y should be numeric vectors of the same length. If y is NULL, x should be a numeric matrix of data.frame.")
    return(median(abs(y-x),na.rm=T))
  }
  m <- matrix(0,ncol=ncol(x),nrow=ncol(x))
  for(i in 1:ncol(x)){
    for(j in 1:ncol(x)){
      if(i!=j)	m[i,j] <- MdAE(x[,i],x[,j])
    }
  }
  colnames(m) <- colnames(x)
  rownames(m) <- colnames(x)
  return(m)
}

#' Calculates median relatibe error (MdRE)
#'
#' @param x Either a numeric vector, or a numeric matrix/data.frame with samples as columns, and features as rows.
#' @param y Should be NULL if x is a matrix/data.frame, or (if x is a vector) a numeric vector of length equal to that of x.
#'
#' @return Either the median relative error (if x and y are vectors), or a symmetric matrix (of size=ncol(x)) of pair-wise MdRE if x is a matrix/data.frame.
#'
#' @export
MdRE <- function(x,y=NULL){
  if(!is.null(y)){
    if(is.matrix(x) | is.data.frame(x) | length(x) != length(y))	stop("If y is defined, x and y should be numeric vectors of the same length. If y is NULL, x should be a numeric matrix of data.frame.")
    return(median(2*abs(y-x)/(y+x),na.rm=T))
  }
  m <- matrix(0,ncol=ncol(x),nrow=ncol(x))
  for(i in 1:ncol(x)){
    for(j in 1:ncol(x)){
      if(i!=j)	m[i,j] <- MdRE(x[,i],x[,j])
    }
  }
  colnames(m) <- colnames(x)
  rownames(m) <- colnames(x)
  return(m)
}

getPeptideConcordance <- function(o, method="RIA", ag.func=mad){
	slot <- paste(method,"kloss",sep=".")
	wr <- intersect(row.names(o[[slot]]),row.names(o$peptides)[which(o$peptides$nbProteins==1)])
	wc <- grep("kloss$",colnames(o[[slot]]))	
	row.names(o$peptides)[which(o$peptides$nbProteins==1)]
	e <- aggregate(o[[slot]][wr,wc], by=list(protein=o$peptides[wr,"protein"]), na.rm=T, FUN=ag.func)
	row.names(e) <- e$protein
	e$protein <- NULL
	return(e)
}


# to double-check
.errProp <- function(sds, values=NULL){
	sds <- sqrt( sum(sds^2,na.rm=T)/(sum(!is.na(sds))^2) )
	if(is.null(values))	return(sds)
	mm <- suppressWarnings(max(sd(values,na.rm=T), sds, na.rm=T))
	if(is.infinite(mm)) return(NA)
	return(mm)
}

