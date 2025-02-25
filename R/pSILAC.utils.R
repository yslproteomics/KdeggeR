.getPeptideDefinition <- function(dataset, isStandardFormat=FALSE){
  if(isStandardFormat){
    d <- data.frame(row.names=c(row.names(dataset),paste(row.names(dataset),"heavy",sep="_")), 
                    nbProteins=rep(sapply(as.character(dataset$Proteins),FUN=function(x){ length(strsplit(x,";",fixed=T)[[1]])}),2), 
                    protein=rep(gsub(";","/",as.character(dataset$Proteins),fixed=T),2),
                    peptide=rep(row.names(dataset),2),
                    charge=ifelse("Charges" %in% colnames(dataset), rep(dataset$Charges,2), rep(NA,nrow(dataset)*2))
                )
  }else{
    d <- data.frame(row.names=row.names(dataset), 
                    nbProteins=sapply(as.character(dataset$Protein),FUN=function(x){ strsplit(x,"/",fixed=T)[[1]][[1]]}), 
                    protein=as.character(sapply(as.character(dataset$Protein),FUN=function(x){ x <- strsplit(x,"/",fixed=T)[[1]]; return(paste(x[2:length(x)],collapse="/"))}))
                    )
    d$peptide <- sapply(row.names(dataset),FUN=function(x){ strsplit(x,"_",fixed=T)[[1]][[2]]})
    d$charge <- sapply(row.names(dataset),FUN=function(x){ as.numeric(strsplit(strsplit(x,"_",fixed=T)[[1]][[3]], "", fixed=T)[[1]][[1]])})
  }
  return(d)
}

.getProtKParam <- function(x, param=NULL){
  param <- match.arg(param,c("removeOutliers","in.all"))
  if(class(x) != "pSILAC")	stop("x should be a pSILAC object.")
  if(is.null(x$info$protk.call))	return(FALSE)
  if(param %in% names(x$info$protk.call)){
    y <- x$info$protk.call[[param]]
    if(y != param)    return(y)
  }
  switch(param, removeOutliers=5,
		in.all=2,
		NULL)
}


#' print method for pSILAC objects.
#'
#' @param x a pSILAC object.
#'
#' @return nothing, but prints info on the object...
#'
#' @export
print.pSILAC <- function(x){
  if(class(x) != "pSILAC")	stop("x should be a pSILAC object.")
  cat(paste("pSILAC dataset with",length(unique(x$design$sample)),"samples and",length(unique(x$design$time)),"timepoints.\n"))
  cat(paste(nrow(x$peptides)," peptides quantified (",length(unique(x$peptides$protein[which(x$peptides$nbProteins==1)]))," unique proteins): ",nrow(x$light)," in the light channel and ",nrow(x$heavy)," in the heavy channel.\n"))
  cat(paste(nrow(x$RIA),"peptides quantified on both channels, with",round(100*sum(is.na(x$RIA))/prod(dim(x$RIA)),2),"% missing values.\n"))
  cat(x$info$protk.method)
  cat("\n")
}

#' summary method for pSILAC objects.
#'
#' @param x a pSILAC object.
#'
#' @return nothing, but prints info on the object...
#'
#' @export
summary.pSILAC <- function(x){
  print.pSILAC(x)
}

#' summary method for pSILAC objects.
#'
#' @param x a pSILAC object.
#'
#' @return the head() of each attribute of x.
#'
#' @export
head.pSILAC <- function(x){
  lapply(x,FUN=head)
}

#' getPeptides
#' 
#' Returns given peptides, or peptides corresponding to a given protein. Only one of 'peptides' or 'protein' should be set.
#'
#' @param x a pSILAC object.
#' @param sample the sample(s) for which to fetch values. NULL (default) fetches all samples. Ignored if returnValues is NA, 'RIA.kloss' or 'hol.kloss'.
#' @param peptides a list of peptide names/sequences.
#' @param protein a protein name.
#' @param returnValues values to return. Either 'light' (light channel intensities), 'heavy', 'RIA', 'hol' (log-transformed heavy/light), 'RIA.kloss', 'NLI.kloss' or 'hol.kloss'.
#' If NA (default), will return the peptides (from x$peptides).
#'
#' @return data.frame.
#'
#' @export
getPeptides <- function(x, peptides=NULL,protein=NULL, returnValues=NA, search=T){
  if(class(x) != "pSILAC")	stop("x should be a pSILAC object.")
  if(!is.na(returnValues))	returnValues <- match.arg(returnValues, c(NULL,NA,"light","heavy","RIA","hol","NLI","NCS","steadystate","hol.kloss","RIA.kloss","NLI.kloss"))
  if( !is.null(peptides) & !is.null(protein) )	stop("Only one of 'peptides' and 'protein' should be given.")
  if(!is.null(peptides)){
    notf <- sum(!(peptides %in% row.names(x$peptides) | peptides %in% x$peptides$peptide))
    w <- which(row.names(x$peptides) %in% peptides | x$peptides$peptide %in% peptides)
    if(length(w)==0)	stop("No peptide found!")
    if(notf > 0)	message(paste("Warning:",notf,"peptides could not be found in the dataset."))
  }else{
    w <- which(x$peptides$protein %in% protein)
    if(length(w)==0){
      if(!search)	return(NULL)
      w <- grep(protein,x$peptides$protein,fixed=T)
      if(length(w)>0){
	stop(paste("Protein not found. However, the following protein groups were found:",paste(head(x$peptides$protein[w]),collapse=", "),"..."))
      }else{
	stop(paste("No protein found!"))
      }
    }
    notf <- sum(!(protein %in% x$peptides$protein))
    if(notf>0)	message(paste("Warning:",notf,"proteins could not be found in the dataset."))
  }
  if(is.null(returnValues) | is.na(returnValues))	return(x$peptides[w,])
  p <- row.names(x$peptides)[w]
  p <- p[p %in% row.names(x[[returnValues]])]
  r <- x[[returnValues]][p,]
  if(is.null(dim(r))){
	r <- as.data.frame(t(r))
	row.names(r) <- p
  }
  return(r)
}

reorderSamples <- function(o, newOrder){
  if(class(o) != "pSILAC")	stop("o should be a pSILAC object.")
  o$design$sample <- as.character(o$design$sample)
  if(!all(newOrder %in% o$design$sample))	stop("Some samples specified in newOrder cannot be found!")
  newOrder2 <- as.numeric(sapply(newOrder,FUN=function(x){ which(o$design$sample == x)}))
  dropping <- sum(length(unique(o$design$sample[!(o$design$sample %in% newOrder)])))
  if(dropping > 0){
	message(paste("Dropping",dropping,"samples..."))
	if(!is.null(o$protein.kloss))	message("Since peptides were selected on the basis of the whole set of samples, if you intend to work only on this subset it might be advisable to recalculate protein-level Kloss.")
  }
  for(f in c("light","heavy","RIA","hol","NCS","NLI")){
    if(!is.null(o[[f]]))	o[[f]] <- o[[f]][,newOrder2]
  }
  newOrder3 <- as.numeric(sapply(newOrder,FUN=function(x){ which(unique(o$design$sample) == x)}))
  o$design <- o$design[newOrder2,]
  o$protein.kloss <- o$protein.kloss[,newOrder]
  for(f in c("RIA.kloss","hol.kloss","NLI.kloss")){
    if(!is.null(o[[f]])){
      o[[f]] <- o[[f]][,as.numeric(sapply(newOrder3*4,FUN=function(x){ x-3:0}))]
    }
  }
  return(o)
}


#' Normalize Data by Median Scaling
#'
#' This function normalizes columns of a matrix or data frame by scaling each column to its median.
#'
#' @param x A numeric matrix or data frame to be normalized.
#' @param na.rm Logical; if `TRUE`, NA values are removed when calculating medians. Default is `TRUE`.
#'
#' @return A data frame with normalized values, where each column is scaled by its median to ensure comparability across samples.
#'
#' @details
#' The function computes a size factor for each column based on the median value, then scales each value by these factors.
#' This normalization method is often used to adjust for variations in sample size or sequencing depth.
#'
#'
#' @export
medianNorm <- function(x, na.rm = TRUE) {
  sizeFactors <- as.numeric(apply(x, 2, na.rm = na.rm, FUN = median))
  as.data.frame(t(mean(sizeFactors) * t(x) / sizeFactors))
}


normalizeProteinKloss <- function(o, dilution.rate=NULL, removeNAs=T){
  if(class(o) != "pSILAC")	stop("o should be a pSILAC object.")
  pk <- o$protein.kloss
  if(removeNAs) pk <- pk[which( apply(pk,1,FUN=function(x){!any(is.na(x))}) ),]
  sizeFactors <- as.numeric(apply(pk,2,na.rm=T,FUN=median))
  o$NPK <- as.data.frame(t( mean(sizeFactors)*t(o$protein.kloss)/sizeFactors ))
  return(o)
}

glimpse <- function(x){
    if(!is.null(dim(x))){
        x[1:5,1:5]
    }else{
        head(x)
    }
}


getSteadystateAmounts <- function(o,normalization.method=NULL){
  if(!("NLI" %in% names(o)) | !is.null(normalization.method)){
	if(is.null(normalization.method)) normalization.method <- "linear"
	message(paste("Light channel not yet normalized; performing",normalization.method,"normalization..."))
	o <- normalizeLightChannel(o,method=normalization.method)
  }
  o$steadystate <- sapply(unique(o$design$sample),ncs=o$NCS,FUN=function(x,ncs){ as.numeric(apply(ncs[,which(o$design$sample==x)],1,FUN=.medianNonNull)) })
  row.names(o$steadystate) <- row.names(o$NCS)
  return(o)
}

summarizeTopN <- function(x, N=3, ag.fun=median){
    if(!is.null(dim(x)) | nrow(x)==1)	return(x)
    if(nrow(x)>N)	x <- x[order(apply(x,1,FUN=function(y){ sum(!is.na(y))}),apply(x,1,na.rm=T,FUN=mean),decreasing=T),]
    return(apply(x,2,na.rm=T,FUN=ag.fun))
}


getHistonesSS <- function(o, summarize=T, histones=list(H1=c("P10412","P16401","P16402","P16403","Q02539"), H4=c("P62805","Q99525"), H3=c("P68431","Q71DI3"), H2=c("Q8IUE6"))){
	h <- lapply(histones, o=o, FUN=function(x,o){ p <- getPeptides(o, protein=x, returnValues="steadystate") })
	if(summarize)	return(t(sapply(h,FUN=function(x){ if(!is.null(dim(x))) return(summarizeTopN(x)); return(x) })))
	h2 <- matrix(unlist(lapply(h,FUN=function(x){ as.numeric(t(x))})), ncol=ncol(h[[1]]), byrow=T)
	row.names(h2) <- unlist(lapply(h,FUN=row.names))
	colnames(h2) <- colnames(h[[1]])
	return(h2)
}


getHistonesK <- function(o, returnValues="protein.kloss", histones=list(H4=c("P62805","Q99525"), H3=c("P68431","Q71DI3"), H2=c("Q8IUE6"))){
	conv <- rep(names(histones),sapply(histones,FUN=length))
	names(conv) <- unlist(histones)
	if(returnValues=="protein.kloss"){
		h <- o[[returnValues]][intersect(row.names(o[[returnValues]]),unlist(histones)),]
		h$name <- conv[row.names(h)]
	}
	return(h)
}