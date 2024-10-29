#' pSILAC object
#'
#' Creates a pSILAC object.
#'
#' @param dataset a data.frame containing peptides (both heavy and light) as row.names and samples and timepoints as columns (with the addition of a 'Protein' column).
#' @param design a data.frame where row.names correspond (optionally with "Intensity_" and "score_" prepended) to columns of 'dataset'. 
#' In addition, the 'design' object should have columns named 'sample' (character) and 'time' (numeric).
#' @param requant whether to 'keep', 'remove' or 'impute' the requantified values with a score greater than 0.05 (default 'remove', or keep if ). 
#' This requires that there be in 'dataset' columns prepended with 'Intensity_' and 'score_' for each sample/timepoint.
#' @param aggregate.replicates Function with which to aggregate replicates ('median' or 'mean'; median recommended), or NA if replicates are absent or should not be aggregated.
#' If there are replicates, they should have the same value in the column 'sample' of the design data.frame.
#' @param ncores the number of cores to use (defaults to 1)
#'
#' @return a pSILAC object.
#'
#' @export
pSILAC <- function(dataset, design, requant="remove", aggregate.replicates=NA, filterPeptides=T, ncores=1, imputeMethod="normD1"){
  if(is.character(dataset) & length(dataset)==1)	dataset <- read.delim(dataset,header=T,row.names=1,stringsAsFactors=F,na.strings=c(""," ","NA","N.A.","NaN","n.a."))
  if(is.character(design) & length(design)==1)	design <- read.delim(design,header=T,row.names=1,stringsAsFactors=F)
  if(!all(c("time","sample") %in% colnames(design)))	stop("The design data.frame should (at least) have 'sample' and 'time' columns.")
  if(!any(c("Proteins","Protein") %in% colnames(dataset)))	stop("The dataset should contain a column called 'Protein' or 'Proteins'.")
  if(!is.numeric(design$time))	stop("The 'time' column of the design data.frame should be numeric.")
  if(!is.na(aggregate.replicates)) aggregate.replicates <- match.arg(aggregate.replicates, c("median","mean"))
  if( is.na(aggregate.replicates) & 
      (length(unique(table(design$sample)))!=1 | length(unique(table(design$time)))!=1 | !all(table(paste(design$sample,design$time))==1)) 
      ){
    stop("Unless there are replicates, each sample of the design data.frame should have the same number of timepoints and each sample-time pair should appear only once. 
    If you are working with replicates, make sure to set 'aggregate.replicates' (see ?pSILAC).")
  }
  design <- design[order(design$sample, design$time),]
  ismq <- c("Proteins" %in% colnames(dataset), length(grep("heavy",row.names(dataset)))==0, length(grep("H.",colnames(dataset),fixed=T))>0 & length(grep("L.",colnames(dataset),fixed=T))>0)
  if(!all(ismq) & any(ismq))    stop("Could not identify the type of data.")
  ismq <- all(ismq)
  if(ismq){
    snh <- paste("Intensity.H.",row.names(design),sep="")
    snl <- paste("Intensity.L.",row.names(design),sep="")
    if(!all( c(snh %in% colnames(dataset), snl %in% colnames(dataset)) )){
        snh <- paste("H.",row.names(design),sep="")
        snl <- paste("L.",row.names(design),sep="")        
        if(!all( c(snh %in% colnames(dataset), snl %in% colnames(dataset)) )) stop("Could not find the light and heavy intensities for all samples of the design data.frame.")
    }
    e1 <- dataset[,snl]
    e2 <- dataset[,snh]
    row.names(e2) <- paste(row.names(e2),"heavy",sep="")
    colnames(e1) <- row.names(design)
    colnames(e2) <- row.names(design)
    e <- rbind(e1,e2)
    rm(e1,e2)
  }else{
    if(!all(row.names(design) %in% colnames(dataset))){
        if(!all(paste("Intensity",row.names(design),sep="_") %in% colnames(dataset))){
        stop(	paste("Could not find the intensity for all samples of the design data.frame. Missing samples:",
                        paste(head(row.names(design)[which(!(paste("Intensity",row.names(design),sep="_") %in% colnames(dataset)))]),"...",collapse=", ")
                    ))
        }
        # rename columns according to usual SWATH output
        e <- dataset[,paste("Intensity",row.names(design),sep="_")]
        colnames(e) <- gsub("^Intensity_","",colnames(e))
        requant <- match.arg(requant, c("remove","keep","impute"))
        if(requant=="impute")	stop("Imputation not yet implemented")  
        if(requant != "keep" & !all(paste("score",row.names(design),sep="_") %in% colnames(dataset))){
        warning("Could not find the scores for all samples. All intensities will be used without filtering.")
        requant == "keep"
        }
        if(requant=="remove"){
        s <- dataset[,paste("score",row.names(design),sep="_")]
        colnames(s) <- gsub("^score_","",colnames(s))
        s <- s[colnames(e)]
        e[s>0.05] <- NA
        message(paste("Replaced",sum(s>0.05, na.rm=T)," (",round(100*sum(s>0.05, na.rm=T)/(ncol(s)*nrow(s))),"% ) badly requantified datapoints with NAs..."))
        rm(s)
        }
        e <- e[,row.names(design)]
    }else{
        e <- dataset[,row.names(design)]
    }
  }

  if(requant=="impute")   e <- imputeAll(e, paste(design$condition,design$time), imputeMethod)
  if(filterPeptides){
    w <- unique(c(grep("K",row.names(e)),grep("R",row.names(e))))
    if(length(w) < nrow(e)){
        message(paste(nrow(e)-length(w),"peptides not containing any K or R were discarded."))
        e <- e[w,]
    }
  }
  if(!is.na(aggregate.replicates)){
    nnames <- paste(design$sample,design$time,sep=".")
    if(length(unique(nnames)) == length(nnames)){
      message("Replicate aggregation was activated, but no replicate was found. If you have replicates, make sure that the 'sample' column of the design data.frame is the same for each set of replicates.")
    }else{
      message(paste("Aggregating replicates using",aggregate.replicates,"..."))
      design <- aggregate(design,by=list(name=nnames),FUN=function(x){ x[[1]] })
      o <- order(design$sample, design$time)
      design <- design[o,]
      row.names(design) <- design$name
      e2 <- matrix(0,nrow=nrow(e),ncol=length(unique(nnames)))
      if(is.null(ncores)){
	library(parallel)
	ncores <- detectCores() - 1
      }else{
	if(ncores>1)	library(parallel)
      }
      if(ncores>1){
	library(parallel)
	cl <- makeCluster(ncores)
	clusterExport(cl, c("nnames","e","aggregate.replicates"), environment())
	e <- parSapply(cl, unique(nnames), FUN=function(x){ apply(e[,which(nnames==x)],1,na.rm=T,FUN=aggregate.replicates) })
	stopCluster(cl)
      }else{
	e <- sapply(unique(nnames), FUN=function(x){ apply(e[,which(nnames==x)],1,na.rm=T,FUN=aggregate.replicates) })
      }
      e <- e[,o]
    }
  }
  if(!("name" %in% colnames(design)))	design$name <- paste(design$sample, design$time, sep=".")
  colnames(e) <- design$name
  if(!("color" %in% colnames(design)))	design$color <- "black"
  object <- list( design=design,
		  info=list(creationDate=date(),requant=requant,aggregate.replicates=aggregate.replicates, protk.method="",pSILAC.call=match.call(),protk.call=NULL),
		  light=e[grep("heavy",row.names(e),invert=T),],
		  heavy=e[grep("heavy",row.names(e),invert=F),],
		  peptides=.getPeptideDefinition(dataset, ismq),
		  RIA=NULL,
		  hol=NULL,
		  NLI=NULL,
		  RIA.kloss=NULL,
		  hol.kloss=NULL,
		  protein.kloss=NULL
		)
  rm(e)
  if(ismq){
    rr <- c(    row.names(object$light)[which(apply(object$light,1,FUN=function(x){ !all(is.na(x)) }))],
                row.names(object$heavy)[which(apply(object$heavy,1,FUN=function(x){ !all(is.na(x)) }))]
            )
    object$light <- object$light[which(row.names(object$light) %in% rr),]
    object$heavy <- object$heavy[which(row.names(object$heavy) %in% rr),]
    object$peptides <- object$peptides[rr,]
  }
  row.names(object$heavy) <- gsub("heavy","",row.names(object$heavy),fixed=T)
  g <- row.names(object$heavy)
  g <- g[which(g %in% row.names(object$light))]
  object$RIA <- as.data.frame(object$light[g,]/(object$heavy[g,]+object$light[g,]))
  object$hol <- as.data.frame(log(object$heavy[g,]/object$light[g,]+1))
  class(object) <- "pSILAC"
  return(object)
}

.getPeptideDefinition <- function(dataset, ismq=FALSE){
  if(ismq){
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

#' Remove outlier peptides
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
  x <- x[which(apply(x,1,FUN=function(z){ !all(is.na(x)) })),]
  if(nrow(x) < min.size)	return(x)
  md <- apply(x,1,na.rm=T,FUN=median)
  f <- ecdf(md)
  x[which(f(md) > pthreshold),]
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

.getChannelSum <- function(o, inBothChannels=T){
    p <- row.names(o$light)
    if(inBothChannels) p <- p[which(p %in% row.names(o$heavy))]
    csu <- o$light[p,]
    tmp <- as.data.frame(o$heavy)[p,]
    for(i in 1:ncol(o$light)){ csu[,i] <- apply(cbind(csu[,i],tmp[,i]),1,na.rm=T,FUN=sum) }
    return(csu)
}

#' normalizeLightChannel
#' 
#' Normalizes the light channel so that the channel sum is constant
#'
#' @param o a pSILAC object.
#' @param method normalization method, either 'linear' (default) or 'geometric'.
#' @param timepoints an optional vector indicating the timepoint for each sample. If NULL (default), falls back to the value of the `time` column of the design.
#' @param removeNAs logical; whether to remove NA values before normalizing (default TRUE)
#'
#' @return data.frame.
#'
#' @export
normalizeLightChannel <- function(o, method="linear", timepoints=NULL, removeNAs=T){
  if(class(o) != "pSILAC")	stop("o should be a pSILAC object.")
  method <- match.arg(method, c("geometric","linear"))
  if(is.null(timepoints)) timepoints <- unique(o$design$time)
  csu <- .getChannelSum(o,T)
  p <- row.names(csu)
  if(removeNAs) csu <- csu[which( apply(o$heavy[p,which(o$design$time %in% timepoints)],1,FUN=function(x){!any(is.na(x))}) & 
                apply(o$light[p,which(o$design$time %in% timepoints)],1,FUN=function(x){!any(is.na(x))})
                ),]
  sizeFactors <- apply(log(csu),2,na.rm=T,FUN=median)
  if(method=="geometric"){
    normFactors <- 1/(sizeFactors/mean(sizeFactors))
    normFactors <- normFactors/exp(mean(log(normFactors)))
    o$NCS <- as.data.frame(exp(t(t(log(.getChannelSum(o,F)))*normFactors)))
    o$NLI <- as.data.frame(exp(t(t(log(o$light))*normFactors)))
  }else{
    normFactors <- exp(max(sizeFactors)-sizeFactors)
    o$NCS <- as.data.frame(t(t(.getChannelSum(o,F))*normFactors))
    o$NLI <- as.data.frame(t(t(o$light)*normFactors))
  }
  return(o)
}

medianNorm <- function(x, na.rm=T){
  sizeFactors <- as.numeric(apply(x,2,na.rm=na.rm,FUN=median))
  as.data.frame(t( mean(sizeFactors)*t(x)/sizeFactors ))
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