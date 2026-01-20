#' plColorMap
#'
#' Maps colors onto numeric values.
#'
#' @param x a numeric vector.
#'
#' @return A vector of colors of the same length as x.
#'
#' @examples
#' plColorMap(1:3)
#'
#' @export
plColorMap <- function(x){
  pal <- colorRampPalette(c("blue","red"),0.5)(100)
  xmin <- min(x,na.rm=T)
  xmax <- max(x,na.rm=T)-xmin
  sapply(x,FUN=function(y){pal[1+floor((99)*(y-xmin)/xmax)]})
}

#' Plot multiple density distributions
#'
#' Maps colors onto numeric values.
#'
#' @param thelist a list or data.frame containing only number values.
#' @param cols a vector of colors, or length=length(thelist)
#' @param xlab passed to the plot function.
#' @param ylab passed to the plot function.
#' @param xlim Either a single value (max x), or a vector of length 2 passed to the plot function.
#' @param lwd line width, passed to the plot function.
#' @param ... addtional arguments passed to the plot function.
#'
#' @return Plots the density distributions.
#'
#' @examples
#' gdensity(list(A=1:3,B=2:4,C=3:5),c("red","green","blue"))
#'
#'
#' Plot multiple density distributions
#'
#' @param thelist a list or data.frame containing only number values.
#' @param cols a vector of colors, or length=length(thelist)
#' @param xlab passed to the plot function.
#' @param ylab passed to the plot function.
#' @param xlim Either a single value (max x), or a vector of length 2 passed to the plot function.
#' @param lwd line width, passed to the plot function.
#' @param main the plot title (default empty)
#' @param ... addtional arguments passed to the plot function.
#'
#' @return Plots the density distributions.
#'
#' @examples
#' gdensity(list(A=1:3,B=2:4,C=3:5),c("red","green","blue"))
#'
#' @export
gdensity <- function(thelist, cols="black", xlab="Value",ylab="Density",xlim=0,lwd=2, main="", ...){
  if(class(thelist)=="matrix"){
    thelist <- as.data.frame(thelist)
    names(thelist) <- paste("V",1:length(thelist),sep="")
  }
  minx <- F
  for(s in names(thelist)){
    d <- density(thelist[[s]],na.rm=T)
    if(class(minx)=='logical'){
      minx=min(d$x)
      maxx=max(d$x)
      miny=min(d$y)
      maxy=max(d$y)
    }else{
      if(min(d$x)<minx) minx<-min(d$x)
      if(min(d$y)<miny) miny<-min(d$y)
      if(max(d$x)>maxx) maxx<-max(d$x)
      if(max(d$y)>maxy) maxy<-max(d$y)
    }
  }
  if(!is.null(xlim)){
    if(length(xlim)==1){
        if(xlim>0) maxx <- xlim
    }else{
        minx <- xlim[1]
        maxx <- xlim[2]
    }
  }
  if(length(cols)!=length(thelist)){
    if(length(cols)==1){
        cols <- rep(cols, length(thelist))
    }else{
        warning("length(cols) is different from the number of items.")
    }
  }
  for(i in 1:length(thelist)){
    if(i==1){
      plot(density(thelist[[i]],na.rm=T),xlim=c(minx,maxx),ylim=c(miny,maxy),col=cols[i],xlab=xlab,ylab=ylab,lwd=lwd,main=main,...)
    }else{
      lines(density(thelist[[i]],na.rm=T),xlim=c(minx,maxx),ylim=c(miny,maxy),col=cols[i],lwd=lwd)
    }
  }
}
 

#' Make transparent
#'
#' Adds transparency to a color.
#'
#' @param x a color specification (see \code{\link[grDevices]{col2rgb}}
#' @param alpha a numeric value between 0 and 255 indicating the degree of transparency.
#'
#' @return A color code.
#'
#' @examples
#' maketrans("blue")
#'
#' @export
maketrans <- function(tcol,alpha=100){
  c <- col2rgb(tcol)
  rgb(c["red",1][[1]],c["green",1][[1]],c["blue",1][[1]],alpha,maxColorValue=255)
}


#' Histogram of missing values per peptide
#'
#' Plots the distribution of missing values per peptide in the desired channel/value.
#'
#' @param x A pSILAC object.
#' @param channel The channel in which to look for missing values. Either 'heavy', 'light', 'RIA' or 'hol' (H/L). 
#' RIA and hol should be equivalent, since both are missing when either the heavy or the light is missing.
#' @param a The starting value of the k_loss for model-fitting.
#'
#' @return A histogram.
#'
#' @export
plotMissingValues <- function(x,channel="RIA", ...){
  if(class(x) != "pSILAC")	stop("x should be a pSILAC object.")
  channel <- match.arg(channel, c("light","heavy","RIA","hol"))
  nna <- apply( x[[channel]], 1, FUN=function(y){ sum(is.na(y))})
  if(!("main" %in% names(list(...)))){
    main <- list(...)[["main"]]
  }else{
    main <- switch(channel, 
		    light="Missing values in the light channel",
		    heavy="Missing values in the heavy channel",
		    hol="Missing H/L values",
		    RIA="Missing RIA values")
  }
  return(hist(nna,main=main,xlab="Number of missing values per peptide",...))
}


#' plotProteinRIA
#'
#' Plots the relative isoform abundance for a given protein over time.
#'
#' @param object A pSILAC object.
#' @param protein A protein name (i.e. uniprot ID)
#' @param samples The samples to plot. If NULL, all samples are plotted.
#' @param peptides.alpha The alpha value of individual data points.
#' @param plot.legend Whether to plot a figure legend (default TRUE).
#' @param main The plot's title, passed to the plot function. Defaults to the protein's name/id.
#' @param xlab Passed to the plot function. Defaults to 'Time (hours)'.
#' @param ... additional arguments passed to the plot function. 
#'
#' @return Nothing, but plots a figure.
#'
#' @export
plotProteinRIA <- function(object, protein, samples=NULL, peptides.alpha=150,plot.legend=T,main=NULL,xlab="Time (hours)",...){
  if(class(object) != "pSILAC")	stop("'object' should be a pSILAC object.")
  if(is.null(object$RIA.kloss))	stop("Peptide-level and protein-level RIA rates of loss must first be calculated. See ?calcAllRates")
  
  if(is.null(samples)){
    
    samples <- as.character(unique(object$design$sample))
  
    }else{
    
    if(!all(samples %in% as.character(object$design$sample)))	stop("Unknown samples specified.")
    
      samples <- unique(object$design$sample[which(object$design$sample %in% samples)])
  }
  
  w <- which(object$design$sample %in% samples)
  
  pria <- getPeptides(object, protein=protein, returnValues="RIA")[,w]
  
  if(is.null(main))	main <- protein
  
  plot(rep(object$design$time[w],nrow(pria)),as.numeric(t(pria)),col=sapply(rep(object$design$color,nrow(pria)),alpha=peptides.alpha,FUN=maketrans),pch=20,cex=1.5,xlim=c(0,max(object$design$time,na.rm=T)),ylim=c(min(pria,na.rm=T),1),main=main,xlab=xlab,ylab="RIA",...)
  
  if(protein %in% row.names(object$protein.kloss)){
    f <- function(x,a){ exp(-a*x) }
    for(s in samples){
      curve(f(x,
            object$protein.kloss[row.names(object$protein.kloss) == protein, paste(s,"kloss",sep=".")]),
            lwd=2,
            col=unique(object$design$color[which(object$design$sample==s)]),
            add=TRUE)
    }
  }
  if(plot.legend)	legend("bottomleft",fill=object$design$color[!duplicated(object$design$sample) & object$design$sample %in% samples],legend=samples,bty="n")
}

#' plotPeptideRIA
#'
#' Plots the relative isoform abundance for a given peptide over time.
#'
#' @param object A pSILAC object.
#' @param peptide A peptide name (i.e. row names of the peptide table)
#' @param samples The samples to plot. If NULL, all samples are plotted.
#' @param peptides.alpha The alpha value of individual data points.
#' @param plot.legend Whether to plot a figure legend (default TRUE).
#' @param main The plot's title, passed to the plot function. Defaults to the protein's name/id.
#' @param xlab Passed to the plot function. Defaults to 'Time (hours)'.
#' @param ... additional arguments passed to the plot function. 
#'
#' @return Nothing, but plots a figure.
#'
#' @export
plotPeptideRIA <- function(object, peptide, samples=NULL, peptides.alpha=150,plot.legend=T,main=NULL,xlab="Time (hours)",...){
  if(class(object) != "pSILAC")	stop("'object' should be a pSILAC object.")
  if(is.null(object$RIA.kloss))	stop("Peptide-level RIA rates of loss must first be calculated. See ?calcAllRates")
  if(is.null(samples)){
    samples <- as.character(unique(object$design$sample))
  }else{
    if(!all(samples %in% as.character(object$design$sample)))	stop("Unknown samples specified.")
    samples <- unique(object$design$sample[which(object$design$sample %in% samples)])
  }
  w <- which(object$design$sample %in% samples)
  pria <- getPeptides(object, peptides=peptide, returnValues="RIA")[,w]
  if(is.null(main))	main <- peptide
  plot(rep(object$design$time[w],nrow(pria)),as.numeric(t(pria)),col=sapply(rep(object$design$color,nrow(pria)),alpha=peptides.alpha,FUN=maketrans),type="b",lty="dashed",cex=1,xlim=c(0,max(object$design$time,na.rm=T)),ylim=c(min(pria,na.rm=T),1),main=main,xlab=xlab,ylab="RIA",...)
  if(peptide %in% row.names(object$RIA.kloss)){
    f <- function(x,a){ exp(-a*x) }
    for(s in samples){
      if(is.na(object$RIA.kloss[peptide,paste(s,"kloss",sep=".")])){next}
      curve(f(x,object$RIA.kloss[peptide,paste(s,"kloss",sep=".")]),lwd=2,col=unique(object$design$color[which(object$design$sample==s)]),add=T)
    }
  }
  if(plot.legend)	legend("bottomleft",fill=object$design$color[!duplicated(object$design$sample) & object$design$sample %in% samples],legend=samples,bty="n")
}

#' plotPeptideHoL
#'
#' Plots the ln(H/L+1) for a given peptide over time.
#'
#' @param object A pSILAC object.
#' @param peptide A peptide name (i.e. row names of the peptide table)
#' @param samples The samples to plot. If NULL, all samples are plotted.
#' @param peptides.alpha The alpha value of individual data points.
#' @param plot.legend Whether to plot a figure legend (default TRUE).
#' @param main The plot's title, passed to the plot function. Defaults to the protein's name/id.
#' @param xlab Passed to the plot function. Defaults to 'Time (hours)'.
#' @param ... additional arguments passed to the plot function. 
#'
#' @return Nothing, but plots a figure.
#'
#' @export
plotPeptideHoL <- function(object, peptide, samples=NULL, peptides.alpha=150,plot.legend=T,main=NULL,xlab="Time (hours)",...){
  if(class(object) != "pSILAC")	stop("'object' should be a pSILAC object.")
  if(is.null(object$hol.kloss))	stop("Peptide-level and protein-level H/L-based rates of loss must first be calculated. See ?calcAllRates")
  if(is.null(samples)){
    samples <- as.character(unique(object$design$sample))
  }else{
    if(!all(samples %in% as.character(object$design$sample)))	stop("Unknown samples specified.")
    samples <- unique(object$design$sample[which(object$design$sample %in% samples)])
  }
  w <- which(object$design$sample %in% samples)
  phol <- getPeptides(object, peptides=peptide, returnValues="hol")[,w]
  design <- object$design[w,]
  if(is.null(main))	main <- peptide
  plot(rep(object$design$time[w],nrow(phol)),as.numeric(t(phol)),col=sapply(rep(design$color,nrow(phol)),alpha=peptides.alpha,FUN=maketrans),type="b",lty="dashed",cex=1,xlim=c(0,max(design$time,na.rm=T)),ylim=c(min(phol,na.rm=T),1),main=main,xlab=xlab,ylab="ln(H/L+1)",...)
  if(peptide %in% row.names(object$hol.kloss)){
    for(s in samples){
      if(is.na(object$hol.kloss[peptide,paste(s,"kloss",sep=".")])){next}
      abline(a=0,b=object$hol.kloss[peptide,paste(s,"kloss",sep=".")],lwd=2,col=unique(design$color[which(design$sample==s)]))
    }
  }
  if(plot.legend)	legend("topleft",fill=object$design$color[!duplicated(object$design$sample) & object$design$sample %in% samples],legend=samples,bty="n")
}



#' plotProteinHol
#'
#' Plots the ln(H/L+1) ratio for a given protein over time.
#'
#' @param object A pSILAC object.
#' @param protein A protein name (i.e. uniprot ID)
#' @param samples The samples to plot. If NULL, all samples are plotted.
#' @param peptides.alpha The alpha value of individual data points.
#' @param plot.legend Whether to plot a figure legend (default TRUE).
#' @param main The plot's title, passed to the plot function. Defaults to the protein's name/id.
#' @param xlab Passed to the plot function. Defaults to 'Time (hours)'.
#' @param ... additional arguments passed to the plot function. 
#'
#' @return Nothing, but plots a figure.
#'
#' @export
plotProteinHol <- function(object, protein, samples=NULL, peptides.alpha=150,plot.legend=T,main=NULL,xlab="Time (hours)",...){
  if(class(object) != "pSILAC")	stop("'object' should be a pSILAC object.")
  if(is.null(object$hol.kloss))	stop("Peptide-level and protein-level H/L-based rates of loss must first be calculated. See ?calcAllRates")
  if(is.null(samples)){
    samples <- as.character(unique(object$design$sample))
  }else{
    if(!all(samples %in% as.character(object$design$sample)))	stop("Unknown samples specified.")
    samples <- unique(object$design$sample[which(object$design$sample %in% samples)])    
  }
  w <- which(object$design$sample %in% samples)
  phol <- getPeptides(object, protein=protein, returnValues="hol")[,w]
  design <- object$design[w,]
  if(is.null(main))	main <- protein
  plot(rep(object$design$time[w],nrow(phol)),as.numeric(t(phol)),col=sapply(rep(object$design$color,nrow(phol)),alpha=peptides.alpha,FUN=maketrans),pch=20,cex=1.5,xlim=c(0,max(object$design$time,na.rm=T)),ylim=c(0,max(phol,na.rm=T)),main=main,xlab=xlab,ylab="ln(H/L+1)")
  if(protein %in% row.names(object$protein.kloss)){
    for(s in samples){
      abline(a=0,b=object$protein.kloss[row.names(object$protein.kloss) == protein, paste(s,"kloss",sep=".")],
             lwd=2,
             col=unique(design$color[which(design$sample==s)]))
    }
  }
  if(plot.legend)	legend("topleft",fill=object$design$color[!duplicated(object$design$sample) & object$design$sample %in% samples],legend=samples,bty="n")
}

#' Plots a heatmap of a protein's peptide-level kloss values.
#'
#' @param object A pSILAC object.
#' @param protein A protein name (i.e. uniprot ID)
#' @param method The method from which kloss values should be taken. Either 'RIA' or 'hol'.
#'
#' @return Nothing, but plots a figure.
#'
#' @export
kloss.heatmap <- function(object,protein,method="RIA"){
  if(class(object) != "pSILAC")	stop("'object' should be a pSILAC object.")
  method <- match.arg(method, c("RIA","hol"))
  library(gplots)
  m2 <- paste(method,"kloss",sep=".")
  if(is.null(object[[m2]]))	stop("Peptide-level rates of loss must first be calculated. See ?calcAllRates")
  pd <- getPeptides(object, protein=protein, returnValues="hol")[,paste(unique(object$design$sample),"kloss",sep=".")]
  colnames(pd) <- gsub(".kloss","",colnames(pd),fixed=T)
  pvar <- apply(object[[m2]][p,grep("stderr",colnames(object[[m2]]),fixed=T)],1,na.rm=T,FUN=mean)
  heatmap.2(	as.matrix(pd),
		RowSideColors=plColorMap(pvar),
		ColSideColors=as.character(object$design$color[!duplicated(object$design$sample)]),
		scale="none",trace="none",key=F,Rowv=NA,Colv=NA,dendrogram="none",
		cellnote=round(pd,3),notecol="black",
		margin=c(6,10),
		main=method,
		cexRow=1,cexCol=1.2,
  lmat=rbind(c(6,0,5),c(0,0,2),c(4,1,3)), 
  lhei=c(0.01,0.3,5.2),
  lwid=c(0.01,0.3,5.2)
  )
}

#' Plots a boxplot of a protein's peptide-level kloss values.
#'
#' @param object A pSILAC object.
#' @param protein A protein name (i.e. uniprot ID)
#' @param method The method from which kloss values should be taken. Either 'RIA' or 'hol'.
#' @param main The plot's title, passed to the plot function. Defaults to the protein's name/id.
#' @param in.all Starting with how many peptides should peptides not quantified in all samples be ignored (default 2), or FALSE to use all peptides, even if not quantified in all samples. 
#' @param removeOutliers Starting with how many peptides should outlier peptides be removed, or FALSE to disable outlier removal (default 5). See ?removeOutlierPeptides
#'
#' @return Nothing, but plots a figure.
#'
#' @export
kloss.boxplot <- function(object,protein,method="RIA",main=NULL,in.all=NULL,removeOutliers=NULL,normalize=FALSE){
  if(class(object) != "pSILAC")	stop("'object' should be a pSILAC object.")
  library(gplots)
  m2 <- paste(method,"kloss",sep=".")
  if(is.null(object[[m2]]))	stop("Peptide-level rates of loss must first be calculated. See ?calcAllRates")
  if(normalize){
    p <- intersect(row.names(getPeptides(object,protein=protein)),row.names(object[[m2]]))
    pd <- medianNorm(object[[m2]][,paste(unique(object$design$sample),"kloss",sep=".")])[p,]
  }else{
    pd <- getPeptides(object,protein=protein,returnValues=m2)[,paste(unique(object$design$sample),"kloss",sep="."), drop = FALSE]
  }
  pd <- pd[which(apply(pd,1,FUN=function(x){ !all(is.na(x)) })), , drop = FALSE]
  colnames(pd) <- gsub(".kloss$","",colnames(pd))
  if(is.null(in.all))	in.all <- .getProtKParam(object, "in.all")
  if(in.all){
    nna <- apply(pd,1,FUN=function(x){ !any(is.na(x)) })
    if(sum(nna)>0 | in.all == 1){
	pd <- pd[which(nna), , drop = FALSE]
    }
  }
  if(is.null(removeOutliers))	removeOutliers <- .getProtKParam(object, "removeOutliers")
  if(removeOutliers > 2 & nrow(pd) > 2)	pd <- removeOutlierPeptides(pd, min.size = as.numeric(removeOutliers))
  ll <- list()
  for(s in unique(object$design$sample)){
    ll[[s]] <- pd[,s]
  }
  boxplot(ll,col=as.character(object$design$color[!duplicated(object$design$sample)]),main=main, ylab=paste("Rate of loss (",ifelse(method=="hol","H/L","RIA"),"-based)",sep=""),las=3, outline=F)
}

#' plotProtein 
#'
#' Creates a layout of 4 plots of a protein's degradation: 
#'  - the protein's relative isotope abundance over time (calls 'plotProteinRIA'), 
#'  - a boxplot of the protein's peptide-level RIA-based kloss (calls 'kloss.boxplot'), 
#'  - the protein's ln(H/L+1) ratio over time (calls 'plotProteinHol'), 
#'  - a boxplot of the protein's peptide-level H/L-based kloss (calls 'kloss.boxplot'), 
#'
#' @param o A pSILAC object.
#' @param protein A protein name (i.e. uniprot ID)
#' @param removeOutliers Starting with how many peptides should outlier peptides be removed, or FALSE to disable outlier removal (default 5). See ?removeOutlierPeptides
#'
#' @return Nothing, but plots a figure.
#'
#' @export
plotProtein <- function(o, protein, removeOutliers=5){
  if(class(o) != "pSILAC")	stop("'object' should be a pSILAC object.")
  if(length(protein)>1)	stop("Please input a unique protein")
  if(is.null(o$protein.kloss)) stop("Protein-level rates of loss must first be calculated. See ?calcAllRates")
  if(!(protein %in% row.names(o$protein.kloss)))	stop("Protein not found!")
  layout(matrix(1:4,nrow=2))
  par(cex=1)
  plotProteinRIA(o,protein)
  plotProteinHol(o,protein)
  kloss.boxplot(o,protein,"RIA",paste(protein,"peptides"),removeOutliers=removeOutliers)
  kloss.boxplot(o,protein,"hol",paste(protein,"peptides"),removeOutliers=removeOutliers)
}

#' Plots each sample's distribution of protein kloss rates.
#'
#' @param o A pSILAC object.
#' @param remove.NAs Whether to remove proteins whose kloss is missing in any of the samples. Otherwise, all proteins will be used for the samples in which they are present.
#' @param xlim Either a single value (max x), or a vector of length 2 passed to the 'gdensity' function. Default c(0,0.1)
#' @param main Plot title, passed to the plot function.
#' @param xlab X-axis label, passed to the plot function.
#'
#' @return Nothing, but plots a figure.
#'
#' @export
plotKlossDensities <- function(o,remove.NAs=T,xlim=c(0,0.1),main="Protein-level rates of loss",xlab="K_loss"){
  if(class(o) != "pSILAC")	stop("'object' should be a pSILAC object.")
  if(is.null(o$protein.kloss)) stop("Protein-level rates of loss must first be calculated. See ?calcAllRates")
  if(remove.NAs){
    k <- o$protein.kloss[which(apply(o$protein.kloss,1,FUN=function(x){ !any(is.na(x)) })),]
  }else{
    k <- o$protein.kloss
  }
  layout(matrix(1:2,nrow=1))
  gdensity(k, o$design$color[!duplicated(o$design$sample)], main=main,xlab=xlab,xlim=xlim)
  legend("topright",bty="n",fill=o$design$color[!duplicated(o$design$sample)], legend=unique(o$design$sample))
  boxplot(log(k),col=o$design$color[!duplicated(o$design$sample)],ylab="ln(Rate of loss)",notch=T,outline=F,las=3)
  abline(h=median(apply(log(k),2,na.rm=T,FUN=median)),lty="dashed",col="grey")
}

#' Plot Intensity Distribution of pSILAC Channels
#'
#' This function plots the log-transformed intensity distribution for a specified channel in a pSILAC object.
#' It allows visualization of intensity distributions across different conditions, such as time points.
#'
#' @param o A `pSILAC` object containing intensity data.
#' @param channel A character string specifying the channel to plot. Options are `"heavy"`, `"light"`, `"sum"`, or `"NLI"`.
#'   Default is `"light"`.
#' @param ... Additional arguments passed to `gdensity` for customizing the plot.
#'
#' @details
#' The function uses the specified channel's intensity data from the pSILAC object `o`. If `channel = "sum"`, the sum of the
#' `heavy` and `light` channels is plotted. If `channel = "NLI"` and `NLI` data is not already in `o`, it will be calculated 
#' using `normalizeLightChannel`. The plot colors are determined by the time points in `o$design$time`.
#'
#' @examples
#' plotIntensityDistribution(o, channel = "light")
#'
#' @export
plotIntensityDistribution <- function(o, channel = "light", ...) {
  if (class(o) != "pSILAC") stop("'o' should be a pSILAC object.")
  channel <- match.arg(channel, c("heavy", "light", "sum", "NLI"))
  
  # Select data based on channel
  dat <- switch(channel,
                light = as.data.frame(o$light),
                heavy = as.data.frame(o$heavy),
                sum   = as.data.frame(.getChannelSum(o)),
                NLI   = as.data.frame(o$NLI)
  )
  
  # Compute NLI if necessary
  if (is.null(dat)) dat <- normalizeLightChannel(o)$NLI
  
  # Set colors for the plot based on time points
  cols <- plColorMap(o$design$time)
  
  # Plot intensity distribution
  gdensity(log(dat), xlab = "log(intensity)", cols = cols, ...)
  
  # Add legend
  legend("topright", bty = "n", fill = unique(cols), legend = paste(unique(o$design$time), "hours"))
}
