###======================================================================================================================================
### Authors of original code of fisher's method: Karl Kugler <karl@eigenlab.net>, and Laurin AJ Mueller and Armin Graber
### Citation: MADAM - An Open Source Toolbox for Meta-Analysis. Source Code for Biology and Medicine 2010, 5:3
### The minor modification here is removing the code associated with 'multicore' package
###======================================================================================================================================
## function to calculate fisher sum
## p: vector of p-values
fisher.sum <- function(p, zero.sub=0.00001, na.rm=FALSE){
  if(any(p>1, na.rm=TRUE)||any(p<0, na.rm=TRUE))
    stop("You provided bad p-values")
  stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
  p[p==0] <- zero.sub
  if(na.rm)
    p<- p[!is.na(p)]
  S= -2*sum(log(p))
  res <- data.frame(S=S, num.p=length(p))
  return(res)
}
## main function of combining p-values by performing Fisher's method
fisher.method <- function(pvals, method=c("fisher"), p.corr=c("bonferroni","BH","none"), zero.sub=0.00001, na.rm=FALSE){
  stopifnot(method %in% c("fisher"))
  stopifnot(p.corr %in% c("none","bonferroni","BH"))
  stopifnot(all(pvals>=0, na.rm=TRUE) & all(pvals<=1, na.rm=TRUE))
  stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
  if(is.null(dim(pvals)))
    stop("pvals must have a dim attribute")
  p.corr <- ifelse(length(p.corr)!=1, "BH", p.corr)
  ##substitute p-values of 0
  pvals[pvals == 0] <- zero.sub
  fisher.sums <- data.frame(do.call(rbind, apply(pvals, 1, fisher.sum, zero.sub=zero.sub, na.rm=na.rm)))
  
  rownames(fisher.sums) <- rownames(pvals)
  fisher.sums$p.value <- 1-pchisq(fisher.sums$S, df=2*fisher.sums$num.p)
  fisher.sums$p.adj <- switch(p.corr,
                              bonferroni = p.adjust(fisher.sums$p.value, "bonferroni"),
                              BH = p.adjust(fisher.sums$p.value, "BH"),
                              none = fisher.sums$p.value)
  return(fisher.sums)
}
###======================================================================================================================================
###other functions shared by 'ciromics' and 'circase'
##adjusting phase with period length, eg. transfering phase from '26h' to '2h', if period length is 24h.
subAdjPha <- function(subpha, subper, adjustV)
{
    ###adjusting phase with expected period length
    if (adjustV == "PERIOD") {
        subper <- as.numeric(subper);
    } else if ( ( is.numeric(adjustV) ) & (adjustV > 0) ) {
        subper <- rep(adjustV, length(subpha));
    } else {
        return(as.numeric(subpha));
    }
    subpha <- as.numeric(subpha);
    pha_fold <- floor(subpha/subper);
    subadj_pha <- subpha - pha_fold*subper;
    return(subadj_pha);
}
##phase average method: circular average
##calculating the mean phase of multiple phases from different methods using the method of 'mean of circular quantities'.
circularMean <- function (z, subper, zweit, meanper, subadj)
{
	if ( length(z) > 0 )  {
        if (is.numeric(subadj))
        {
            if (subadj > 0)
            {   
                subper <- rep(subadj, length(subper)); 
            } else if (!subadj) {                                         ##subadj = 0 (indicating no adjustment of phase); can not use 'circularMean' method, and return linear mean values.
                zmean <- sum(z*zweit)/sum(zweit);
                return(zmean);
            } else {
                cat("There is unknown bug associated with 'adjustedPhase'. Please contact the author. Thanks. \n");
            }
        }
        zpolar <- z/subper*2*pi;                                          ##transform phase value to angle value; for 'predictedPer', 'subper' here is the period length predicted by each method; 
                                                                          ##another strategy may use the same 'meanper' value for different methods, currently it is not used based on the consideration that using predicted period length by each method seems more reasonable; 
                                                                          ##eg. if 'meanper'=25, the max per for JTK is 24 (one cycle sampling), it seems more reasonable to use predicted period by JTK for transferring its phase values to polar values than using the 'meanper'. 
        siny <- sum(sin(zpolar)*zweit);                                   ##do not have to divided by 'sum(zweit)', since for 'y/x' will get the same result if divided or not divided by the same value
        cosx <- sum(cos(zpolar)*zweit);
        meanpolar <- atan2(siny, cosx);                                   ##get mean of circular quantities
        meanpha <- meanpolar/(2*pi)*meanper;
        if (meanpha < 0)
        { meanpha <- meanpha +  meanper;  } 
        return(meanpha);
    }  else {
        return(NA);
    }
}
###======================================================================================================================================
