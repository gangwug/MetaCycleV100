###======================================================================================================================================
### Authors of original R code of Lomb-Scargle: Earl F. Glynn, Jie Chen and Arcady R. Mushegian
### Email: chenj@umkc.edu
### Associated literature: Earl F. Glynn, Jie Chen and Arcady R. Mushegian. Bioinformatics. 22(3):310-6 (2006).
### Website: http://research.stowers-institute.org/mcm/efg/2005/LombScargle/R/index.htm
### Lomb-Scargle Normalized Periodogram:  
    ## "Fast Algorithm for Spectral Analysis of Unevenly Sampled Data". William H. Press and George B. Rybicki. Astrophysical Journal, 338:277-280, March 1989.
    ## Also appeared in Section 13.8, "Spectral Analysis of Unevenly Sampled Data" in Numerical Recipes in C (2nd Ed). William H. Press, et al, Cambridge University Press, 1992.
### Revision of Lomb-Scargle script: Gang Wu
### Email: wggucas@gmail.com
### Lab: John Hogenesch's lab in Perelman School of Medicine at University of Pennsylvania (http://hogeneschlab.org/)
###======================================================================================================================================
###--------------------------------------------------amplitude calculation
###------------------------amplitude calculation method: ordinary least square method
###If Y is expression values corresponding to time points(tim),
###Thus set a general expression function as simple as: Y = Baseline + Trend*( tim - mean(tim) ) + Amplitude*cos(2*pi/Period*(tim - Phase)); Period and Phase can be taken from average period and phase calculated above, respectively.
###Using ordinary least square method to estimate three unknown parameters - Baseline, Trend and Amplitude; Baseline and Relative Amplitude will be outputed in the results of 'CirOmics'. The Relative Amplitude is definated as Amplitude/Baseline.
getAMP <- function(exprperpha, amptim)
{
    z <- exprperpha;
    tim <- amptim;                                                                 ##amptim is a variable defining the time points
    trendt <- tim - mean(tim[!is.na(tim) & !is.nan(tim)]);
    expr <- z[1:(length(z) - 2)];
    per <- z[length(z) - 1];
    pha <- z[length(z)];
    if ( (is.na(per)) | is.na(pha) | is.nan(per) | is.nan(pha) )
    {
        exprv <- expr[!is.na(expr) & !is.nan(expr)];
        basev <- median(exprv);
        out <- c(basev, NaN, NaN);
    }  else  {
        cost <- cos(2*pi/per*(tim - pha));
        fit <- lm(expr~trendt + cost);
        fitcoef <- fit$coefficients;
        basev <- fitcoef[1];
        ampv <- fitcoef[3];
        exprv <- expr[!is.na(expr) & !is.nan(expr)];
        if ( (basev < min(exprv)) | (basev > max(exprv)) ) {
            basev <- median(exprv);
            timv <- tim[!is.na(expr) & !is.nan(expr)];
            if (length(exprv) > length(unique(timv)) )                                             ##calculating the mean value if with replicates; if '3.1,2.9' are two replicate values in one time point, detrend_linear() will give different results when changing the order to '2.9, 3.1'
            {
                timv_uni <- sort(unique(timv));
                timv <- factor(timv, levels=timv_uni);
                exprv_mean <- tapply(exprv, timv, mean);
                exprv_mean <- as.numeric(exprv_mean);
                timv <- timv_uni;
                exprv <- exprv_mean;
            }
            exprv_dt <- detrend_linear(exprv);                                                    ##detrend_linear(source code in 'ARS.R') function is used to detrend the obvious trend in the expression profile
            costv <- cos(2*pi/per*(timv - pha));
            fitv <- lm(exprv_dt~costv);
            fitcoefv <- fitv$coefficients;
            ampv <- fitcoefv[2];
        } 
        if (abs(basev) >= 1) { 
            out <- c( basev, abs(ampv), abs(ampv/basev) );
        } else {
            out <- c( basev, abs(ampv), abs(ampv) );
        }
    }
    names(out) <- c("ExpBaseline", "Amplitude", "RelativeAmplitude");
    return(out);
}
###======================================================================================================================================
runCirCaseLS <- function(indata,LStime,minper=20,maxper=28)
{
	RawData <- indata;
    outID <- dimnames(RawData)[[1]];
    libID <- dimnames(RawData)[[2]];
	Expression <- data.matrix(RawData);
	Time <- LStime;                             
	###-----------------------
	N <- ncol(Expression);
	stopifnot(length(Time) == N);
	M <- 4*N;
	MinFrequency <- 1/maxper;
	MaxFrequency <- 1/minper;
	if (MaxFrequency > 1/(2*mean(diff(Time))))
	{  cat("MaxFrequency may be above Nyquist limit.\n"); }
	TestFrequencies <- MinFrequency +
					  (MaxFrequency - MinFrequency) * (0:(M-1) / (M-1))
	###-----------------------
	header <- c("CirID", "p", "Period", "PhaseShift", "PhaseShiftHeight", "CirCase_baselineEXP", "CirCase_AMP", "CirCase_relativeAMP");
	LSoutM <- header;
	for (j in 1:length(outID))                                                     
	{
		if (.Platform$OS.type == "windows")                                       
		{ flush.console(); }                                                       # display immediately in Windows 
		Nindependent <- NHorneBaliunas(sum(!is.na(Expression[j,])));               # NHorneBaliunas(source code in 'LS.R')
		LS <- ComputeAndPlotLombScargle(Time, Expression[j,], TestFrequencies, Nindependent);  #ComputeAndPlotLombScargle(source code in 'LS.R')
        lsamp <- getAMP(exprperpha=c(Expression[j,], LS$PeakPeriod, LS$h.peak$maximum), amptim=Time);
		LSoutM <- rbind(LSoutM, c(outID[j], LS$PeakPvalue, LS$PeakPeriod, LS$h.peak$maximum, LS$h.peak$objective, lsamp));
	}
	LSoutM <- LSoutM[2:nrow(LSoutM),];
    colnames(LSoutM) <- header;
	bhq <- p.adjust(as.numeric(LSoutM[,"p"]),"BH");
	LSoutM <- cbind(LSoutM,bhq,Expression);	       
	dimnames(LSoutM) <-list("r"=outID, "c"=c(header, "BH.Q", libID));
    LSoutM <- LSoutM[,c(1:2, 9, 3:8, 10:ncol(LSoutM))];
	return(LSoutM);
}
###============================================================================================  
