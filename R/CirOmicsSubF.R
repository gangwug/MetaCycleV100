###======================================================================================================================================
### This files contains functions used in 'CirOmicsMainF.R'
### Author: Gang Wu
### Email: wggucas@gmail.com
### Lab: John Hogenesch's lab in Perelman School of Medicine at University of Pennsylvania (http://hogeneschlab.org/)
###======================================================================================================================================
##adjusting output phase values from ARS with defined period length, and also extracting the phase and period value corresponding to the max amplitude for those profiles with two or more period and phase value reported by ARS.
adjPhaARS <- function(arsM,adjustV)
{
	## 'arsM' is a matrix containing all the output values from ARS, adjustV is the defined period length.
    if (adjustV == 0) {
		return(arsM[,c("period","phase","amplitude")]);                                            ##if the 'adjustV' is set '0', output the original period, phase and amplitude value from ARS
	} else {
		arsid <- dimnames(arsM)[[1]];
		arsnum <- as.numeric(arsM[,"period_number"]);
		names(arsnum) <- arsid;
		arsnum[is.na(arsnum) | is.nan(arsnum)] <- 0;
		uni_id <- names(arsnum[arsnum <= 1]);
		multi_id <- names(arsnum[arsnum > 1]);
        ## selecting the period and phase information corresponding to the largest amplitude if multiple periods for one probeset are exist.
		if (length(multi_id) > 0) {
            catM <- arsM[multi_id,c("period","phase","amplitude")];
            if (length(multi_id) == 1)
            { catM <- matrix(catM, nrow=1, ncol=length(catM));  }
            outM <- apply(catM,1,
						  function(z){
							zper <- unlist(strsplit(z[1],","));
							names(zper) <- 1:length(zper);
							zpha <- unlist(strsplit(z[2],","));
							names(zpha) <- 1:length(zpha);
							zamp <- unlist(strsplit(z[3],","));
							zamp <- as.numeric(zamp);
							names(zamp) <- 1:length(zamp);
							zamp<- sort(zamp,decreasing=TRUE);
							maxid<- names(zamp[1]);
							return(c(zper[maxid],zpha[maxid],zamp[maxid]));                        ##if one profile have two or more phase or period values, return only one phase and period corresponding to the max amplitude.
						  });
			outM <- t(outM);                                                                        ##each column corresponding to one profile in the results from 'apply(,1,)' 
			outM <- rbind(outM,arsM[uni_id,c("period","phase","amplitude")]);
			dimnames(outM) <- list("r"=c(multi_id,uni_id),"c"=c("period","phase","amplitude"));
			outM <- outM[arsid,];
		} else {
			outM <- arsM[arsid,c("period","phase","amplitude")];
		}
		#----------------adjusting phase value with defined period length
		adj_outM <- outM;
		pha <- as.numeric(adj_outM[,"phase"]);
		per <- as.numeric(adj_outM[,"period"]);
		adj_pha <- subAdjPha(pha, subper=per, adjustV);
		adj_outM <- cbind(adj_outM[,"period"],adj_pha,adj_outM[,"amplitude"]);
		dimnames(adj_outM) <- list("r"=arsid,"c"=c("period","phase","amplitude"));
		return(adj_outM);
	}
}
##-------------------------------------
##adjusting output phase values from JTK with defined period length.
adjPhaJTK <- function(jtkM,adjustV)
{
	## 'jtkM' is a matrix containing all the output values from JTK, adjustV is the defined period length.
    jtkid <- dimnames(jtkM)[[1]];
	pha <- as.numeric(jtkM[,"LAG"]);
    if (adjustV == 0) {
        adjpha <- pha;
    }  else  {
        pha <- pha + START_TIME;
        per <- as.numeric(jtkM[,"PER"]);
        adjpha <- subAdjPha(pha, subper=per, adjustV);
    }
    names(adjpha) <- jtkid;
	return(adjpha);
}
##-------------------------------------
##adjusting output phase value from LS with defined period length.
adjPhaLS <- function(lsM,adjustV)
{
	## 'lsM' is a matrix containing all the output values from LS, adjustV is the defined period length.
    lsid <- dimnames(lsM)[[1]];
	pha <- as.numeric(lsM[,"PhaseShift"]);
    if (adjustV == 0) {
        adjpha <- pha;
    }  else  {
        per <- as.numeric(lsM[,"Period"]);
        adjpha <- subAdjPha(pha, subper=per, adjustV);
    }
    names(adjpha) <- lsid;
	return(adjpha);
}
###======================================================================================================================================
##outputting the integrated P-value from multiple P-values calculated by multiple methods by "bonferroni" or "fisher method".
getCirOmicsPVA <- function(pvalueM, method="bonferroni")
{
    ##the first column in pvalueM is id name, and other columns are p-values
    pvalueID <- pvalueM[,1];
    pvaM <- matrix(rep(NA,nrow(pvalueM)),ncol=1);
    if ( (method == "bonferroni") | (method == "Bonferroni") ) {
		pvalue_meta <- NA;
		pvalue_meta <- apply(pvalueM[,-1],1,function(z) {
							z <- as.numeric(z);
                            z[is.na(z) | is.nan(z)] <- 1;
							z_meta <- p.adjust(z, method="bonferroni");
							z <- min(z_meta);
							return (z);
						   } )
		qvalue <- p.adjust(pvalue_meta, method="BH");
		pvaM <- cbind(pvalue_meta, qvalue);
	} else if ( (method == "Fisher") | (method == "fisher") ) {
		charM <- pvalueM[,-1];
		char2M <- matrix(as.numeric(charM),nrow=nrow(charM),ncol=ncol(charM));
		char2M[is.na(char2M) | is.nan(char2M)] <- 1;
		pvalue_meta <- fisher.method(char2M, p.corr="BH", zero.sub = 1e-50);
		pvalue_meta <- as.numeric(pvalue_meta$p.value);
		qvalue <- p.adjust(pvalue_meta, method="BH");
		pvaM <- cbind(pvalue_meta, qvalue);
	} else {
		cat("If the method is correlctly set as 'bonferroni', 'Bonferroni', 'bonf', 'fisher' or 'Fisher', there is unknown bug in combining P-values. Please contat the author of CirOmics. Thanks. \n");
	}
	rownames(pvaM)<- pvalueID;
	return (pvaM);
}
##-------------------------------------Period and phase average
##------------------------calculate the average period and phase
##outputting the integrated period and phase value from multiple period and phase values calculated by multiple methods.
getCirOmicsPerPha <- function(periodM, phaseM, pvalueM, adjustV, weightedPerPha)
{
    ##the first column of periodM, phaseM and pvalueM is id name, and other columns are period, phase and p-value information
    ##------------preparing weight matrix
    perpha_ID <- periodM[,1];
    rankM <- apply(pvalueM[perpha_ID,-1], 2, function(z) {
                                 z <- as.numeric(z);
                                 z[is.na(z) | is.nan(z)] <- 1;
                                 z[!z] <- 1e-300;
                                 z <- -log10(z);
                                 return (rank(z));
                                 });
    perphaM <- cbind(periodM[,-1], phaseM[perpha_ID,-1], rankM);
    ##------------calculating average period length and phase value
    ADJ <<- adjustV;
    perpha_avgM <- apply(perphaM, 1, function(z) {
                        z <- as.numeric(z);
                        zper <- z[1:(length(z)/3)];
                        zpha <- z[(length(z)/3+1):(length(z)/3*2)];
                        zwei <- z[(length(z)/3*2+1):length(z)];
                        if (!weightedPerPha)
                        {   zwei <- rep(1, length(zwei));    }
                        per_index <- which(!is.na(zper) & !is.nan(zper));
                        zper_mean <- sum(zper[per_index]*zwei[per_index])/sum(zwei[per_index]);
                        pha_index <- which(!is.na(zpha) & !is.nan(zpha));
                        zpha_mean <- circularMean(zpha[pha_index], subper=zper[pha_index], zweit=zwei[pha_index], meanper=zper_mean, subadj=ADJ);
                        return(c(zper_mean, zpha_mean));
                     });
    perpha_avgM <- t(perpha_avgM);
    per_avg <- perpha_avgM[,1];
    pha_avg <- perpha_avgM[,2];
    ##adjusting the average phase
    pha_avg <- subAdjPha(pha_avg, subper=per_avg, ADJ);
    ##------------
    outM <- cbind(per_avg, pha_avg);
    rownames(outM) <- perpha_ID;
    return(outM);
}
##-------------------------------------Amplitude calculation
##------------------------amplitude calculation method: ordinary least square method
##calculating the estimated baseline and relative amplitude based on the ordinary least square method.
##If Y is expression values corresponding to time points(tim),
##Thus set a general expression function as simple as: Y = Baseline + Trend*tim + Amplitude*cos(2*pi/Period*(tim - Phase)); Period and Phase can be taken from average period and phase calculated above, respectively.
##Using ordinary least square method to estimate three unknown parameters - Baseline, Trend and Amplitude; Baseline, Amplitude and Relative Amplitude will be outputted in the results of 'CirOmics'. The Relative Amplitude is defined as Amplitude/Baseline.
getCirOmicsAMP <- function(exprperphaM)
{
    outM <- apply(exprperphaM, 1, function(z) {
                    z <- as.numeric(z);
                    tim <- AMPTIM;                                                                 ##AMPTIM is a global variables defined in the script file of 'CirOmicsMainF.R'
                    trendt <- tim - mean(tim[!is.na(tim) & !is.nan(tim)]);
                    expr <- z[1:(length(z) - 2)];
                    per <- z[length(z) - 1];
                    pha <- z[length(z)];
                    if ( (is.na(per)) | (is.na(pha)) | is.nan(per) | is.nan(pha) )
                    {
                        exprv <- expr[!is.na(expr) & !is.nan(expr)];
                        basev <- median(exprv);
                        out <- c(basev, NA, NA);
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
                            if (length(exprv) > length(unique(timv)) )                             ##calculating the mean value if with replicates; if '3.1,2.9' are two replicate values in one time point, detrend_linear() will give different results when changing the order to '2.9, 3.1'
                            {
                                timv_uni <- sort(unique(timv));
                                timv <- factor(timv, levels=timv_uni);
                                exprv_mean <- tapply(exprv, timv, mean);
                                exprv_mean <- as.numeric(exprv_mean);
                                timv <- timv_uni;
                                exprv <- exprv_mean;
                            }
                            exprv_dt <- detrend_linear(exprv);                                     ##detrend_linear() is in 'ARS.R' file, detrending the obvious trend in the expression profile
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
                    names(out) <- c("Baseline", "Amplitude", "RelativeAmplitude");
                    return(out);
                  });
    outM <- t(outM);                                                                               ##each column corresponding to one profile in the results from 'apply(,1,)' 
    return(outM);
}
###======================================================================================================================================
###End
