######======================================================================================================================================
###### This files contains functions used in 'CirCaseMainF.R'
###### Author: Gang Wu
###### Email: wggucas@gmail.com
###### Lab: John Hogenesch's lab in Perelman School of Medicine at University of Pennsylvania (http://hogeneschlab.org/)
######======================================================================================================================================
##extracting numeric number from strings
getNumber <- function(strvector)
{
    ###extracting numeric number from strings
    strvector <- as.character(strvector);
    strL <- strsplit(strvector, "");
    numvector <- sapply(strL, function(z) {
                            zindex <- grep("\\d", z);
                            if (length(zindex) > 1) {
                                zout <- paste(z[zindex[1]:zindex[length(zindex)]], collapse="");
                                return(as.numeric(zout));
                            }   else if (length(zindex) == 1)  {
                                return(as.numeric(z[zindex]));
                            }   else {
                                return(NA);
                            }
                        });
    return(numvector);
}
######---------------------------------------
##integrating multiple P-values by "fisher method".
getCirCasePVA <- function(pvalueM, method="Fisher")
{
    ### the first column in pvalueM is id name, and other columns are p-values
    pvalueID <- pvalueM[,1];
    pvaM <- matrix(rep(NA,nrow(pvalueM)),ncol=1);
    if ( (method == "Fisher") | (method == "fisher") ) {
		charM <- pvalueM[,-1];
		char2M <- matrix(as.numeric(charM),nrow=nrow(charM),ncol=ncol(charM));
		char2M[is.na(char2M) | is.nan(char2M)] <- 1;
		pvalue_meta <- fisher.method(char2M, p.corr="BH", zero.sub = 1e-50);
		pvalue_meta <- as.numeric(pvalue_meta$p.value);
		qvalue <- p.adjust(pvalue_meta, method="BH");
		pvaM <- cbind(pvalue_meta, qvalue);
	} else {
		cat("If the method is correctly set as 'Fisher' or 'fisher', there is unknown bug in combining P-values. Please contact the author of 'circase'. Thanks. \n");
	}
	rownames(pvaM) <- pvalueID;
	return(pvaM);
}
######---------------------------------------
##transferring p-value to weight value
getweightM <- function(pvalueM)
{
    weitID <- pvalueM[,1];
    weitM <- apply(pvalueM[,-1], 1, function(z) {
                                 z <- as.numeric(z);
                                 z[is.na(z) | is.nan(z)] <- 1;
                                 z[!z] <- 1e-300;
                                 z <- -log10(z);
                                 return(z);
                                 });
    weitM <- t(weitM);
    dimnames(weitM)[[1]] <- weitID;
    return(weitM);
}
######---------------------------------------
##averaging multiple period and phase values.
getCirCasePerPha <- function(periodM, phaseM, pvalueM, adjustV)
{
    ###the first column of periodM, phaseM and pvalueM is id name (id order), and other columns are period, phase and p-value information
    perpha_ID <- periodM[,1];
    weitM <- getweightM(pvalueM);
    perphaM <- cbind(periodM[,-1], phaseM[perpha_ID,-1], weitM[perpha_ID,]);
    ###------------
    ADJL <<- adjustV;
    perpha_avgM <- apply(perphaM, 1, function(z) {
                        z <- as.numeric(z);
                        zper <- z[1:(length(z)/3)];
                        zpha <- z[(length(z)/3+1):(length(z)/3*2)];
                        zwei <- z[(length(z)/3*2+1):length(z)];
                        if (!WEIT)                                                                 ##WEIT is a global variable, which is defined by 'weightS' in 'integrateLSOUT'.
                        {   zwei <- rep(1, length(zwei));    }
                        per_index <- which(!is.na(zper) & !is.nan(zper));
                        zper_mean <- sum(zper[per_index]*zwei[per_index])/sum(zwei[per_index]);
                        pha_index <- which(!is.na(zpha) & !is.nan(zpha));
                        zpha_mean <- circularMean(zpha[pha_index], subper=zper[pha_index], zweit=zwei[pha_index], meanper=zper_mean, subadj=ADJL);
                        return(c(zper_mean,zpha_mean));
                     } );
    perpha_avgM <- t(perpha_avgM);
    per_avg <- perpha_avgM[,1];
    pha_avg <- perpha_avgM[,2];
    pha_avg <- subAdjPha(pha_avg, subper=per_avg, adjustV=ADJL);                                    ###adjusting the average phase
    ###------------
    outM <- cbind(per_avg, pha_avg);
    rownames(outM) <- perpha_ID;
    return(outM);
}
######-------------------------------------
##baseline and relative amplitude calculation
getCirCaseBaseAMP <- function(baseM, ampM, rampM, pvalueM)
{
    baseID <- baseM[,1];
    weitM <- getweightM(pvalueM);
    baseweitM <- cbind(baseM[,-1], ampM[baseID, -1], rampM[baseID,-1], weitM[baseID,]);
    base_avgM <- apply(baseweitM, 1, function(z) {
                                z <- as.numeric(z);
                                zbase <- z[1:(length(z)/4)];
                                zamp <- z[(length(z)/4 + 1): (length(z)/4*2)];
                                zramp <- z[(length(z)/4*2 + 1): (length(z)/4*3)];
                                zwei <- z[(length(z)/4*3 + 1):length(z)];
                                if (!WEIT)                                                         ##WEIT is a global variable, which is defined by 'weightS' in 'integrateLSOUT'.
                                {   zwei <- rep(1, length(zwei));    }
                                base_index <- which(!is.na(zbase) & !is.nan(zbase));
                                base_mean <- sum(zbase[base_index]*zwei[base_index])/sum(zwei[base_index]);
                                amp_index <- which(!is.na(zamp) & !is.nan(zamp));
                                amp_mean <- sum(zamp[amp_index]*zwei[amp_index]) / sum(zwei[amp_index]);
                                ramp_index <- which(!is.na(zramp) & !is.nan(zramp));
                                ramp_mean <- sum(zramp[ramp_index]*zwei[ramp_index])/sum(zwei[ramp_index]);
                                return(c(base_mean, amp_mean, ramp_mean));
                            });
    base_avgM <- t(base_avgM);
    dimnames(base_avgM) <- list("r"=baseID,"c"=c("base", "amp", "ramp"));
    return(base_avgM);
}
######======================================================================================================================================
######integrating multiple LS output results, and outputting integrated P-value, period, phase, expression baseline and relative amplitude values
integrateLSOUT <- function(lsoutL, adjustV, weightS)
{
    generalID <- NA;
    lsLen <- length(lsoutL);
    WEIT <<- weightS;
    if (lsLen)
    {
        initM <- lsoutL[[1]];
        generalID <- initM[,"CirID"];
        outM <- initM[,c("CirID", "p", "BH.Q", "Period", "PhaseShift", "CirCase_baselineEXP", "CirCase_AMP", "CirCase_relativeAMP")];
        dimnames(outM)[[1]] <- generalID;
        dimnames(outM)[[2]] <- c("CirID", "CirCase_Pvalue", "CirCase_BH.Q", "CirCase_Period", "CirCase_Phase", "CirCase_baselineEXP", "CirCase_AMP", "CirCase_relativeAMP");
        if (lsLen >= 2)
        {
            integrateM <- initM[,c("p", "Period", "PhaseShift", "CirCase_baselineEXP", "CirCase_AMP", "CirCase_relativeAMP")];
            dimnames(integrateM)[[1]] <- generalID;
            flag=0;
            for (i in 2:lsLen)
            {
                processM <- lsoutL[[i]];
                if (all(generalID == processM[,"CirID"]))
                {
                    integrateM <- cbind(integrateM, processM[, c("p", "Period", "PhaseShift", "CirCase_baselineEXP", "CirCase_AMP", "CirCase_relativeAMP")]);
                }  else  {
                    flag <- c(flag, i);
                }
            }
            if (length(flag) > 1)
            {   cat("Not all subjects have the same row names. Please check the input 'datafile'.\n");    }
            ###------------integrating p-value, period, phase, expression baseline and relative amplitude
            if (ncol(integrateM) > 6)
            {
                subjectNum <- ncol(integrateM)/6;
                pvalueM <- periodM <- phaseM <- baseM <- ampM <- rampM <- generalID;
                for (k in 1:subjectNum)
                {
                    pvalueM <- cbind(pvalueM, integrateM[,(k-1)*6+1]);
                    periodM <- cbind(periodM, integrateM[,(k-1)*6+2]);
                    phaseM <- cbind(phaseM, integrateM[,(k-1)*6+3]);
                    baseM <- cbind(baseM, integrateM[,(k-1)*6+4]);
                    ampM <- cbind(ampM, integrateM[,(k-1)*6+5]);
                    rampM <- cbind(rampM, integrateM[,(k-1)*6+6]);
                }
                dimnames(pvalueM)[[1]] <- dimnames(periodM)[[1]] <- dimnames(phaseM)[[1]] <- dimnames(baseM)[[1]] <- dimnames(ampM)[[1]] <- dimnames(rampM)[[1]] <- generalID;
                ###calling associated integration function
                cirhomoPvaM <- getCirCasePVA(pvalueM, method="Fisher");                                      ##'pvalue_meta', 'qvalue'
                cirhomoPerPhaM <- getCirCasePerPha(periodM, phaseM, pvalueM, adjustV);                       ##'per_avg', 'pha_avg'
                cirhomoBaseAmpM <- getCirCaseBaseAMP(baseM, ampM, rampM, pvalueM);                           ##'base', 'amp'
                outM[, c("CirCase_Pvalue", "CirCase_BH.Q")] <- cirhomoPvaM[generalID,];
                outM[, c("CirCase_Period", "CirCase_Phase")] <- cirhomoPerPhaM[generalID,];
                outM[, c("CirCase_baselineEXP", "CirCase_AMP", "CirCase_relativeAMP")] <- cirhomoBaseAmpM[generalID,];
            }  else  {
                cirhomoPha <- subAdjPha(subpha=integrateM[,3], subper=integrateM[,2], adjustV);
                outM[,"CirCase_Phase"] <- cirhomoPha;
            }
        }  else if  (lsLen == 1)  {
            cirhomoPha <- subAdjPha(subpha=initM[,"PhaseShift"], subper=initM[,"Period"], adjustV);
            outM[,"CirCase_Phase"] <- cirhomoPha;
        }  else  {
            cat("There is unknown bug associated with 'adjustedPhase'. Please contact the author. Thanks. \n");
        }
        return(outM);
    }  else  {
        cat("The input parameter of 'lsoutL' in 'integrateLSOUT' is empty, please check whether each parameter in 'circase' is correctly setted according to instruction. \n");
        return(NA);
    }
}
######======================================================================================================================================
######End
