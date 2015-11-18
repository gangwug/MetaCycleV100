##' Detect rhythimic signals from time-series datasets with individual
##'   information
##'
##' This is a function that takes use of Lomb-Scargle to detect rhythmic
##'   signals from time-series datasets containing individual information.
##'
##' This function is originally aimed to analyze large scale perodic data with
##'   individual information. Please pay attention to the data format of
##'   \code{datafile} and \code{designfile}(see \code{Examples} part).
##'   Time-series experimental values(missing values as \code{NA}) from
##'   all individuals shoud be stored in \code{datafile}, with the first row
##'   containing all library ID(unique identification number for each sample)
##'   and the first column containing all detected molecular names(eg.
##'   transcript ID or gene name). The \code{designfile} should at least have
##'   three columns-library ID, subject ID and sampling time collumn.
##'   Experimental group information of each subject ID may be in another
##'   column. In addition, sampling time information may be stored in multiple
##'   columns instead of one column. For example, sampling time-"36 hours" may
##'   be recorded as "day 2"(sampling day column, \code{design_dayColm}) plus
##'   "12 hours"(sampling hour column, \code{design_hrColm}). The library ID
##'   in \code{datafile} and \code{designfile} should be same. If there are
##'   different characters between library ID in these two files, try
##'   \code{design_libIDrename} to keep them same.
##'
##'   \href{http://research.stowers-institute.org/efg/2005/LombScargle/}{
##'   Lomb-Scargle}(Glynn, 2006) is used to analyze time-series profiles
##'   individual by individual. The original code is modified and added with
##'   calculating baseline value and amplitude(using similar method
##'   as \code{\link{meta2d}}). Then p-values, period, phase, baseline
##'   value, amplitude and relative amplitude from multiple individuals are
##'   integrated group by group. P-values from different individuals are
##'   integrated with Fisher's method(\code{"fisher"}) (Fisher,1925;
##'   implementation code from \pkg{MADAM}). The integrated period, baseline
##'   value, amplitude and relative amplitude are arithmetic mean of multiple
##'   individuals, respectively. The integrated phase is
##'   \href{https://en.wikipedia.org/wiki/Mean_of_circular_quantities}{
##'   mean of circular quantities}(\code{adjustedPhase = "predictedPer"}) or
##'   a arithmetic mean(\code{adjustedPhase = "notAdjusted"}) of multiple
##'   individual phases. If \code{weightedMethod = TRUE} is selected, weighted
##'   scores(\code{-log10(p-values)}) will be taken into account in calculation
##'   of integrated period, phase, baseline, amplitude and relative amplitude.
##'
##'   For completly removing the potential problem of averaging phases with
##'   quite different period length(also mentioned in \code{\link{meta2d}}),
##'   setting \code{minper} and \code{maxper} to a same value may be the only
##'   known way. For short time-series profiles(eg. 10 time points or less),
##'   p-values given by Lomb-Scargle may be over conservative, which will also
##'   lead to conservative integrated p-values. In such case, selecting a
##'   proper p-value cut-off based on the p-value distribution is suggested.
##'
##' @param datafile a character string. The name of data file containing
##'   time-sereis experimental values of all individuals.
##' @param designfile a character string. The name of experimental design file,
##'   at least containing the library ID(column names of \code{datafile}),
##'   subject ID(the individual corresponding to each library ID), and
##'   sampling time information of each library ID.
##' @param outdir a character string. The name of directory used to store
##'   output files.
##' @param filestyle a character vector(length 1 or 3). The data format of
##'   input files, must be \code{"txt"}, or \code{"csv"}, or a character
##'   vector containing field separator character(\code{sep}), quoting
##'   character(\code{quote}), and the character used for decimal
##'   points(\code{dec}, for details see \code{\link[utils]{read.table}}).
##' @param design_libColm a numeric value. The order index(from left to right)
##'   of the column storing library ID in \code{designfile}.
##' @param design_subjectColm a numeric value. The order index(from left to
##'   right) of the column storing subject ID in \code{designfile}.
##' @param minper a numeric value. The minimum period length of interested
##'   rhythms. The default is \code{20} for circadian rhythms.
##' @param maxper a numeric value. The maximum period length of interested
##'   rhythms. The default is \code{28} for circadian rhythms.
##' @param timeUnit a character string. The basic time-unit, must be one of
##'   \code{"day"}, \code{"hour"}(default for circadian study),
##'   \code{"minute"}, or \code{"second"} depending on specific experimental
##'   design.
##' @param design_hrColm a numeric value. The order index(from left to right)
##'   of the column storing time point value-sampling hour information in
##'   \code{designfile}. If there is no such column in \code{designfile},
##'   set it as \code{NULL}.
##' @param design_dayColm a numeric value. The order index(from left to right)
##'   of the column storing time point value-sampling day information in
##'   \code{designfile}. If there is no such column in \code{designfile},
##'   set it as \code{NULL}(default).
##' @param design_minColm a numeric value. The order index(from left to right)
##'   of the column storing time point value-sampling minute information in
##'   \code{designfile}. If there is no such column in \code{designfile},
##'   set it as \code{NULL}(default).
##' @param design_secColm a numeric value. The order index(from left to right)
##'   of the column storing time point value-sampling second information in
##'   \code{designfile}. If there is no such column in \code{designfile},
##'   set it as \code{NULL}(default).
##' @param design_groupColm a numeric value. The order index(from left to
##'   right) of the column storing experimental group information of each
##'   individual in \code{designfile}. If there is no such column in
##'   \code{designfile}, set it as \code{NULL}(default) and take all
##'   individuals as one group.
##' @param design_libIDrename a character vector(length 2) containing  a
##'   matchable character string in each library ID of \code{designfile}, and
##'   a replacement character string. If it is not necessary to replace
##'   characters in library ID of \code{designfile}, set it as \code{NULL}(
##'   default).
##' @param adjustedPhase a character string. The method used to adjust each
##'   phase calculated by Lomb-Scargle before getting integrated phase,
##'   must be one of \code{"predictedPer"}(adjust phase with predicted
##'   period length) or \code{"notAdjusted"}(not adjust phase).
##' @param weightedMethod logical. If \code{TRUE}(default), weighted score
##'   based on p-value of each individual will be used to integrate period,
##'   phase and amplitude values of multiple individuals.
##' @param combinePvalue a character string. The method used to integrate
##'   p-values of multiple individuals, currently only \code{"fisher"}(Fisher's
##'   method) could be selected.
##' @param outIntegrationFile a character string. This parameter controls what
##'   kinds of analysis results will be outputted, must be one of \code{"both"}
##'   (default), \code{"onlyIntegration"}(only output integrated analysis
##'   results of each experimental group), or \code{"noIntegration"}(only
##'   output analysis results of each individual).
##' @param dayZeroBased logical. If \code{TRUE}, the first sampling day is
##'   recorded as day zero in the \code{designfile}.
##' @param outSymbol a character string. A common prefix exists in the names of
##'   output files.
##' @return
##' \code{meta3d} will write analysis results to \code{outdir} instead of
##'   returning them as objects. Output files with "LSresultSubjectID" in
##'   the file name are Lomb-Scargle analysis results for each individual.
##'   Files named with "CirCaseIntegrationGroupID" store integrated p-values,
##'   period, phase, baseline, amplitude and relative amplitude values
##'   from multiple individuals of each group and calculated FDR values based
##'   on integrated p-values.
##' @references
##' Glynn E. F., Chen J., and Mushegian  A. R. (2006). Detecting periodic
##'   patterns in unevenly spaced gene expression time series using
##'   Lomb-Scargle periodograms. \emph{Bioinformatics}, \bold{22(3)},
##'   310--316
##'
##' Fisher, R.A. (1925). \emph{Statistical methods for research workers}.
##'   Oliver and Boyd (Edinburgh).
##'
##' Kugler K. G., Mueller L.A., and Graber A. (2010). MADAM - an open source
##'   toolbox for meta-analysis. \emph{Source Code for Biology and Medicine},
##'   \bold{5}, 3.
##'
##' @examples
##' # write 'cycHumanBloodData' and 'cycHumanBloodDesign' into two 'csv' files
##' write.csv(cycHumanBloodData, file="cycHumanBloodData.csv",
##'   row.names=FALSE)
##' write.csv(cycHumanBloodDesign, file="cycHumanBloodDesign.csv",
##'   row.names=FALSE)
##'
##' # detect circadian transcripts in studied individuals
##' meta3d(datafile="cycHumanBloodData.csv",
##'   designfile="cycHumanBloodDesign.csv", outdir="example", filestyle="csv",
##'   design_libColm=1, design_subjectColm=2, design_hrColm=4,
##'   design_groupColm=3)
##' @export

meta3d <- function(datafile, designfile, outdir="cirout", filestyle,
                    design_libColm, design_subjectColm, minper=20, maxper=28,
                    timeUnit="hour", design_hrColm, design_dayColm=NULL,
                    design_minColm=NULL, design_secColm=NULL,
                    design_groupColm=NULL, design_libIDrename=NULL,
                    adjustedPhase="predictedPer", weightedMethod=TRUE,
                    combinePvalue="fisher", outIntegrationFile="both",
                    dayZeroBased=FALSE, outSymbol="")
{
    ####exist objects before running this function
    obst <- ls();
    run_start=proc.time();

    ####create the directory of storing output files, if it is not exist
    outdir2 <- unlist(strsplit(outdir,.Platform$file.sep));
    OUTDIR <<- paste(outdir2,collapse=.Platform$file.sep);
    if (! file.exists(OUTDIR) )
    { dir.create(OUTDIR); }

    ####extract the field separator character(FILE_SEP), the set of quoting
    ####characters(FILE_QUOTE), the character used for decimal points(FILE_DEC).
    FILE_SEP <<- "";
    FILE_QUOTE <<- "";
    FILE_QUOTE2 <<- FALSE;
    FILE_DEC <<- "";
    if (length(filestyle) == 1) {
        if (filestyle=="csv") {
            FILE_SEP <<- ",";
            FILE_QUOTE <<- "\"";
            FILE_QUOTE2 <<- TRUE;
            FILE_DEC <<- ".";
        } else if (filestyle=="txt") {
            FILE_SEP <<- "\t";
            FILE_QUOTE <<- "";
            FILE_DEC <<- ".";
        } else {
            stop(c("Please set 'filestyle' before running this function ",
                   "(the 'filesyle' could be set as 'txt' or 'csv').\n") );
        }
    } else if (length(filestyle) > 1) {
        if (length(filestyle) == 3) {
            FILE_SEP <<- filestyle[1];
            FILE_QUOTE <<- filestyle[2];
            FILE_QUOTE2 <<- TRUE;
            FILE_DEC <<- filestyle[3];
        } else {
            stop(c("Please set 'filestyle' before running this function ",
                   "(the 'filesyle' should be assigned a vector containing ",
                   "three characters, which corresponding to symbols used to ",
                   "separate columns, quote values and used for decimal points).\n") );
        }
    } else {
            stop(c("Please set 'filestyle' before running this function ",
                   "(the 'filesyle' could be set as 'txt' or 'csv', or could be ",
                   "assigned a vector containing three characters, ",
                   "which corresponding to symbols used to separate columns, ",
                   "quote values and used for decimal points).\n") );
    }

    ####extract the essential information in the designfile
    ##check the 'designfile'
    designD <- read.table(file=designfile, header=TRUE, sep=FILE_SEP, quote=FILE_QUOTE, dec=FILE_DEC, stringsAsFactors=FALSE);
    IDS <<- dimnames(designD)[[1]];
    if (ncol(designD) < 3)
    {
        stop(c("The data in 'designfile' is stored in a data frame with less ",
               "than three columns ('design_libColm', 'design_subjectColm', ",
               "'design_hrColm' are necessary), please check whether ",
               "'filestyle' is correctly set.\n") );
    }
    if (nrow(designD) == 0)
    {
        stop(c("The 'designfile' contains only one line which is taken as ",
               "column names, please check the input file.\n") );
    }
    ##extract time information
    timeT <- dayT <- hrT <- minT <- secT <- rep(0, length(IDS));
    ##hour information
    noteA <- c("The setted column number of 'design_hrColm' is larger ",
               "than the total column number of designfile.",
               "The total column number of designfile is ");
    noteB <- c(". The reason may be mistakenly setting ",
               "this parameter, or mistakenly setting 'filestyle'.\n");
    if (length(design_hrColm) > 0)
    {
        if ( is.numeric(design_hrColm) & (length(design_hrColm) == 1) & (design_hrColm > 0) )
        {
            if ( design_hrColm > ncol(designD) )
            {
                stop(c(noteA, ncol(designD), noteB) );
            }
            hrT <- designD[IDS, as.numeric(design_hrColm)];
        }  else  {
            stop(c("Please set 'design_hrColm' as a positive numeric value ",
                   "corresponding to the column number storing ",
                   "timepoints-hour information in the 'designfile'.\n") );
        }
    }
    ##day information
    if (length(design_dayColm) > 0)
    {
        if ( is.numeric(design_dayColm) & (length(design_dayColm) == 1) & (design_dayColm > 0) )
        {
            if (design_dayColm > ncol(designD) )
            {
                stop(c(noteA, ncol(designD), noteB) );
            }
            dayT <- designD[IDS, design_dayColm];
        }  else  {
            stop(c("Please set 'design_dayColm' as a positive numeric value ",
                   "corresponding to the column number storing timepoints-day ",
                   "information in the 'designfile'.\n") );
        }
    }
    ##minute information
    if (length(design_minColm) > 0)
    {
        if ( is.numeric(design_minColm) & (length(design_minColm) == 1) & ( design_minColm > 0) )
        {
            if (design_minColm > ncol(designD))
            {
                stop(c(noteA, ncol(designD), noteB) );
            }
            minT <- designD[IDS, design_minColm];
        }  else  {
            stop(c("Please set 'design_minColm' as a positive numeric value ",
                   "corresponding to the column number storing ",
                   "timepoints-minute information in the 'designfile'.\n") );
        }
    }
    ##second information
    if (length(design_secColm) > 0)
    {
        if ( is.numeric(design_secColm) & (length(design_secColm) == 1) & ( design_secColm > 0) )
        {
            if (design_secColm > ncol(designD))
            {
                stop(c(noteA, ncol(designD), noteB) );
            }
            secT <- designD[IDS, design_secColm];
        }  else  {
            stop(c("Please set 'design_secColm' as a positive numeric value ",
                   "corresponding to the column number storing ",
                   "timepoints-second information in the 'designfile'.\n") );
        }
    }
    ##the function of 'getNumber()' is in the file of 'CirCaseSubF.R'.
    if ( !all(is.numeric(dayT)) )
    {    dayT <- getNumber(dayT);   }
    if ( !all(is.numeric(hrT)) )
    {    hrT <- getNumber(hrT);   }
    if ( !all(is.numeric(minT)) )
    {    minT <- getNumber(minT);   }
    if ( !all(is.numeric(secT)) )
    {    secT <- getNumber(secT);   }

    ##transfer multiple columns into one time column,
    ##and the default time unit is 'hour'
    if ( (timeUnit != "day") & (timeUnit != "hour") & (timeUnit != "minute") & (timeUnit != "second") )
    {
        stop("Please set 'timeUnit' as 'day', 'hour', 'minute', or 'second'.\n");
    }
    if ( !is.logical(dayZeroBased) )
    {
        stop("The 'dayZeroBased' should be set as a logical value (TRUE or FALSE).\n");
    }
    unitT <- c(24, 1, 1/60, 1/3600);
    if (timeUnit == "day")
    {    unitT <- c( 1, 1/24, 1/(60*24), 1/(3600*24) );   }
    if (timeUnit == "minute")
    {   unitT <- c(24*60, 60, 1, 1/60);    }
    if (timeUnit == "second")
    {   unitT <- c(24*3600, 3600, 60, 1);    }
    if ( (!dayZeroBased) & (sum(dayT)) )
    { dayT <- dayT - 1;  }
    timeT <- cbind(dayT*unitT[1], hrT*unitT[2], minT*unitT[3], secT*unitT[4]);
    negti <- NULL;
    for (ti in 1:nrow(timeT))
    {
        teptm <- timeT[ti,];
        if (length(teptm[teptm < 0]) > 0)
        {    negti <- c(negti, ti);    }
    }
    if (length(negti) > 0)
    {
        stop(c("The below lines contains negative time information",
                paste(", NegLine", negti, sep=""), ", please modify time ",
                "information in these lines before running 'meta3d'.\n") );
    }
    timeT <- apply(timeT, 1, sum);

    ##extract libraryID, subjectID and modifying the libraryID
    ##if 'design_libIDrename' is not NULL.
    if ( is.numeric(design_libColm) & (length(design_libColm) == 1) & (design_libColm > 0) )
    {
        libID <- designD[IDS, design_libColm];
    }  else  {
        stop(c("Please set 'design_libColm' as a positive numeric value ",
               "corresponding to the column number storing libraryID in ",
               "the 'designfile'.\n") );
    }
    if (length(design_libIDrename))
    {
        if (length(design_libIDrename) == 2) {
            libID <- gsub(design_libIDrename[1], design_libIDrename[2], libID);
        }  else  {
            stop(c("For trying to replace some characters in the library names ",
                   "in design file, please set 'design_libIDrename' as a ",
                   "character vector containing two elements, the first element ",
                   "are characters before replacement, and the second element ",
                   "are characters after replacement.") );
        }
    }
    if ( is.numeric(design_subjectColm) & (length(design_subjectColm) == 1) & (design_subjectColm > 0) )
    {
        subjectID <- designD[IDS, design_subjectColm];
        subjectID <- gsub("\\s+", "", subjectID);
    }  else  {
        stop(c("Please set 'design_subjectColm' as a positive numeric value ",
               "corresponding to the column number storing subjectID in the ",
               "'designfile'.\n") );
    }
    if ( ( design_libColm > ncol(designD) ) | (design_subjectColm > ncol(designD) ) )
    {
        stop(c("The setted column number of 'design_libColm' or ",
               "'design_subjectColm' is larger than the total column number ",
               "of designfile. The total column number of designfile is ",
               ncol(designD), ". The reason may be mistakenly setting one or ",
               "both of these two parameters, or mistakenly setting the 'filestyle'.\n") );
    }

    ####divide subjectID into different groups if design_groupColm is not set as '0'.
    groupID <- NULL;
    UNI_GROUPID <<- NULL;
    GROUP_SUBJECTL <<- list();
    if (length(design_groupColm) > 0)
    {
        if ( is.numeric(design_groupColm) & (length(design_groupColm) == 1) & (design_groupColm > 0) ) {
            if (design_groupColm > ncol(designD))
            {
                stop(c("The setted column number of 'design_groupColm' is ",
                       "larger than the total column number of designfile. ",
                       "The total column number of designfile is ",
                       ncol(designD), ". The reason may be mistakenly setting ",
                       "this parameter, or mistakenly setting the 'filestyle'.\n") );
            }
            groupID <- designD[IDS, design_groupColm];
            groupID <- gsub("\\s+", "", groupID);
            subjectID <- paste(subjectID, groupID, sep="");
            names(groupID) <- subjectID;
            UNI_GROUPID <- unique(groupID);
            for (k in 1:length(UNI_GROUPID))
            {
                group_subjectID <- groupID[groupID == UNI_GROUPID[k]];
                group_subjectID <- unique(names(group_subjectID));
                GROUP_SUBJECTL[[k]] <- group_subjectID;
            }
        }  else  {
            stop(c("Please set 'design_groupColm' as a positive numeric value ",
                   "corresponding to the column number storing groupID in the ",
                   "'designfile'.\n") );
        }
    }

    ####store the extracted information into a dataframe
    if (sum(timeT) > 0)
    {
        if ( length(timeT) == nrow(designD) )  {
            if ( length(groupID) )
            {
                DESIGND <<- data.frame("libID"=libID, "subjectID"=subjectID, "groupID"=groupID, "timeT"=timeT, row.names=IDS, stringsAsFactors=FALSE);
            }  else  {
                DESIGND <<- data.frame("libID"=libID, "subjectID"=subjectID, "timeT"=timeT, row.names=IDS, stringsAsFactors=FALSE);
            }
        }  else {
            stop("There is an unknown bug in the script of 'meta3d', please contact the author.\n");
        }
    }  else if (sum(timeT) == 0) {
        stop(c("There is no time information in the design file ",
               "or has not set any time column before running 'meta3d'.\n") );
    }  else  {
        stop("There is an unknown bug in the script of 'meta3d', please contact the author.\n");
    }

    ####check the input parameters
    ##check 'minper' and 'maxper'
    if ( (length(minper) != 1) | (length(maxper) != 1) )
    {   stop("The 'minper' or 'maxper' contains more than one element.\n");  }
    if ( (!is.numeric(minper)) | (minper <= 0 ) )
    {	stop("The 'minper' should be a positive numeric value.\n");	}
    if ( (!is.numeric(maxper)) | (maxper <= 0 ) )
    {	stop("The 'maxper' should be a positive numeric value.\n");	}
    if (minper > maxper)
    {	stop("The 'minper' should not be larger than 'maxper'.\n");	}

    ##check whether output the integration results at the end of analysis
    INTEGRATION <<- "";
    if (outIntegrationFile == "both") {
        INTEGRATION <- "both";
    } else if (outIntegrationFile == "onlyIntegration") {
        INTEGRATION <- "onlyIntegration";
    } else if (outIntegrationFile == "noIntegration") {
        INTEGRATION <- "noIntegration";
    } else {
        stop(c("Please check the parameter of 'outIntegrationFile', it should ",
               "be set as 'both', 'onlyIntegration' or 'noIntegration'.\n") );
    }

    ##check the expected period length used to adjust phase
    ADPHAL <<- "";
    if (adjustedPhase == "predictedPer") {
        ADPHAL <- "PERIOD";
    } else if (adjustedPhase == "notAdjusted") {
        ADPHAL <- 0;
    } else {
        stop("The 'adjustedPhase' should be set as 'predictedPer' or 'notAdjusted'.\n");
    }

    ###check 'combinePvalue' and 'weightedMethod'
    if ( (combinePvalue != "Fisher") & (combinePvalue != "fisher") )
    {
        stop("The 'combinePvalue' should be set as 'Fisher' (equal to 'fisher').\n");
    }
    if ( !is.logical(weightedMethod) )
    {
        stop("The 'weightedMethod' should be set as a logical value (TRUE or FALSE).\n");
    }

    ####divide libraryID into different subgroups according to the subjectID
    UNI_SUBJECTID <<- unique(subjectID);
    SUBJECT_DGL <<- list();
    for (i in 1:length(UNI_SUBJECTID))
    {
        subject_index <- which(subjectID == UNI_SUBJECTID[i]);
        SUBJECT_DGL[[i]] <- DESIGND[subject_index,];
    }

    ####run LS subject by subject
    EXPD <<- read.table(file=datafile, header=TRUE, sep=FILE_SEP, quote=FILE_QUOTE, dec=FILE_DEC, stringsAsFactors=FALSE);
    EXP_IDS <<- dimnames(EXPD)[[1]];
    SIG_row <<- FALSE;
    if (ncol(EXPD) == 1)
    {
        stop(c("The whole input file is stored in a one-column dataframe, ",
               "please check whether 'filestyle' is correctly set.\n") );
    }
    if (nrow(EXPD) == 0)
    {
        stop(c("The whole input file contains only one line which is taken ",
               "as column names by 'meta3d', please check the input file.\n") );
    }
    if ( nrow(DESIGND) != ( ncol(EXPD) - 1) )
    {
        stop(c("Please check the 'designfile' and 'datafile', the column number ",
               "of data file should be one more than the number of total samples ",
               "in 'designfile', which the first column in datafile is name of ",
               "each row, and other columns are expression values corresponding to ",
               "sample name set in 'designfile'.\n") );
    }
    if (nrow(EXPD) == 1)
    {
        EXPD <- rbind(EXPD, EXPD);
        dimnames(EXPD)[[1]] <- c("1","2");
        EXP_IDS <- c("1","2");
        SIG_row <- TRUE;
    }
    EXP_IDNAME <<- EXPD[,1];
    names(EXP_IDNAME) <- EXP_IDS;
    # explibID <- dimnames(EXPD)[[2]];
    # explibID <- gsub("\\s+", "", explibID);
    # dimnames(EXPD)[[2]] <- explibID;
    SUBJECT_LSL <<- list();
    for (i in 1:length(UNI_SUBJECTID))
    {
        subject_designD <- SUBJECT_DGL[[i]];
        subject_libID <- subject_designD$libID;
        subject_timeT <- subject_designD$timeT;
        names(subject_timeT) <- 1:length(subject_timeT);
        order_timeT <- names(sort(subject_timeT));
        order_timeT <- as.numeric(order_timeT);
        expdataD <- EXPD[,subject_libID[order_timeT]];
        dimnames(expdataD)[[1]] <- EXP_IDS;
        dimnames(expdataD)[[2]] <- subject_libID[order_timeT];
        cat(paste("The LS is processing ", UNI_SUBJECTID[i],
                  ", the ", i, " in total ", length(UNI_SUBJECTID),
                  " subjects.\n", sep="") );
        lsoutM <- runCirCaseLS(indata=expdataD, LStime=subject_timeT[order_timeT], minper=minper, maxper=maxper);
        SUBJECT_LSL[[i]] <- lsoutM;
    }

    ####integrate and output results
    if (!length(SUBJECT_LSL))
    {
        cat(paste("Warning: No output result is from the analysis about '",
                  datafile, "', please carefully check the 'datafile' ",
                  "'designfile' and self-setted parameters in 'meta3d' ",
                  "function.\n", sep="") );
    }  else  {
        ##group the subjects according to groupID and then integrate the results
        ##with the function named as 'integrateLSOUT' in 'CirCaseSubF.R'.
        INTEG_LSL <<- list();
        subject_index <- 1:length(UNI_SUBJECTID);
        names(subject_index) <- UNI_SUBJECTID;
        if (length(GROUP_SUBJECTL))
        {
            for (k in 1:length(UNI_GROUPID))
            {
                subject_lsID <- as.character(GROUP_SUBJECTL[[k]]);
                subject_lsindex <- as.numeric(subject_index[subject_lsID]);
                i <- 1;
                subjects_integL <- list();
                for (m in subject_lsindex)
                {
                    subjects_integL[[i]] <- SUBJECT_LSL[[m]];
                    i <- i + 1;
                }
                INTEG_LSL[[k]] <- integrateLSOUT(lsoutL=subjects_integL, adjustV=ADPHAL, weightS=weightedMethod);
            }
        }  else  {
            INTEG_LSL[[1]] <- integrateLSOUT(lsoutL=SUBJECT_LSL, adjustV=ADPHAL, weightS=weightedMethod);
        }

        ##output separate analysis result for each subject-'SUBJECT_LSL'
        filename <- unlist(strsplit(designfile,.Platform$file.sep,fixed=TRUE));
        filename <- filename[length(filename)];
        if ( grepl("\\", filename, fixed=TRUE) )
        {
            filename <- unlist(strsplit(designfile,"\\",fixed=TRUE));
            filename <- filename[length(filename)];
        }  else if ( grepl("//", filename, fixed=TRUE) )  {
            filename <- unlist(strsplit(designfile,"//",fixed=TRUE));
            filename <- filename[length(filename)];
        }
        if ( INTEGRATION != "onlyIntegration" )
        {
            for ( i in 1:length(SUBJECT_LSL) )
            {
                subject_outname <- paste(OUTDIR, .Platform$file.sep, outSymbol, "LSresultSubjectID_", UNI_SUBJECTID[i], "_", filename, sep="");
                subject_outputM <- SUBJECT_LSL[[i]];
                subject_nameorder <- subject_outputM[,"CirID"];
                subject_outputM[,"CirID"] <- EXP_IDNAME[subject_nameorder];
                if (!SIG_row)
                {
                    write.table(subject_outputM, file=subject_outname, row.names=FALSE, sep=FILE_SEP, quote=FILE_QUOTE2, dec=FILE_DEC);
                }  else  {
                    sigsubject_outputM <- subject_outputM;
                    sigsubject_outputM[,"BH.Q"] <- sigsubject_outputM[,"p"];
                    sigsubject_outputM <- data.frame(sigsubject_outputM, stringsAsFactors=FALSE);
                    write.table(sigsubject_outputM[1,], file=subject_outname, row.names=FALSE, sep=FILE_SEP, quote=FILE_QUOTE2, dec=FILE_DEC);
                }
            }
        }

        ##output integrated results group by group-'INTEG_LSL'
        if ( INTEGRATION != "noIntegration" )
        {
            for ( i in 1:length(INTEG_LSL) )
            {
                if (length(GROUP_SUBJECTL))
                {
                    group_outname <- paste(OUTDIR, .Platform$file.sep, outSymbol, "CirCaseIntegrationGroupID_", UNI_GROUPID[i], "_", filename, sep="");
                }  else  {
                    group_outname <- paste(OUTDIR, .Platform$file.sep, outSymbol, "CirCaseIntegration_", filename, sep="");
                }
                group_outputM <- INTEG_LSL[[i]];
                group_nameorder <- group_outputM[,"CirID"];
                group_outputM[,"CirID"] <- EXP_IDNAME[group_nameorder];
                if (!SIG_row)
                {
                    write.table(group_outputM, file=group_outname, row.names=FALSE, sep=FILE_SEP, quote=FILE_QUOTE2, dec=FILE_DEC);
                }  else  {
                    group_outputM[, "CirCase_BH.Q"] <- group_outputM[, "CirCase_Pvalue"];
                    group_outputM <- data.frame(group_outputM, stringsAsFactors=FALSE);
                    write.table(group_outputM[1,], file=group_outname, row.names=FALSE, sep=FILE_SEP, quote=FILE_QUOTE2, dec=FILE_DEC);
                }
            }
        }
        cat(paste("DONE! The analysis about '", datafile, "' and '",
                   designfile, "'" , " has been finished.\n", sep=""));
    }

    ####output the analysis note
    print( c("Time used:",(proc.time()-run_start)) );
    cat("\n\n");
    obend <- ls();
    obend_index <- 1:length(obend);
    names(obend_index) <- obend;
    rm(list=obend[-obend_index[obst]]);
}
####End
