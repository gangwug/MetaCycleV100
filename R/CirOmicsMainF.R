##' Detect rhythmic signals from time-series datasets with multiple methods
##'
##' This is a function that incorporates ARSER, JTK_CYCLE and Lomb-Scargle to
##'   detect rhythmic signals from time-series datasets.
##'
##' \href{https://github.com/cauyrd/ARSER}{ARSER}(Yang, 2010),
##'   \href{http://openwetware.org/wiki/HughesLab:JTK_Cycle}{JTK_CYCLE}(
##'   Hughes, 2010), and
##'   \href{http://research.stowers-institute.org/efg/2005/LombScargle/}{
##'   Lomb-Scargle}(Glynn, 2006) are three popular methods of detecting
##'   rhythmic signals. \code{ARS} can not analyze unevenly sampled datasets,
##'   or enenly sampled datasets but with missing values, or with replicate
##'   samples, or with non-integer sampling interval. \code{JTK} is not
##'   suitable to analyze unevenly sampled datasets or evenly sampled datasets
##'   but with non-integer sampling interval. If set \code{analysisStrategy}
##'   as \code{"auto"}(default), \code{meta2d} will automatically select
##'   proper method from \code{cirMethod} for each input dataset. If the user
##'   clearly know that the dataset could be analyzed by each method defined
##'   by \code{cirMethod} and do not hope to output integrated values,
##'   \code{analysisStrategy} can be set as \code{"selfUSE"}.
##'
##'   \code{ARS} used here is translated from its python version which always
##'   uses \code{"yule-walker"}, \code{"burg"}, and \code{"mle"} methods(see
##'   \code{\link[stats]{ar}}) to fit autoregressive models to time-series
##'   data. Fitting by \code{"mle"} will be very slow for datasets
##'   with many time points. If \code{ARSmle = "auto"} is used,
##'   \code{meta2d} will only include \code{"mle"} when number of time points
##'   is smaller than 24. In addition, one evalution work(Wu, 2014) indicates
##'   that \code{ARS} shows relative high false positive rate in analyzing
##'   high-resolution datasets (1h/2days and 2h/2days). \code{JTK}(version 1)
##'   used here is minorly modified for easing its p-value calculation problem
##'   in analyzing datasets with missing values(reporting more significant
##'   p-values). \pkg{MetaCycle} will add a new version of \code{JTK}(may
##'   provide a better solution for this problem) when it is released.
##'
##'   The power of detecting rhythmic signals for an  algorithm is associated
##'   with the nature of data and interested periodic pattern(Deckard, 2013),
##'   which indicates that integrating analysis results from multiple methods
##'   may be helpful to rhythmic detection. For integrating p-values,
##'   Bonferroni correction(\code{"bonferroni"}) and Fisher's method(
##'   \code{"fisher"}) (Fisher, 1925; implementation code from \pkg{MADAM})
##'   could be selected, and \code{"bonferroni"} is usually more conservative
##'   than \code{"fisher"}. The integrated period is arithmetic mean of
##'   multiple periods. For integrating phase, \code{meta2d} takes use of
##'   \href{https://en.wikipedia.org/wiki/Mean_of_circular_quantities}{
##'   mean of circular quantities}. Integrated period and phase is further
##'   used to calculate the baseline value and amplitude through fitting a
##'   constructed periodic model.
##'
##'   Phases given by \code{JTK} and \code{LS} need to be adjusted with their
##'   predicted period (\code{adjustedPhase = "predictedPer"}) before
##'   integration. If \code{adjustedPhas = "notAdjusted"} is selected, no
##'   integrated phase will be calculated. If set \code{weightedPerPha} as
##'   \code{TRUE}, weighted scores will be used in averaging periods and
##'   phases. Weighted scores for one method are based on all its reported
##'   p-values, which means a weighted score assigned to any one profile will
##'   be affected by all other profiles. It is always a problem of averaging
##'   phases with quite different period lengths(eg. averaging two phases
##'   with 16-hours' and 30-hours' period length). Currently, setting
##'   \code{minper}, \code{maxper} and \code{ARSdefaultPer} to a same value
##'   may be the only way of completely eliminating such problem.
##'
##'   This function is originally aimed to analyze large scale perodic data(
##'   eg. circadian transcriptome data) without individual information.
##'   Please pay attention to data format of input file(see \code{Examples}
##'   part). Except the first column and first row, others are time-series
##'   experimental values(setting missing values as \code{NA}).
##'
##' @param infile a character string. The name of input file containing
##'   time-series data.
##' @param outdir a character string. The name of directory used to store
##'   output files.
##' @param filestyle a character vector(length 1 or 3). The data format of
##'   input file, must be \code{"txt"}, or \code{"csv"}, or a character vector
##'   containing field separator character(\code{sep}), quoting character
##'   (\code{quote}), and the character used for decimal points(\code{dec},
##'   for details see \code{\link[utils]{read.table}}).
##' @param timepoints a numeric vector corresponding to sampling time points
##'   of input time-series data; if sampling time points are in the first line
##'   of input file, it could be set as a character sting-"Line1".
##' @param minper a numeric value. The minimum period length of interested
##'   rhythms. The default is \code{20} for circadian rhythms.
##' @param maxper a numeric value. The maximum period length of interested
##'   rhythms. The default is \code{28} for circadian rhythms.
##' @param cirMethod a character vector(length 1 or 2 or 3). User-defined
##'   methods for detecting rhythmic signals, must be selected as any one, any
##'   two or all three methods(default) from \code{"ARS"}(ARSER),
##'   \code{"JTK"}(JTK_CYCLE) and \code{"LS"}(Lomb-Scargle).
##' @param analysisStrategy a character string. The strategy used to select
##'   proper methods from \code{cirMethod} for analyzing input time-series
##'   data, must be \code{"auto"}(default), or \code{"selfUSE"}. See
##'   \code{Details} part for more information.
##' @param outIntegrationFile a character string. This parameter controls what
##'   kinds of analysis results will be outputted, must be one of \code{"both"}
##'   (default), \code{"onlyIntegration"}(only output integration file), or
##'   \code{"noIntegration"}(not output integration file).
##' @param adjustedPhase a character string. The method used to adjust original
##'   phase calculated by each method in integration file, must be one of
##'   \code{"predictedPer"}(adjust phase with predicted period length) or
##'   \code{"notAdjusted"}(not adjust phase).
##' @param combinePvalue a character string. The method used to integrate
##'   multiple p-values, must be one of \code{"bonferroni"}(Bonferroni
##'   correction), or \code{"fisher"}(Fisher's method).
##' @param weightedPerPha logical. If \code{TRUE}, weighted scores based on
##'   p-value given by each method will be used to calculate the integrated
##'   period length and phase.
##' @param ARSmle a character string. The stragegy of using MLE method in
##'   \code{\link[stats]{ar}} fit of \code{"ARS"}, must be one of
##'   \code{"auto"}(use MLE depending the number of time points), \code{"mle"}
##'   (always use MLE), or \code{"nomle"}(never use MLE).
##' @param ARSdefaultPer a numeric value. The expected period length of
##'   interested rhythm, which is a necessary parameter for \code{ARS}. The
##'   default is \code{24}(for circadian rhythms). Set it to another proper
##'   numeric value for other rhythms.
##' @param outRawExpValue logical. If \code{TRUE}, time-series data will be
##'   added in the output files.
##' @param outSymbol a character string. A common prefix exists in the names of
##'   output files.
##' @return
##' \code{meta2d} will store analysis results in different files under
##'   \code{outdir} instead of returning them as objects. Files named with
##'   "ARSresult", "JTKresult" and "LSreult" store analysis results from
##'   \code{ARS}, \code{JTK} and \code{LS} respectively. The file named with
##'   "CirOmicsresult" is integration file, and it stores integrated values in
##'   columns with a common name-"CirOmics". The integration file also
##'   contains p-value, FDR value, period, phase(adjusted phase if
##'   \code{adjustedPhase = "predictedPer"}) and amplitude values calculated
##'   by each method.
##' @references
##' Yang R. and  Su Z. (2010). Analyzing circadian expression data by
##'   harmonic regression based on autoregressive spectral estimation.
##'   \emph{Bioinformatics}, \bold{26(12)}, i168--i174.
##'
##' Hughes M. E., Hogenesch J. B. and Kornacker K. (2010). JTK_CYCLE: an
##'   efficient nonparametric algorithm for detecting rhythmic components in
##'   genome-scale data sets. \emph{Journal of Biological Rhythms},
##'   \bold{25(5)}, 372--380.
##'
##' Glynn E. F., Chen J. and Mushegian A. R. (2006). Detecting periodic
##'   patterns in unevenly spaced gene expression time series using
##'   Lomb-Scargle periodograms. \emph{Bioinformatics}, \bold{22(3)},
##'   310--316.
##'
##' Wu G., Zhu J., Yu J., Zhou L., Huang J. Z. and  Zhang Z. (2014). Evaluation
##'   of five methods for genome-wide circadian gene identification.
##'   \emph{Journal of Biological Rhythms}, \bold{29(4)}, 231--242.
##'
##' Deckard A., Anafi R. C., Hogenesch J. B., Haase S.B. and Harer J. (2013).
##'   Design and analysis of large-scale biological rhythm studies:
##'   a comparison of algorithms for detecting periodic signals in biological
##'   data. \emph{Bioinformatics}, \bold{29(24)}, 3174--3180.
##'
##' Fisher, R.A. (1925). \emph{Statistical methods for research workers}.
##'   Oliver and Boyd (Edinburgh).
##'
##' Kugler K. G., Mueller L.A. and Graber A. (2010). MADAM - an open source
##'   toolbox for meta-analysis. \emph{Source Code for Biology and Medicine},
##'   \bold{5}, 3.
##' @examples
##' # write 'cycSimu4h2d', 'cycMouseLiverRNA' and 'cycYeastCycle' into three
##' # 'csv' files
##' write.csv(cycSimu4h2d, file="cycSimu4h2d.csv", row.names=FALSE)
##' write.csv(cycMouseLiverRNA, file="cycMouseLiverRNA.csv", row.names=FALSE)
##' write.csv(cycYeastCycle, file="cycYeastCycle.csv", row.names=FALSE)
##'
##' # write 'cycMouseLiverProtein' into a 'txt' file
##' write.table(cycMouseLiverProtein, file="cycMouseLiverProtein.txt",
##'   sep="\t", quote=FALSE, row.names=FALSE)
##'
##' # analyze 'cycMouseLiverRNA.csv' with JTK_CYCLE
##' meta2d(infile="cycMouseLiverRNA.csv", filestyle="csv", outdir="example",
##'   timepoints=18:65, cirMethod="JTK", outIntegrationFile="noIntegration")
##'
##' # analyze 'cycMouseLiverProtein.txt' with JTK_CYCLE and Lomb-Scargle
##' meta2d(infile="cycMouseLiverProtein.txt", filestyle="txt",
##'   outdir="example", timepoints=rep(seq(0, 45, by=3), each=3),
##'   cirMethod=c("JTK","LS"), outIntegrationFile="noIntegration")
##'
##' # analyze 'cycSimu4h2d.csv' with ARSER, JTK_CYCLE and Lomb-Scargle and
##' # output integration file with analysis results from each method
##' meta2d(infile="cycSimu4h2d.csv", filestyle="csv", outdir="example",
##'   timepoints="Line1")
##'
##' # analyze 'cycYeastCycle.csv' with ARSER, JTK_CYCLE and Lomb-Scargle to
##' # detect transcripts associated with cell cycle, and only output
##' # integration file
##' meta2d(infile="cycYeastCycle.csv",filestyle="csv", outdir="example",
##'   minper=80, maxper=96, timepoints=seq(2, 162, by=16),
##'   outIntegrationFile="onlyIntegration", ARSdefaultPer=85,
##'   outRawExpValue=TRUE)
##' @export

meta2d <- function(infile, outdir="cirout", filestyle, timepoints,
                   minper=20, maxper=28, cirMethod=c("ARS","JTK","LS"),
                   analysisStrategy="auto", outIntegrationFile="both",
                   adjustedPhase="predictedPer", combinePvalue="fisher",
                   weightedPerPha=FALSE, ARSmle="auto", ARSdefaultPer=24,
                   outRawExpValue=FALSE, outSymbol="")
{
    ##exist objects before running 'CirOmics.R'
    obst <- ls();
    run_start=proc.time();

    ####extract 'infile', 'outdir', 'minper', 'maxper',
    ####set by users and store them as global variables
    INFILE <<- infile;
    outdir2 <- unlist(strsplit(outdir,.Platform$file.sep));
    OUTDIR <<- paste(outdir2,collapse=.Platform$file.sep);
    ##create the directory if it is not exist
    if (! file.exists(OUTDIR) )
    { dir.create(OUTDIR); }
    ##the maximum and minimum period length should be positive number,
    ##and maximum period should not smaller than minimum period length
    if ( (!is.numeric(minper)) | (minper <= 0 ) )
    {	stop("The 'minper' should be a positive numeric value.\n");	}
    if ( (!is.numeric(maxper)) | (maxper <= 0 ) )
    {	stop("The 'maxper' should be a positive numeric value.\n");	}
    if (minper > maxper)
    {	stop("The 'minper' should not be larger than 'maxper'.\n");	}
    MINPER <<- minper;
    MAXPER <<- maxper;

    ####check whether the 'timepoints' is a effective numeric vector,
    ####or the 'timepoints' information is stored in the first line of 'infile'
    NUMT <<- TRUE;
    if (is.character(timepoints))
    {
        if ( (timepoints == "Line1") | (timepoints == "line1") )
        {   NUMT <<- FALSE; }
        else
        {
            stop(paste("If the time points values are in the first line of ",
                       "infile, please set a character value to 'timepoints' ",
                       "as 'Line1' or 'line1'.\n", sep="") );
        }
    }  else if ( is.numeric(timepoints) )  {
         if ( any(is.na(timepoints)) | any(is.nan(timepoints)) )
         {  stop("The 'timepoints' should not contain 'NA' or 'NaN'.\n");    }
         if (min(timepoints) < 0 )
         {
            stop("The 'timepoints' should be a non-negative numeric vector.\n");
         }
    }  else {
        stop(paste("The 'timepoints' should be set as a non-negative numeric ",
                   "vector value. If the first line of infile contains ",
                   "timepoints information, 'timepoints' could be set as a ",
                   "character value-'Line1' or 'line1'.\n", sep="") );
    }

    ####extract field separator character(FILE_SEP), set of quoting characters(FILE_QUOTE)
    ####and the character used for decimal points(FILE_DEC).
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
            stop(paste("Please set the input and output 'filestyle' before ",
                       "running 'meta2d'(the 'fileSyle' could be set as ",
                       "'txt' or 'csv'.\n", sep="") );
        }
    } else if (length(filestyle) > 1) {
        if (length(filestyle) == 3) {
            FILE_SEP <<- filestyle[1];
            FILE_QUOTE <<- filestyle[2];
            FILE_QUOTE2 <<- TRUE;
            FILE_DEC <<- filestyle[3];
        } else {
            stop(paste("Please set the input and output 'filestyle' before ",
                       "running meta2d(the 'filestyle' should be assigned a ",
                       "vector containing three characters, which corresponding",
                       " to symbols used to separating the columns,quoting the ",
                       "values and used for decimal points).\n", sep="") );
        }
    } else {
        stop(paste("Please set the input and output filestyle before running ",
                   "meta2d(the 'filestyle' should be set as 'txt', 'csv', or",
                   " assigned a vector containing three characters, which ",
                   "corresponding to symbols used to separating the columns, ",
                   "quoting the values and used for decimal points).\n",sep="") );
    }

    ####extract user-defined method for analyzing expression profiles
    CIRM <<- rep(FALSE,3);
    names(CIRM) <- c("ARS","JTK","LS");
    cirMethod <- unique(cirMethod);
    for (i in 1:length(cirMethod))
    {
        if ( (cirMethod[i] != "ARS") & (cirMethod[i] != "JTK") & (cirMethod[i] != "LS") )
        {
            stop(paste("Please check the parameter of 'cirMethod', one or multiple methods",
                       " could only be selected from 'ARS', 'JTK' and 'LS'.\n", sep="") );
        }
    }
    CIRM[cirMethod] <- TRUE;

    ####set default period for ARS
    ARS_PER <- "";
    if (length(ARSdefaultPer) == 1)
    {
        if (is.numeric(ARSdefaultPer))
        {
           ARS_PER <-  ARSdefaultPer;
           if (length(grep("ARS", cirMethod)) )
           {
                if ( (ARS_PER < minper) | (ARS_PER > maxper) )
                {
                    stop(paste("If hope to use ARS for this analysis, please set",
                               " 'ARSdefaultPer' no smaller than 'minper' and no ",
                               "larger than 'maxper'. By default, 'ARSdefaultPer' ",
                               "is set as expected period length in search.\n", sep="") );
                }
           }
        }  else  {
            stop(paste("If hope to use ARS for this analysis, please give one numeric",
                       "value to 'ARSdefaultPer'. By default, 'ARSdefaultPer' is set ",
                       "as expected period length in search.\n") );
        }
    } else {
        stop(paste("If hope to use ARS for this analysis, please give one numeric ",
                   "value to 'ARSdefaultPer'. By default, 'ARSdefaultPer' is set ",
                   "as expected period length in search.\n", sep="") );
    }

    ####read time-series values and extract each row's name
    if (NUMT)
    {
        EXP_dataframe <<- read.table(INFILE, header=TRUE, sep=FILE_SEP, quote=FILE_QUOTE, dec=FILE_DEC, stringsAsFactors=FALSE);
        ID_ORDER <<- dimnames(EXP_dataframe)[[1]];
        LIB_ID <<- dimnames(EXP_dataframe)[[2]];
    }  else  {
        EXP_dataframe <<- read.table(INFILE, header=FALSE, sep=FILE_SEP, quote=FILE_QUOTE, dec=FILE_DEC, stringsAsFactors=FALSE);
        dimnames(EXP_dataframe)[[2]] <- as.character(EXP_dataframe[1,]);
        timepoints <- as.numeric(EXP_dataframe[1, -1]);
        if ( any(is.na(timepoints)) | any(is.nan(timepoints)) )
        {
            stop(paste("Please check the time points information in the first line of",
                 " input file. Please remove non-numeric character of each time",
                 " points or replacing 'NA' or 'NaN' with numeric value.\n", sep="") );
        }
        if (min(timepoints) < 0 )
        {
            stop("The minimal time points value should not be a negative value.\n");
        }
        EXP_dataframe <- EXP_dataframe[-1,];
        dimnames(EXP_dataframe)[[1]] <- 1:nrow(EXP_dataframe);
        ID_ORDER <<- dimnames(EXP_dataframe)[[1]];
        LIB_ID <<- dimnames(EXP_dataframe)[[2]];
    }
    SIG_row <<- FALSE;
    if (ncol(EXP_dataframe) == 1)
    {
        stop(paste("The whole input file is stored in a one-column dataframe,",
                   "please check whether 'filestyle' is correctly set.\n", sep="") );
    }
    if (nrow(EXP_dataframe) == 0)
    {
        stop(paste("The whole input file contains only one line which is ",
                   "taken as column names by 'meta2d',",
                   " please check the input file.\n", sep="") );
    }
    if ( length(timepoints) != ( ncol(EXP_dataframe) - 1) )
    {
        stop(paste("Please check the input file or the 'timepoints',",
                   "the column number of input file should be one more ",
                   "than the length of time points, which the first ",
                   "column in the input file is the ID column",
                   " of each row, and other columns are expression values",
                   "corresponding to the setted time points.\n", sep="") );
    }
    if (nrow(EXP_dataframe) == 1)
    {
        EXP_dataframe <- rbind(EXP_dataframe, EXP_dataframe);
        dimnames(EXP_dataframe)[[1]] <- c("1","2");
        ID_ORDER <- c("1","2");
        SIG_row <- TRUE;
    }

    ####sort 'timepoints' value with increasing order
    names(timepoints) <- 1:length(timepoints);
    timepoints <- sort(timepoints);
    time_orderIndex <- as.numeric(names(timepoints));
    EXPM <<- as.matrix(EXP_dataframe[,-1]);
    dimnames(EXPM)[[1]] <- ID_ORDER;
    EXPM <- EXPM[,time_orderIndex];
    expm_ID <- LIB_ID[-1];
    expm_ID <- expm_ID[time_orderIndex];
    dimnames(EXPM)[[2]] <- expm_ID;
    LIB_ID <- c(LIB_ID[1], expm_ID);
    EXP_dataframe[ID_ORDER, 2:ncol(EXP_dataframe)] <- EXPM[ID_ORDER,];
    dimnames(EXP_dataframe)[[2]] <- LIB_ID;
    OUT_ID <<- as.character(EXP_dataframe[,1]);
    names(OUT_ID) <- ID_ORDER;

    ####When total number of time points is >= 24, 'mle' will not be
    ####used for estimating AR coefficients('mle' method is too slow).
    uni_timepoints <- unique(timepoints);
    if ( length(uni_timepoints) < 4 )
    {
        cat(paste("Warning: the number of time points is too small, ",
                  "it may be not a good choice to analyze this expression file ",
                  "with 'ARS','JTK' or 'LS'. If possible, please carefully ",
                  "check the output results and try other useful methods.\n", sep="") );
    }
    ARS_MET <- "";
    if ( (ARSmle == "auto") & ( length(uni_timepoints) < 24) ) {
        ARS_MET <- c("yule-walker","mle","burg");
    } else if ( (ARSmle == "auto") & (length(uni_timepoints) >= 24) ) {
        ARS_MET <- c("yule-walker","burg");
    } else if (ARSmle == "mle") {
        ARS_MET <- c("yule-walker","mle","burg");
    } else if (ARSmle == "nomle") {
        ARS_MET <- c("yule-walker","burg");
    } else {
        stop(paste("Please check the parameter of 'ARSmle', it ",
                   "should be set as 'auto', 'mle' or 'nomle'.\n", sep="") );
    }

    ####check whether output integration results at the end of analysis
    INTEGRATION <- "";
    if (outIntegrationFile == "both") {
        INTEGRATION <- "both";
    } else if (outIntegrationFile == "onlyIntegration") {
        INTEGRATION <- "onlyIntegration";
    } else if (outIntegrationFile == "noIntegration") {
        INTEGRATION <- "noIntegration";
    } else {
        stop(paste("Please check the parameter of 'outIntegrationFile', ",
                   "it should be set as 'both', 'onlyIntegration' ",
                   "or 'noIntegration'.\n", sep="") );
    }

    ####check defined period length used to adjust phase
    ADPHA <<- "";
    if (adjustedPhase == "predictedPer") {
        ADPHA <- "PERIOD";
    } else if (adjustedPhase == "notAdjusted") {
        ADPHA <- 0;
    } else {
        stop("The 'adjustedPhase' should be set as 'predictedPer', or 'notAdjusted'.\n");
    }

    ####check which method is used to integrate multiple P-values
    if ( (combinePvalue != "Bonferroni") & (combinePvalue != "bonferroni") & (combinePvalue != "Fisher") & (combinePvalue != "fisher") )
    {
        stop(paste("The 'combinePvalue' should be set as ",
                   "'bonferroni' (equal to 'Bonferroni'), ",
                   "or 'Fisher' (equal to 'fisher').\n", sep="") );
    }
    if ( !is.logical(weightedPerPha) )
    {
        stop("The 'weightedPerPha' should be set as a logical value (TRUE or FALSE).\n");
    }
    if ( !is.logical(outRawExpValue) )
    {
        stop("The 'outRawExpValue' should be set as a logical value (TRUE or FALSE).\n");
    }

    ####extract key features of input dataset, including with/without non-integer interval,
    ####even/uneven sampling, with/without missing values, with/without replicates
    START_TIME <<- uni_timepoints[1];
    END_TIME <<- uni_timepoints[length(uni_timepoints)];
    MISSING_VALUE <<- FALSE;
    WITH_REPLICATE <<- FALSE;
    non_integerInterval <- FALSE;
    uneven_interval <- FALSE;
    if ( !all( round(diff(uni_timepoints)) == diff(uni_timepoints) ) )
    {	non_integerInterval <- TRUE;	}
    if ( length( unique(diff(uni_timepoints)) ) > 1 )
    {	uneven_interval <- TRUE;	}
    if ( (!all( !is.na(EXPM) )) | (!all( !is.nan(EXPM) )) )
    {	MISSING_VALUE <<- TRUE;	}
    if ( length(timepoints) != length(uni_timepoints) )
    {	WITH_REPLICATE <<- TRUE;	}

    ####set the output file name
    ARS_OUTM <<- "None_output";
    JTK_OUTM <<- "None_output";
    LS_OUTM <<- "None_output";
    filename <- unlist(strsplit(INFILE,.Platform$file.sep,fixed=TRUE));
    filename <- filename[length(filename)];
    if ( grepl("\\", filename, fixed=TRUE) )
    {
        filename <- unlist(strsplit(INFILE,"\\",fixed=TRUE));
        filename <- filename[length(filename)];
    }  else if ( grepl("//", filename, fixed=TRUE) )  {
        filename <- unlist(strsplit(INFILE,"//",fixed=TRUE));
        filename <- filename[length(filename)];
    }
    if (length(outSymbol) > 1)
    {
        outSymbol <- as.character(outSymbol);
        outSymbol <- paste(outSymbol, collapse="");
    }
    ars_outname <- paste(OUTDIR, .Platform$file.sep, outSymbol, "ARSresult_", filename,sep="");
    jtk_outname <- paste(OUTDIR, .Platform$file.sep, outSymbol, "JTKresult_", filename,sep="");
    ls_outname <- paste(OUTDIR, .Platform$file.sep,  outSymbol, "LSresult_", filename,sep="");
    outfile_name <- c(ars_outname, jtk_outname, ls_outname);
    names(outfile_name) <- c("ARS","JTK","LS");
    outfile_tag <<- rep(0,3);
    names(outfile_tag) <<- c("ARS","JTK","LS");
    CirOmics_outname <- paste(OUTDIR, .Platform$file.sep, outSymbol, "CirOmicsresult_", filename,sep="");

    ####select proper method to analyze profiles depending on sampling pattern
    ####ARS (even sampling & without non-integer intervals & without missing values & without replicates)
    ####JTK (even sampling & without non-integer intervals)
    ####LS is not restricted in analyzing profiles in current design
    if (analysisStrategy == "auto") {
        if (non_integerInterval) {
            if (CIRM["LS"])
            {
                LS_OUTM <<- runLS(EXP_dataframe,timepoints,minper=MINPER,maxper=MAXPER);
                outfile_tag["LS"] <<- 1;
            }
        } else if (uneven_interval) {
            if (CIRM["LS"])
            {
                LS_OUTM <<- runLS(EXP_dataframe,timepoints,minper=MINPER,maxper=MAXPER);
                outfile_tag["LS"] <<- 1;
            }
        } else if (MISSING_VALUE) {
            if (CIRM["JTK"])
            {
                JTK_OUTM<<- runJTK(EXP_dataframe,timepoints,minper=MINPER,maxper=MAXPER);
                outfile_tag["JTK"] <<- 1;
            }
            if (CIRM["LS"])
            {
                LS_OUTM<<- runLS(EXP_dataframe,timepoints,minper=MINPER,maxper=MAXPER);
                outfile_tag["LS"] <<- 1;
            }
        } else if (WITH_REPLICATE) {
            if (CIRM["JTK"])
            {
                JTK_OUTM <<- runJTK(EXP_dataframe,timepoints,minper=MINPER,maxper=MAXPER);
                outfile_tag["JTK"] <<- 1;
            }
            if (CIRM["LS"])
            {
                LS_OUTM <<- runLS(EXP_dataframe,timepoints,minper=MINPER,maxper=MAXPER);
                outfile_tag["LS"] <<- 1;
            }
        } else if ( (!non_integerInterval)&(!uneven_interval)&(!MISSING_VALUE)&(!WITH_REPLICATE) ) {
            if (CIRM["ARS"])
            {
                ARS_OUTM <<- runARS(EXP_dataframe, timepoints, minper=MINPER, maxper=MAXPER, arsper=ARS_PER, arsmet=ARS_MET);
                outfile_tag["ARS"] <<- 1;
            }
            if (CIRM["JTK"])
            {
                JTK_OUTM <<- runJTK(EXP_dataframe,timepoints,minper=MINPER,maxper=MAXPER);
                outfile_tag["JTK"] <<- 1;
            }
            if (CIRM["LS"])
            {
                LS_OUTM <<- runLS(EXP_dataframe,timepoints,minper=MINPER,maxper=MAXPER);
                outfile_tag["LS"] <<- 1;
            }
        } else {
            cat("Sorry for this bug, please contact the author.\n");
        }
    } else if (analysisStrategy == "selfUSE") {
        cat(paste("Warning: the parameter of 'analysisStrategy' is set as 'selfUSE', ",
                  "please be sure that the 'infile' could be analyzed by all methods ",
                  "that defined by 'cirMethod'. If not sure about that, please ",
                  "set the parameter of 'analysisStrategy' as 'auto'.\n", sep="") );
        if (CIRM["ARS"])
        {
            ARS_OUTM <<- runARS(EXP_dataframe, timepoints, minper=MINPER, maxper=MAXPER, arsper=ARS_PER, arsmet=ARS_MET);
            outfile_tag["ARS"] <<- 1;
        }
        if (CIRM["JTK"])
        {
            JTK_OUTM <<- runJTK(EXP_dataframe,timepoints,minper=MINPER,maxper=MAXPER);
            outfile_tag["JTK"] <<- 1;
        }
        if (CIRM["LS"])
        {
            LS_OUTM <<- runLS(EXP_dataframe,timepoints,minper=MINPER,maxper=MAXPER);
            outfile_tag["LS"] <<- 1;
        }
    } else {
        stop(paste("Please check the parameter of 'analysisStrategy', ",
                   "it should be set as 'auto' or 'selfUSE'.\n", sep="") );
    }

    ####output separate result file from each used method
    integration_header <- "CirID";
    if (INTEGRATION != "onlyIntegration")
    {
        ##output ARS analysis results
        if (outfile_tag["ARS"])
        {
            if (length(ARS_OUTM) > 1) {
                if (!SIG_row)
                {
                    if (!outRawExpValue)
                    {
                        write.table(ARS_OUTM,file=outfile_name["ARS"],row.names=FALSE,sep=FILE_SEP,quote=FILE_QUOTE2,dec=FILE_DEC);
                    }  else  {
                        write.table(cbind(ARS_OUTM[ID_ORDER,],EXPM[ID_ORDER,]),file=outfile_name["ARS"],row.names=FALSE,sep=FILE_SEP,quote=FILE_QUOTE2,dec=FILE_DEC);
                    }
                }  else  {
                    if (!outRawExpValue)
                    {
                        sigoutarsM <- data.frame(ARS_OUTM, stringsAsFactors=FALSE);
                    }  else  {
                        sigoutarsM <- data.frame(cbind(ARS_OUTM[ID_ORDER,], EXPM[ID_ORDER,]), stringsAsFactors=FALSE);
                    }
                    sigoutarsM[,"fdr_BH"] <- sigoutarsM[,"pvalue"];
                    write.table(sigoutarsM[1,],file=outfile_name["ARS"],row.names=FALSE,sep=FILE_SEP,quote=FILE_QUOTE2,dec=FILE_DEC);
                }
                integration_header <- c(integration_header,"ARS_pvalue","ARS_BH.Q","ARS_period","ARS_adjphase","ARS_amplitude");
            } else {
                cat(paste("Warning: no output result from ARS method. ",
                          "If 'analysisStrategy' was set as 'selfUSE', ",
                          "please check whether ARS is suitable for analyzing ",
                          "the input file. If 'analysisStrategy' was set as 'auto', ",
                          "please contact the author of this package.\n", sep="") );
            }
        }
        ##output JTK analysis results
        if (outfile_tag["JTK"])
        {
            if (length(JTK_OUTM) > 1) {
                if (!SIG_row)
                {
                    if (!outRawExpValue)
                    {
                        write.table(JTK_OUTM,file=outfile_name["JTK"],row.names=FALSE,sep=FILE_SEP,quote=FILE_QUOTE2,dec=FILE_DEC);
                    }  else  {
                        write.table(cbind(JTK_OUTM[ID_ORDER,],EXPM[ID_ORDER,]),file=outfile_name["JTK"],row.names=FALSE,sep=FILE_SEP,quote=FILE_QUOTE2,dec=FILE_DEC);
                    }
                }  else  {
                    if (!outRawExpValue)
                    {
                        sigoutjtkM <- data.frame(JTK_OUTM, stringsAsFactors=FALSE);
                    }  else  {
                        sigoutjtkM <- data.frame(cbind(JTK_OUTM[ID_ORDER,], EXPM[ID_ORDER,]) , stringsAsFactors=FALSE);
                    }
                    sigoutjtkM[,"BH.Q"] <- sigoutjtkM[,"ADJ.P"];
                    write.table(sigoutjtkM[1,], file=outfile_name["JTK"], row.names=FALSE, sep=FILE_SEP, quote=FILE_QUOTE2, dec=FILE_DEC);
                }
                integration_header <- c(integration_header,"JTK_pvalue","JTK_BH.Q","JTK_period","JTK_adjphase","JTK_amplitude");
            } else {
                cat(paste("Warning: no output result from JTK method. ",
                          "If 'analysisStrategy' was set as 'selfUSE', ",
                          "please check whether JTK is suitable for analyzing ",
                          "the input file. If 'analysisStrategy' was set as 'auto', ",
                          "please contact the author of this package.\n", sep="") );
            }
        }
        ##output LS analysis results
        if (outfile_tag["LS"])
        {
            if (length(LS_OUTM) > 1) {
                if (!SIG_row)
                {
                    if (!outRawExpValue)
                    {
                        write.table(LS_OUTM, file=outfile_name["LS"],row.names=FALSE,sep=FILE_SEP,quote=FILE_QUOTE2,dec=FILE_DEC);
                    }  else  {
                        write.table(cbind(LS_OUTM[ID_ORDER,], EXPM[ID_ORDER,]), file=outfile_name["LS"],row.names=FALSE,sep=FILE_SEP,quote=FILE_QUOTE2,dec=FILE_DEC);
                    }
                }  else  {
                    if (!outRawExpValue)
                    {
                        sigoutlsM <- data.frame(LS_OUTM, stringsAsFactors=FALSE);
                    }  else  {
                        sigoutlsM <- data.frame(cbind(LS_OUTM[ID_ORDER,], EXPM[ID_ORDER,]), stringsAsFactors=FALSE);
                    }
                    write.table(sigoutlsM[1,],file=outfile_name["LS"],row.names=FALSE,sep=FILE_SEP,quote=FILE_QUOTE2,dec=FILE_DEC);
                }
                integration_header <- c(integration_header,"LS_pvalue","LS_BH.Q","LS_period","LS_adjphase","LS_amplitude");
            } else {
                cat(paste("Warning: no output result from LS method. ",
                          "If 'analysisStrategy' was set as 'selfUSE', ",
                          "please check whether LS is suitable for analyzing ",
                          "the input file. If 'analysisStrategy' was set as 'auto', ",
                          "please contact the author of this package.\n", sep="") );
            }
        }
    }
    ####integration step and output integrated results
    AMPTIM <<- timepoints;
    if (  INTEGRATION != "noIntegration" )
    {
        ##integration step
        integration_outM <- "None_output";
        integration_num <- length(outfile_tag[outfile_tag > 0]);
        if (integration_num > 1) {
            ##'infile' is analyzed by two or three methods
            outLIST <- list("ARS"=ARS_OUTM,"JTK"=JTK_OUTM,"LS"=LS_OUTM);
            integration_outM <- ID_ORDER;
            pvaM <- ID_ORDER;
            phaM <- ID_ORDER;
            perM <- ID_ORDER;
            out_index <- as.numeric(which(outfile_tag == 1));
            out_index <- sort(out_index);
            for (index in out_index)
            {
                sep_outM <- outLIST[[index]];
                if (length(sep_outM) > 1) {
                    if (index == 1) {
                        ars_adjM <- adjPhaARS(sep_outM[ID_ORDER,],adjustV=ADPHA);
                        integration_outM <- cbind(integration_outM,sep_outM[ID_ORDER,c("pvalue","fdr_BH")],ars_adjM);
                        pvaM <- cbind(pvaM,sep_outM[ID_ORDER,"pvalue"]);
                        phaM <- cbind(phaM,ars_adjM[ID_ORDER,"phase"]);
                        perM <- cbind(perM,ars_adjM[ID_ORDER,"period"]);
                        if (INTEGRATION == "onlyIntegration")
                        {
                            integration_header <- c(integration_header,"ARS_pvalue","ARS_BH.Q","ARS_period","ARS_adjphase","ARS_amplitude");
                        }
                    } else if (index == 2) {
                        jtk_adjpha <- adjPhaJTK(sep_outM[ID_ORDER,],adjustV=ADPHA);
                        integration_outM <- cbind(integration_outM, sep_outM[ID_ORDER,c("ADJ.P","BH.Q","PER")], jtk_adjpha[ID_ORDER], sep_outM[ID_ORDER,"AMP"]);
                        pvaM <- cbind(pvaM, sep_outM[ID_ORDER,"ADJ.P"]);
                        phaM <- cbind(phaM, jtk_adjpha[ID_ORDER]);
                        perM <- cbind(perM, sep_outM[ID_ORDER,"PER"]);
                        if (INTEGRATION == "onlyIntegration")
                        {
                            integration_header <- c(integration_header,"JTK_pvalue","JTK_BH.Q","JTK_period","JTK_adjphase","JTK_amplitude");
                        }
                    } else if (index == 3) {
                        ls_adjpha <- adjPhaLS(sep_outM[ID_ORDER,],adjustV=ADPHA);
                        integration_outM <- cbind(integration_outM, sep_outM[ID_ORDER,c("p","BH.Q","Period")], ls_adjpha[ID_ORDER], sep_outM[ID_ORDER,"PhaseShiftHeight"]);
                        pvaM <- cbind(pvaM, sep_outM[ID_ORDER,"p"]);
                        phaM <- cbind(phaM, ls_adjpha[ID_ORDER]);
                        perM <- cbind(perM, sep_outM[ID_ORDER,"Period"]);
                        if (INTEGRATION == "onlyIntegration")
                        {
                            integration_header <- c(integration_header,"LS_pvalue","LS_BH.Q","LS_period","LS_adjphase","LS_amplitude");
                        }
                    }
                } else {
                    cat(paste("Warning: During integration process, no result from ", names(outfile_tag[index]),
                               " method. If 'analysisStrategy' was set as 'selfUSE', ", "please check whether ",
                               names(outfile_tag[index]), " is suitable for analyzing the input file. ",
                               "If 'analysisStrategy' was set as 'auto', ",
                               "please contact the author of package.\n", sep="") );
                }
            }
            ##pvalue adjust, period average, phase average, amplitude calculation
            if (ncol(pvaM) > 1)
            {
                if (nrow(pvaM) >= 2) {
                    cirpvaM <- getCirOmicsPVA(pvaM, method=combinePvalue);
                    if (adjustedPhase == "notAdjusted")
                    {
                        integration_outM <- cbind(integration_outM, cirpvaM[ID_ORDER,]);
                    }  else  {
                        cirperphaM <- getCirOmicsPerPha(periodM=perM, phaseM=phaM, pvalueM=pvaM, adjustV=ADPHA, weightedPerPha);
                        EXPMPERPHA <- cbind(EXPM[ID_ORDER,],cirperphaM[ID_ORDER,]);
                        ampM <- getCirOmicsAMP(exprperphaM=EXPMPERPHA);
                        rownames(ampM) <- ID_ORDER;
                        integration_outM <- cbind(integration_outM, cirpvaM[ID_ORDER,], cirperphaM[ID_ORDER,], ampM[ID_ORDER,]);
                    }
                } else {
                    cat(paste("Warning: there is a bug associated with ",
                              "'INTEGRATION' step, please contact the ",
                              "author of this package.\n", sep="") );
                }
            }
            if (adjustedPhase == "notAdjusted")
            {
                integration_header <- c(integration_header, "CirOmics_pvalue","CirOmics_BH.Q");
            }  else  {
                integration_header <- c(integration_header, "CirOmics_pvalue","CirOmics_BH.Q","CirOmics_period",
                                        "CirOmics_phase", "CirOmics_baselineEXP", "CirOmics_AMP", "CirOmics_relativeAMP");
            }
        } else {
            ##'infile' is analyzed by only one method, ARS or JTK or LS
            if (length(ARS_OUTM) > 1)
            {
                ars_adjM <- adjPhaARS(ARS_OUTM[ID_ORDER,],adjustV=ADPHA);
                arsEXPMPERPHA <- cbind(EXPM[ID_ORDER,],ars_adjM[ID_ORDER,1:2]);
                arsampM <- getCirOmicsAMP(exprperphaM=arsEXPMPERPHA);
                rownames(arsampM) <- ID_ORDER;
                integration_outM <- cbind(integration_outM,ARS_OUTM[ID_ORDER,c("pvalue","fdr_BH")],ars_adjM[ID_ORDER,], arsampM[ID_ORDER,]);
                if (INTEGRATION == "onlyIntegration")
                {
                    integration_header <- c(integration_header,"ARS_pvalue","ARS_BH.Q","ARS_period","ARS_adjphase","ARS_amplitude");
                }
                integration_header <- c(integration_header, "CirOmics_baselineEXP", "CirOmics_AMP", "CirOmics_relativeAMP");
            }
            if (length(JTK_OUTM) > 1)
            {
                jtk_adjpha <- adjPhaJTK(JTK_OUTM[ID_ORDER,],adjustV=ADPHA);
                jtkEXPMPERPHA <- cbind(EXPM[ID_ORDER,], JTK_OUTM[ID_ORDER,"PER"], jtk_adjpha[ID_ORDER]);
                jtkampM <- getCirOmicsAMP(exprperphaM=jtkEXPMPERPHA);
                rownames(jtkampM) <- ID_ORDER;
                integration_outM <- cbind(integration_outM,JTK_OUTM[ID_ORDER,c("ADJ.P","BH.Q","PER")],jtk_adjpha[ID_ORDER],JTK_OUTM[ID_ORDER,"AMP"], jtkampM[ID_ORDER,]);
                if (INTEGRATION == "onlyIntegration")
                {
                    integration_header <- c(integration_header,"JTK_pvalue","JTK_BH.Q","JTK_period","JTK_adjphase","JTK_amplitude");
                }
                integration_header <- c(integration_header, "CirOmics_baselineEXP", "CirOmics_AMP", "CirOmics_relativeAMP");
            }
            if (length(LS_OUTM) > 1)
            {
                ls_adjpha <- adjPhaLS(LS_OUTM[ID_ORDER,],adjustV=ADPHA);
                lsEXPMPERPHA <- cbind(EXPM[ID_ORDER,], LS_OUTM[ID_ORDER, "Period"], ls_adjpha[ID_ORDER]);
                lsampM <- getCirOmicsAMP(exprperphaM=lsEXPMPERPHA);
                rownames(lsampM) <- ID_ORDER;
                integration_outM <- cbind(integration_outM,LS_OUTM[ID_ORDER,c("p","BH.Q","Period")],ls_adjpha[ID_ORDER],LS_OUTM[ID_ORDER,"PhaseShiftHeight"], lsampM[ID_ORDER,]);
                if (INTEGRATION == "onlyIntegration")
                {
                    integration_header <- c(integration_header,"LS_pvalue","LS_BH.Q","LS_period","LS_adjphase","LS_amplitude");
                }
                integration_header <- c(integration_header, "CirOmics_baselineEXP", "CirOmics_AMP", "CirOmics_relativeAMP");
            }
        }
        ##output integrated results
        if (!is.vector(integration_outM)) {
            if ( ncol(integration_outM) > 1 )
            {
                integration_outM <- cbind(OUT_ID[ID_ORDER],integration_outM[,-1]);
                rownames(integration_outM) <- ID_ORDER;
                colnames(integration_outM) <- integration_header;
                if (!SIG_row)
                {
                    if (!outRawExpValue)
                    {
                        write.table(integration_outM,file=CirOmics_outname,row.names=FALSE,sep=FILE_SEP,quote=FILE_QUOTE2,dec=FILE_DEC);
                    }  else  {
                        write.table(cbind(integration_outM[ID_ORDER,], EXPM[ID_ORDER,]),file=CirOmics_outname,row.names=FALSE,sep=FILE_SEP,quote=FILE_QUOTE2,dec=FILE_DEC);
                    }
                }  else  {
                    if (!outRawExpValue)
                    {
                        integration_outM <- data.frame(integration_outM,stringsAsFactors=FALSE);
                    }  else  {
                        integration_outM <- data.frame(cbind(integration_outM[ID_ORDER,], EXPM[ID_ORDER,]),stringsAsFactors=FALSE);
                    }
                    pva_index <- grep("_pvalue", integration_header);
                    qva_index <- grep("_BH.Q", integration_header);
                    integration_outM[,qva_index] <- integration_outM[,pva_index];
                    integration_outM <- integration_outM[1,];
                    write.table(integration_outM,file=CirOmics_outname,row.names=FALSE,sep=FILE_SEP,quote=FILE_QUOTE2,dec=FILE_DEC);
                }
            } else {
                cirMethod <- paste(cirMethod, collapse=",");
                cat(paste("Warning: no integration results from ", cirMethod, " methods. ",
                           "It seems that there is a bug in the integration step, ",
                           "please contact the author of this package.\n", sep="") );
            }
        }
    }
    ####output analysis note on the screen
    if (sum(outfile_tag) == 0)
    {
        cirMethod <- paste(cirMethod, collapse=",");
        cat(paste("Warning: no output result from", cirMethod, "methods. ",
                   "Please check whether at least one of ",
                   cirMethod, " is suitable for analyzing the input file, ",
                   "and whether all parameters are correctly set as suggest. ",
                   "If it is quite sure that at least one of " ,
                   cirMethod, "can analyze the input file and ",
                   "parameters are correctly set as suggest, ",
                   "please contact the author of 'meta2d'.\n",sep="") );
    } else {
        cat(paste("DONE! The analysis about '",INFILE,"'"," has been finished.\n",sep=""));
    }
    print( c("Time used:",(proc.time()-run_start)) );
    cat("\n\n");
    obend <- ls();
    obend_index <- 1:length(obend);
    names(obend_index) <- obend;
    rm(list=obend[-obend_index[obst]]);
}
####End
