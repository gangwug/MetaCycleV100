% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/CirCaseMainF.R
\name{meta3d}
\alias{meta3d}
\title{Detect rhythimic signals from time-series datasets with individual
  information}
\usage{
meta3d(datafile, designfile, outdir = "cirout", filestyle, design_libColm,
  design_subjectColm, minper = 20, maxper = 28, timeUnit = "hour",
  design_hrColm, design_dayColm = NULL, design_minColm = NULL,
  design_secColm = NULL, design_groupColm = NULL,
  design_libIDrename = NULL, adjustedPhase = "predictedPer",
  weightedMethod = TRUE, combinePvalue = "fisher",
  outIntegrationFile = "both", dayZeroBased = FALSE, outSymbol = "")
}
\arguments{
\item{datafile}{a character string. The name of data file containing
time-sereis experimental values of all individuals.}

\item{designfile}{a character string. The name of experimental design file,
at least containing the library ID(column names of \code{datafile}),
subject ID(the individual corresponding to each library ID), and
sampling time information of each library ID.}

\item{outdir}{a character string. The name of directory used to store
output files.}

\item{filestyle}{a character vector(length 1 or 3). The data format of
input files, must be \code{"txt"}, or \code{"csv"}, or a character
vector containing field separator character(\code{sep}), quoting
character(\code{quote}), and the character used for decimal
points(\code{dec}, for details see \code{\link[utils]{read.table}}).}

\item{design_libColm}{a numeric value. The order index(from left to right)
of the column storing library ID in \code{designfile}.}

\item{design_subjectColm}{a numeric value. The order index(from left to
right) of the column storing subject ID in \code{designfile}.}

\item{minper}{a numeric value. The minimum period length of interested
rhythms. The default is \code{20} for circadian rhythms.}

\item{maxper}{a numeric value. The maximum period length of interested
rhythms. The default is \code{28} for circadian rhythms.}

\item{timeUnit}{a character string. The basic time-unit, must be one of
\code{"day"}, \code{"hour"}(default for circadian study),
\code{"minute"}, or \code{"second"} depending on specific experimental
design.}

\item{design_hrColm}{a numeric value. The order index(from left to right)
of the column storing time point value-sampling hour information in
\code{designfile}. If there is no such column in \code{designfile},
set it as \code{NULL}.}

\item{design_dayColm}{a numeric value. The order index(from left to right)
of the column storing time point value-sampling day information in
\code{designfile}. If there is no such column in \code{designfile},
set it as \code{NULL}(default).}

\item{design_minColm}{a numeric value. The order index(from left to right)
of the column storing time point value-sampling minute information in
\code{designfile}. If there is no such column in \code{designfile},
set it as \code{NULL}(default).}

\item{design_secColm}{a numeric value. The order index(from left to right)
of the column storing time point value-sampling second information in
\code{designfile}. If there is no such column in \code{designfile},
set it as \code{NULL}(default).}

\item{design_groupColm}{a numeric value. The order index(from left to
right) of the column storing experimental group information of each
individual in \code{designfile}. If there is no such column in
\code{designfile}, set it as \code{NULL}(default) and take all
individuals as one group.}

\item{design_libIDrename}{a character vector(length 2) containing  a
matchable character string in each library ID of \code{designfile}, and
a replacement character string. If it is not necessary to replace
characters in library ID of \code{designfile}, set it as \code{NULL}(
default).}

\item{adjustedPhase}{a character string. The method used to adjust each
phase calculated by Lomb-Scargle before getting integrated phase,
must be one of \code{"predictedPer"}(adjust phase with predicted
period length) or \code{"notAdjusted"}(not adjust phase).}

\item{weightedMethod}{logical. If \code{TRUE}(default), weighted score
based on p-value of each individual will be used to integrate period,
phase and amplitude values of multiple individuals.}

\item{combinePvalue}{a character string. The method used to integrate
p-values of multiple individuals, currently only \code{"fisher"}(Fisher's
method) could be selected.}

\item{outIntegrationFile}{a character string. This parameter controls what
kinds of analysis results will be outputted, must be one of \code{"both"}
(default), \code{"onlyIntegration"}(only output integrated analysis
results of each experimental group), or \code{"noIntegration"}(only
output analysis results of each individual).}

\item{dayZeroBased}{logical. If \code{TRUE}, the first sampling day is
recorded as day zero in the \code{designfile}.}

\item{outSymbol}{a character string. A common prefix exists in the names of
output files.}
}
\value{
\code{meta3d} will write analysis results to \code{outdir} instead of
  returning them as objects. Output files with "LSresultSubjectID" in
  the file name are Lomb-Scargle analysis results for each individual.
  Files named with "CirCaseIntegrationGroupID" store integrated p-values,
  period, phase, baseline, amplitude and relative amplitude values
  from multiple individuals of each group and calculated FDR values based
  on integrated p-values.
}
\description{
This is a function that takes use of Lomb-Scargle to detect rhythmic
  signals from time-series datasets containing individual information.
}
\details{
This function is originally aimed to analyze large scale perodic data with
  individual information. Please pay attention to the data format of
  \code{datafile} and \code{designfile}(see \code{Examples} part).
  Time-series experimental values(missing values as \code{NA}) from
  all individuals shoud be stored in \code{datafile}, with the first row
  containing all library ID(unique identification number for each sample)
  and the first column containing all detected molecular names(eg.
  transcript ID or gene name). The \code{designfile} should at least have
  three columns-library ID, subject ID and sampling time collumn.
  Experimental group information of each subject ID may be in another
  column. In addition, sampling time information may be stored in multiple
  columns instead of one column. For example, sampling time-"36 hours" may
  be recorded as "day 2"(sampling day column, \code{design_dayColm}) plus
  "12 hours"(sampling hour column, \code{design_hrColm}). The library ID
  in \code{datafile} and \code{designfile} should be same. If there are
  different characters between library ID in these two files, try
  \code{design_libIDrename} to keep them same.

\href{http://research.stowers-institute.org/efg/2005/LombScargle/}{
  Lomb-Scargle}(Glynn, 2006) is used to analyze time-series profiles
  individual by individual. The original code is modified and added with
  calculating baseline value and amplitude(using similar method
  as \code{\link{meta2d}}). Then p-values, period, phase, baseline
  value, amplitude and relative amplitude from multiple individuals are
  integrated group by group. P-values from different individuals are
  integrated with Fisher's method(\code{"fisher"}) (Fisher,1925;
  implementation code from \pkg{MADAM}). The integrated period, baseline
  value, amplitude and relative amplitude are arithmetic mean of multiple
  individuals, respectively. The integrated phase is
  \href{https://en.wikipedia.org/wiki/Mean_of_circular_quantities}{
  mean of circular quantities}(\code{adjustedPhase = "predictedPer"}) or
  a arithmetic mean(\code{adjustedPhase = "notAdjusted"}) of multiple
  individual phases. If \code{weightedMethod = TRUE} is selected, weighted
  scores(\code{-log10(p-values)}) will be taken into account in calculation
  of integrated period, phase, baseline, amplitude and relative amplitude.

For completly removing the potential problem of averaging phases with
  quite different period length(also mentioned in \code{\link{meta2d}}),
  setting \code{minper} and \code{maxper} to a same value may be the only
  known way. For short time-series profiles(eg. 10 time points or less),
  p-values given by Lomb-Scargle may be over conservative, which will also
  lead to conservative integrated p-values. In such case, selecting a
  proper p-value cut-off based on the p-value distribution is suggested.
}
\examples{
# write 'cycHumanBloodData' and 'cycHumanBloodDesign' into two 'csv' files
write.csv(cycHumanBloodData, file="cycHumanBloodData.csv",
  row.names=FALSE)
write.csv(cycHumanBloodDesign, file="cycHumanBloodDesign.csv",
  row.names=FALSE)

# detect circadian transcripts in studied individuals
meta3d(datafile="cycHumanBloodData.csv",
  designfile="cycHumanBloodDesign.csv", outdir="example", filestyle="csv",
  design_libColm=1, design_subjectColm=2, design_hrColm=4,
  design_groupColm=3)
}
\references{
Glynn E. F., Chen J., and Mushegian  A. R. (2006). Detecting periodic
  patterns in unevenly spaced gene expression time series using
  Lomb-Scargle periodograms. \emph{Bioinformatics}, \bold{22(3)},
  310--316

Fisher, R.A. (1925). \emph{Statistical methods for research workers}.
  Oliver and Boyd (Edinburgh).

Kugler K. G., Mueller L.A., and Graber A. (2010). MADAM - an open source
  toolbox for meta-analysis. \emph{Source Code for Biology and Medicine},
  \bold{5}, 3.
}

