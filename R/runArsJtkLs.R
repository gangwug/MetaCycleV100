###======================================================================================================================================
### This files contains three functions used to analyze expression profiles with ARS, JTK and LS, respectively.
### Author: Gang Wu
### Email: wggucas@gmail.com
### Lab: John Hogenesch's lab in Perelman School of Medicine at University of Pennsylvania (http://hogeneschlab.org/)
###======================================================================================================================================
runARS <- function(indata,ARStime,minper=20,maxper=28, arsper=24, arsmet="")
{
  #-----------------------
  cat("The ARS is in process from ",format(Sys.time(), "%X %m-%d-%Y"),"\n");
  start=minper;
  end=maxper;
  pvalues <-NA
  #-----------------------
  self.data=indata;
  self.delta=ARStime[2]-ARStime[1];                
  time_points=ARStime;
  idorder<-dimnames(self.data)[[1]];
  dataM=as.matrix(self.data[,2:ncol(self.data)]);
  outID=as.character(self.data[,1]);
  names(outID)<-idorder;
  dimnames(dataM)<-list("r"=idorder,"c"=paste("T",self.delta*(0:(ncol(dataM)-1)),sep="" ) );
  #-----------------------
  expSD<-apply(dataM,1,sd);
  expMEAN<-apply(dataM,1,mean);
  constantID<-names(expSD[expSD == 0]);
  flagV<-rep(1,length(idorder));
  names(flagV)<-idorder;
  flagV[constantID]<-0;
  #-----------------------
  run_start=proc.time();
  set.seed(run_start["elapsed"]);
  header<-c("filter_type","ar_method","period_number","period","amplitude","phase","mean","R_square","R2_adjust","coef_var","pvalue");
  ars.outM<<-header;
  ori.op<-options();
  #-----------------------try apply in latter version, it may improve the Compuatational Efficiency
  for (line in idorder )
  {
    if (flagV[line]) {
		d=evaluate(eva.x=time_points,eva.y=dataM[line,],eva.delta=self.delta,eva.pID=line,T_start=start, T_end=end, T_default=arsper,arsmethods=arsmet);
		ars.filter=0;
		if (d$filter)
		{ ars.filter=1; }
		ars.period.num=length(d$period);
		ars.period=paste(d$period,collapse=",");
		ars.amplitude=paste(d$amplitude,collapse=",");
		ars.phase=paste(d$phase,collapse=",");
		ars.out=c(ars.filter,d$armethod,ars.period.num,ars.period,ars.amplitude,ars.phase,mean(dataM[line,]),d$R2,d$R2adj,d$coefvar,d$pvalue);
		pvalues=c(pvalues,d$pvalue);
	} else {
		ars.out=c(rep(NA,4),0,NA,expMEAN[line],rep(NA,3),1);                           
		pvalues=c(pvalues,1);                              #assign it as '1' insead of 'NA' for avoiding error report in 'pi0.est' step
	}
    ars.outM=rbind(ars.outM,ars.out);
  }
  pvalues=pvalues[2:length(pvalues)];
  ars.outM=ars.outM[2:nrow(ars.outM),];
  options(ori.op);
  dimnames(ars.outM)[[1]]=idorder;
  names(pvalues)=idorder;
  #-----------------------
  #header=c("CirID",header,"qvalue","fdr_BH");               ## qvalue is not used in the new version of ARSER
  # pi0=pi0.est(pvalues);                                    #It will report error when analyzing expression profile with equal values among all time points.
  # if (pi0$p0 == 0)                                         #p0 is required for calculating qvalues by qvalue.cal(); p0 will define the largest q-value calculated by qvalue.cal
  # { pi0$p0=0.95; }                                         #If p0 is '0', the calculated q-value is '0' for all pvalues, thus change it to 0.95 is a stategy to reduce false positive?
  # qvalues=qvalue.cal(pvalues[idorder],pi0$p0);
  header=c("CirID",header,"fdr_BH");
  qvalues_BH=p.adjust(pvalues[idorder],"BH");
  ARSoutM=cbind(outID[idorder],ars.outM[idorder,],qvalues_BH);
  dimnames(ARSoutM)<-list("r"=idorder,"c"=header);
  cat("The analysis by ARS is finished at ",format(Sys.time(), "%X %m-%d-%Y"),"\n");
  return(ARSoutM);
}
###======================================================================================================================================
runJTK <- function(indata,JTKtime,minper=20,maxper=28)
{
	cat("The JTK is in process from ",format(Sys.time(), "%X %m-%d-%Y"),"\n");
	perTK<<-"";
	uni_JTKtime<<-unique(JTKtime);	
	freq=uni_JTKtime[2] - uni_JTKtime[1];
	timefactor<<-factor(JTKtime,levels=sort(uni_JTKtime));
	data_endtime<- length(uni_JTKtime)*freq;
	if ( (data_endtime >= maxper) & ( round(maxper/freq) >= 2) ) {
		if (round(minper/freq) >= 2){
			perTK=seq(round(minper/freq),round(maxper/freq),by=1);
		} else {
			perTK=seq(2,round(maxper/freq),by=1);
			#cat(paste("Warning: the input 'minper' is too small for JTK, it was reset as ", 2*freq ,"\n",sep=""));
		}
	} else if ( (data_endtime < maxper) & (data_endtime >= minper) & ( round(data_endtime/freq) >= 2 ) ) {
		if (round(minper/freq) >= 2) {
			perTK=seq(round(minper/freq),round(data_endtime/freq),by=1);
			#cat(paste("Warning: the input 'maxper' is too large for JTK, it was reset as ", data_endtime ,"\n",sep=""));
		} else {
			perTK=seq(2,round(data_endtime/freq),by=1);
			#cat(paste("Warning: the input 'minper' is too small for JTK, it was reset as ", 2*freq ,"\n",sep=""));
		}
	} else {
		stop(paste("The input 'minper' and 'maxper' is out of the range that JTK can detect.\n", "If hope to use JTK for this analysis, please reset the 'minper' and 'maxper' between ", 2*freq, " and ", data_endtime, ".\n", sep=""));
	}
    if (min(perTK)*freq != minper)
    {   cat(paste("Warning: the input 'minper' is not suitable for JTK, it was reset as ", min(perTK)*freq, "\n",sep=""));  }
    if (max(perTK)*freq != maxper)
    {   cat(paste("Warning: the input 'maxper' is not suitable for JTK, it was reset as ", max(perTK)*freq, "\n",sep=""));    }
	options(stringsAsFactors=FALSE);
	data=indata;
	idorder=dimnames(data)[[1]];
	outID<-as.character(data[,1]);           
	names(outID)<-idorder;
	data <- data[,-1];
	#------------------------
	if (!MISSING_VALUE) {	
		jtk_replicates<-summary(timefactor);
		jtkdist(length(uni_JTKtime),as.numeric(jtk_replicates));
		periods <- perTK;                  
		jtk.init(periods,freq) ;          
		flush.console();
		res <- apply(data,1,function(z) {
		  jtkx(z);
		  c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP);
		})
	} else {
		regular_missing<-apply(data,2,function(z) {
							  z[is.nan(z)] <- NA;
							  if ( all(is.na(z)) | all(!is.na(z)) ) {
								return(TRUE);
							  } else {
								return(FALSE);
							  }
						})
		#--------------
		if ( all(regular_missing) ) {
			tepva<-as.numeric(data[1,]);
			tepvb<-rep(1,length(tepva));
			tepvb[is.na(tepva) | is.nan(tepva)]<-0;
			jtk_replicates<-tapply(tepvb,timefactor,sum);
			#--------
			jtkdist(length(uni_JTKtime),as.numeric(jtk_replicates));
			periods <- perTK;                  
			jtk.init(periods,freq) ;          
			flush.console();
			res <- apply(data,1,function(z) {
			  z<-z[(!is.na(z))&(!is.nan(z))];
			  jtkx(z);
			  c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP);
			})
		} else if ( !all(regular_missing) ) {
			res<-apply(data,1,function(z) {
                    tepva<-as.numeric(z);
                    tepvb<-rep(1,length(tepva));
                    tepvb[is.na(tepva) | is.nan(tepva)]<-0;
                    jtk_replicates<-tapply(tepvb,timefactor,sum);
                    #--------
                    jtkdist(length(uni_JTKtime),as.numeric(jtk_replicates));
                    periods <- perTK;                  
                    jtk.init(periods,freq) ;          
                    flush.console();
                    z<-z[(!is.na(z))&(!is.nan(z))];
                    jtkx(z);
                    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP);
			   })
		} else {
			cat("Sorry for this unknown bug in running JTK, please contact the author of 'ciromics'. Thank you.");
		}
	}
	#------------------------
	JTKoutM <- as.data.frame(t(res));
	bhq <- p.adjust(unlist(JTKoutM[,1]),"BH");
	JTKoutM <- cbind(outID,bhq,JTKoutM);	       
	rownames(JTKoutM)<-idorder;
	colnames(JTKoutM)<-c("CirID","BH.Q","ADJ.P","PER","LAG","AMP");
	cat("The analysis by JTK is finished at ",format(Sys.time(), "%X %m-%d-%Y"),"\n");
	return(JTKoutM);
}
###======================================================================================================================================
runLS <- function(indata,LStime,minper=20,maxper=28)
{
	cat("The LS is in process from ",format(Sys.time(), "%X %m-%d-%Y"),"\n");
	RawData <- indata;
	outID <- as.character(RawData[,1]);
	Expression <- data.matrix(RawData[,2:ncol(RawData)]);
	Time <- LStime;
	stopifnot( length(outID)   == nrow(Expression) );                               
	#-----------------------
	N <- ncol(Expression);
	stopifnot(length(Time) == N);
	M <- 4*N;
	MinFrequency <- 1/maxper;
	MaxFrequency <- 1/minper;
	if (MaxFrequency > 1/(2*mean(diff(Time))))
	{  cat("MaxFrequency may be above Nyquist limit.\n"); }
	TestFrequencies <- MinFrequency +
					  (MaxFrequency - MinFrequency) * (0:(M-1) / (M-1))
	#-----------------------
	header<-c("CirID","PhaseShift","PhaseShiftHeight","PeakIndex","PeakSPD","Period","p","N","Nindependent","Nyquist");
	LSoutM<-header;
	for (j in 1:length(outID))                                                     
	{
		if (.Platform$OS.type == "windows")                                       
		{ flush.console(); }                                                       # display immediately in Windows 
		Nindependent <- NHorneBaliunas(sum(!is.na(Expression[j,]))); 
		LS <- ComputeAndPlotLombScargle(Time, Expression[j,], TestFrequencies, Nindependent);
		LSoutM <- rbind(LSoutM,c(outID[j], LS$h.peak$maximum,LS$h.peak$objective, LS$PeakIndex,LS$PeakSPD, 
							  LS$PeakPeriod, LS$PeakPvalue, LS$N, LS$Nindependent,  LS$Nyquist));
	}
	LSoutM<-LSoutM[2:nrow(LSoutM),];
    colnames(LSoutM) <- header;
	bhq <- p.adjust(as.numeric(LSoutM[,"p"]),"BH");
	LSoutM <- cbind(LSoutM,bhq);	       
	dimnames(LSoutM)<-list("r"=1:length(outID),"c"=c(header,"BH.Q"));
	cat("The analysis by LS is finished at ",format(Sys.time(), "%X %m-%d-%Y"),"\n");
	return(LSoutM);
}
###======================================================================================================================================
#End
