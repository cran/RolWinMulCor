######################################################################
#:: rolwinmulcor_1win function - R package RolWinMulCor              #
#:: Multi-variate case! 
#:: Programmed by Josué M. Polanco-Martinez a.k.a jomopo             #
#:: Email: josue.m.polanco@gmail.com                                 #
######################################################################
#   Copyright (C) 2020 by Josué M. Polanco-Martínez 	             #
#   This file/code is part of the R package RolWinMulCor             #
######################################################################
#								     
#   RolWinMulCor (Rolling Window Multiple Correlation) is free software: 
#   otu can redistribute it and/or modify it under the terms of the GNU 
#   General Public License as published by the Free Software 
#   Foundation, either version 3 of the License, or (at your option) 
#   any later version.
#
#   RolWinMulCor is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with RolWinMulCor If not, see <http://www.gnu.org/licenses/>.
#
#####################################################################

 rolwinmulcor_1win <- function(inputdata, varnametsY="", 
              varnametsX="", units="", coltsY="black", 
              coltsX="", colCOEF="black", colPVAL="gray", 
              CEXLAB=1.15, CEXAXIS=1.05, LWDtsY=1, LWDtsX="", 
              LWDcoef=1, LWDpval=1, CorMethod="spearman", 
	      widthwin=3, Align="center", pvalcorectmethod="BH", 
              rmltrd="Y", Scale="Y", device="screen", Hfig=1200, 
	      Wfig=900, resfig=150, Hpdf=5.15, Wpdf=13.5, 
	      HeigWin1=2.05, HeigWin2=2.75, ofilename) { 

 #---------------------------------------------------------------------#
 # rolwinmulcor_1win function estimates the rolling (running) 
 # window correlation among multiple time series (multi-variate case) 
 # sampled on identical time points for ONLY ONE window-length and plots 
 # the correlation coefficients and their respective p-values.
 #---------------------------------------------------------------------#
 # Description of variables (INPUTS): 
 # 1.  inputdata:  matrix of P columns (time, dependent variable ("Y"), 
 #                 independent variables ("X1", "X2",... "Xp-1")
 # 2.  varnametsY: name of the dependent variable: "Y", please note that 
 #                 you MUST be define (the name) this variable
 # 3.  varnametsX: name of the independent variables "X1", "X2", ... please
 #  		   note that you MUST define (the names) these variables 
 # 4.  units:      time's unit (e.g. days, weeks, years, etc), you MUST 
 #	   	   define the units. 
 # 5.  coltsY:     the color to be used to plot the dependent variable "Y", 
 # 		   by default the color is "black"
 # 6.  coltsX:     the color to be used to plot the independent variables 
 #                 ("X1", "X2",... "Xp-1"), please note that you MUST 
 # 		   supply this information (e.g. coltsX=c("red","blue",...))
 # 7.  colCOEF:    the color to be used to plot correlation coefficients 
 # 8.  colPVAL:    the color to be used to plot p-values 
 # 9.  CEXLAB:     size for labels
 # 10. CEXAXIS:    size for axis 
 # 11. LWDtsY:     the line-width for the dependent variable ("Y"), by 
 #		   default it has a value of 1, but other values can be used  
 # 12. LWDtsX:     the line width for the independent variables ("X1", 
 #		   "X2",... "Xp-1"), please note that you MUST supply this 
 #		   information (e.g. LWDtsX = c(1,1,2,...))
 # 13. LWDcoef:    the line width to be used to plot cor. coef. 
 # 14. LWDpval:    the line width to be used to plot p-values 
 # 15. CorMethod:  method used to estimate the correlations, by default 
 #                 is "spearman" but other options ("pearson" and "kendall")
 #                 are possible (please see at in R>?cor.test)
 # 16. widthwin:   window's size to compute the rolling window correlations, 
 #   		   this MUST be odd (i.e. of the form 2m + 1)
 # 17. Align:      to align the rolling object, RolWinMulCor ONLY uses 
 #                 the "center" option (please look at in R>?running) 
 # 18. pvalcorectmethod: p value correction method (please look at: 
 #                 R>?p.adjust), by default we use the method of Benjamini 
 # 		   & Hochberg (BH), but other 6 methods can be used
 # 19. rmltrd:  remove (by default is "Y" or "y"; please use "N" or "n" 
 # 	        otherwise) linear trend in the data analysed  
 # 20. Scale:   scale  (by default is "Y" or "y"; please use "N" or "n" 
 # 		otherwise) or "normalized" or "standardized" the data analysed
 # 21. device:  the kind of output plot (please look at: R?device), 
 #              you can use: "png", "jpeg/jpg", "eps", "pdf" and "screen" 
 # 22. Hfig:    plot's height (device) in "png" and "jpg" format (>R?pdf)
 # 23. Wfig:    plot's width (device) in "png" and "jpg" format (>R?postscript)
 # 24. resfig:  resolution for the plot in "png" and "jpg" format (>R?jpeg)
 # 25. Hpdf:    plot's height (device) in "pdf" or "eps" format (>R?pdf)
 # 26. Wpdf:    plot's width (device) in "pdf" or "eps" format (>R?postscript)
 # 27. HeigWin1:   proportion of window to plot the variables (>R?layout)
 # 28: HeigWin2:   proportion of window to plot cor. coef. and p-values
 # 29. ofilename:  output file name 
 #------------------------------------------------------------------------#

 # To reset the par options modified within this function!
 oldpar <- par(no.readonly = TRUE)
 on.exit(par(oldpar))   

 #:: Devices options: png, jpg, eps, pdf & screen! 
 path <- tempdir()
 if (device=="png") {
  fileout <- file.path(path, paste("plot_multivariate_", ofilename, 
              ".png", sep=""))
  png(fileout, height=Hfig, width=Wfig, res=resfig) 
  par(mar=c(3.90, 4.05, 2.2, 3.3) + 0.1) 
  layout(matrix(c(1,2), 2, 1, byrow=FALSE), heights=c(HeigWin1, HeigWin2))
 }

 if (device=="jpeg" || device=="jpg") {
  fileout <- file.path(path, paste("plot_multivariate_", ofilename, 
              ".jpg", sep=""))
  jpeg(fileout, height=Hfig, width=Wfig, res=resfig)
  par(mar=c(3.90, 4.05, 2.2, 3.3) + 0.1) 
  layout(matrix(c(1,2), 2, 1, byrow=FALSE), heights=c(HeigWin1, HeigWin2))
 }

 if (device=="pdf") {
  fileout <- file.path(path, paste("plot_multivariate_", ofilename, 
	      ".pdf", sep=""))
  pdf(fileout, height=Hpdf, width=Wpdf)
  par(mar=c(3.90, 4.05, 2.2, 3.3) + 0.1) 
  layout(matrix(c(1,2), 2, 1, byrow=FALSE), heights=c(HeigWin1, HeigWin2))
 }

 if (device=="eps") {
  fileout <- file.path(path, paste("plot_multivariate_", ofilename, 
	      ".eps", sep="")) 
  postscript(file=fileout, height=Hpdf, width=Wpdf)
  par(mar=c(3.90, 4.05, 2.2, 3.3) + 0.1) 
  layout(matrix(c(1,2), 2, 1, byrow=FALSE), heights=c(HeigWin1, HeigWin2))
 }
 
  # ----------------------------------------------------------------------- 
  # Check 1: number of columns, inputdata MUST contain at least four
  #          columns: time, dependent variable Y, and two (or more)
  # 	     independent variables X1, X2,...  
  if (dim(inputdata)[2] < 4) 
   stop("The input MUST be an array/matrix with at least 4 columns 
    (first column the time, second the dependent variable (Y) and 
    the others the independent variables (X1, X2,...). Thank you for 
    using RolWinMulCor package.)")

  # Check 2: the time steps MUST be regular/evenly - no gaps! 
  Deltat <- diff(inputdata[,1])  # Deltat is the temporal resolution! 
  if (length(unique(Deltat)) != 1)
   stop("The data have some gaps (it's irregular), please, consider 
    to face this before use RolWinMulCor. Thank you for using 
    RolWinCor package.")

  # Check 3: widthwin MUST be odd (of the form 2m + 1, where m is 1, 2,...)
  if (widthwin %% 2 == 0) { 
   stop(paste(widthwin, "is EVEN and this MUST be ODD. Thank you 
    for using RolWinMulCor package."))
  }

  # Check 4: removing linear trend - if rmltrd="Y" 
  if (rmltrd == "Y" || rmltrd == "y") {
   dtrd_tmp  <- pracma::detrend(inputdata[,-1])
   inputdata <- cbind(inputdata[,1], dtrd_tmp) 
  } 

  # Check 5: scaling data: [X_i - mean(X_i)]/sd(X-i), mean=0 & sd=1
  if (Scale == "Y" || Scale == "y") { 
   inputdata <- cbind(inputdata[,1], apply(inputdata[,-1], 2, scale))
  }
  # ----------------------------------------------------------------------- 

  # Transforming input data to time series object! 
  NL  <- dim(inputdata)[1]
  NP  <- dim(inputdata)[2]
  ts1 <- ts(inputdata[,2], start=inputdata[1,1], end=inputdata[NL,1], 
     	    deltat=unique(Deltat)) # ts1 = Y
  ts2 <- ts(inputdata[,3:NP], start=inputdata[1,1], end=inputdata[NL,1], 
	    deltat=unique(Deltat)) #ts2 = X (X1, X2, ..., Xp) 

  time.runcor <- time(ts1)
  Nrun        <- length(time.runcor) 

  # ------------------------------------------------------------------
  # Computing the rolling window (running) correlation and p-values  

  #############	        Auxiliary   function              ############
  #Function to estimate the squared multiple correlation coefficient 
  cor.def  <- function(Y, X) { 
   # where X (dependent variables) is a matrix with P columns and M rows! 
   lm_estimate <- lm(Y ~ X)
   summary_lm  <- summary(lm_estimate) 
   adjRsqu     <- summary_lm$adj.r.squared
   Fstat       <- summary_lm$fstatistic[1]
   pvalueF     <- pf(summary_lm$fstatistic[1], summary_lm$fstatistic[2], 
                    summary_lm$fstatistic[3], lower.tail=F)
   namesLS     <- c("Squared_multiple_cor_coef", "F-statistics", "P-value")
   LIST        <- list(adjRsqu, Fstat, pvalueF) 
  return(LIST)
  }

  rc.ts1_ts2_tmp <- gtools::running(ts1, ts2, fun=cor.def,
		     width=widthwin, align=Align)
  rc.ts1_ts2     <- rc.ts1_ts2_tmp[1,]
  pvalrc.ts1_ts2 <- rc.ts1_ts2_tmp[3,]

  # p-value correction (p.adjust)
  ncompa 	      <- length(rc.ts1_ts2) 
  CORTD_pval_ts1_ts2  <- round(p.adjust(pvalrc.ts1_ts2, 
                          method=pvalcorectmethod, n=ncompa), 4)
  ##### 
  left_win <-  (widthwin - 1)/2 + 1 
  righ_win <-  (widthwin - 1)/2 

 # Setting up the graphical parameters
 if (device == "screen") { 
  par(mar=c(3.90, 4.05, 2.2, 3.3) + 0.1) 
  layout(matrix(c(1,2), 2, 1, byrow=FALSE), heights=c(HeigWin1, HeigWin2))
 }
 # Plot data 
 # To be used in "main"
 strspl_vnam <- unlist(strsplit(varnametsX, ","))
 compos_vnam <- paste(strspl_vnam, rep("(", length(strspl_vnam)), 
                      coltsX, rep(")", length(strspl_vnam)), sep="", 
		      collapse = ", ") 
 plot(ts1, t="l", col=coltsY, las=1, xlab="", ylab="", xaxt="n", 
  yaxt="n", xaxs="i", yaxs="i", xlim=c(time(ts1)[1], time(ts1)[NL]), 
  main=paste(varnametsY, "(", coltsY, ")", "<-", compos_vnam, sep=""), 
  lwd=LWDtsY)
 for(j in 1:(NP-2)){
 points(ts2[,j], t="l", col=coltsX[j], las=1, xlab="", ylab="", 
  lwd=LWDtsX[j])
 }
 axis(1, at=pretty(time.runcor), labels=pretty(time.runcor), 
  cex.axis=CEXAXIS)
 axis(2, at=pretty(ts1), labels=pretty(ts1), col=coltsY, las=1, 
  cex.axis=CEXAXIS)
 mtext(1, text="Time", line=2.5, cex=CEXLAB)
 # Plot coef. cor. values 
 plot(time.runcor[left_win:(Nrun-righ_win)], rc.ts1_ts2, t="l", yaxt="n", 
  xaxs="i", yaxs="i", xlim=c(range(time.runcor)[1], range(time.runcor)[2]),  
  ylim=c(-1,1), cex.lab=CEXLAB, cex.axis=CEXAXIS, lwd=LWDcoef,  
  xlab="Time", ylab="Dynamic correlation coefficient", las=1, 
   col="black", main=paste(varnametsY, " <- ", 
   varnametsX, " (", "window size = ", widthwin, " ", units, ")", sep="")) 
 axis(2, seq(-1, 1, 0.20), col=colPVAL, las=1, cex.axis=CEXAXIS)
 abline(h=0, col=1, lty=1)
 par(new=T)
 plot(time.runcor[left_win:(Nrun-righ_win)], CORTD_pval_ts1_ts2, t="l", 
  col=colPVAL, xlab="", ylab="", xaxs="i", yaxs="i", lwd=LWDpval,  
  xlim=c(range(time.runcor)[1], range(time.runcor)[2]), ylim=c(0,1), axes=FALSE) 
 axis(4, seq(0,1,0.1), col=colPVAL, las=1, cex.axis=CEXAXIS)
 mtext(4, text="P-value", col=colPVAL, line=2.25, las=3, cex=CEXLAB)
 abline(h=0.05, col=colPVAL, lty=1)
 abline(h=0.10, col=colPVAL, lty=2)

 if (device != "screen") 

 dev.off()

 # Numerical output 
 namesLS <- c("Correlation_coefficients", "P_values_corrected", 
              "P_values_not_corrected") 
 LIST    <- list(cbind(time.runcor[left_win:(Nrun-righ_win)], 
             rc.ts1_ts2), cbind(time.runcor[left_win:(Nrun-righ_win)], 
             pvalrc.ts1_ts2), cbind(time.runcor[left_win:(Nrun-righ_win)], 
             CORTD_pval_ts1_ts2))
 names(LIST) <- namesLS

 return(LIST)
 
} # End Main function 


