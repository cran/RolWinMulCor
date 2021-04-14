######################################################################
#:: plot_1win function - R package RolWinMulCor                      #
#:: Programmed by Josué M. Polanco-Martinez a.k.a jomopo             #
#:: Email: josue.m.polanco@gmail.com                                 #
######################################################################
#   Copyright (C) 2020 by Josué M. Polanco-Martínez 	             #
#   This file/code is part of the R package RolWinMulCor             #
######################################################################
#								     
#   RolWinMulCor (Rolling Window Multiple Correlation) is free software: 
#   you can redistribute it and/or modify it under the terms of the GNU 
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

 plot_1win <- function(inputdata, corcoefs, pvalues, left_win, righ_win, 
	       widthwin, KCASE="", varX="", varY="", coltsX=c("black"), 
               coltsY="blue", rmltrd=TRUE, Scale=TRUE, HeigWin1=2.05, 
               HeigWin2=2.75, colCOEF="black", colPVAL="gray", 
               CEXLAB=1.15, CEXAXIS=1.05, LWDtsX=c(1), LWDtsY=1, 
               LWDcoef=1, LWDpval=1, NUMLABX=5, parcen=c(0.5,25))
 { 

 #---------------------------------------------------------------------#
 # plot_1win uses the outputs of the functions rolwincor_1win and 
 # rolwinmulcor_1win to plot in ONLY ONE window-length the correlation 
 # coefficients and p-values (corrected or not corrected) for the 
 # bi-variate and multi-variate cases. 
 #---------------------------------------------------------------------#
 #########################################################################
 # Parameters (from 1 to 11) that MUST be provided: 
 # 1. inputdata: the same data that were used in the rolwincor_1win function. 
 # 2. corcoefs:  rolwincor_1win's output named "Correlation_coefficients"
 # 3. pvalues:   rolwincor_1win's output named "P_values_corrected" or 
 #               P_values_not_corrected. 
 # 4. left_win:  rolwincor_1win's output named "left_win", which is used 
 #     	  	 in time domain when the cor. coef and p-values are plotted. 
 # 5. righ_win:  rolwincor_1win's output named "righ_win", which is used   
 #     	  	 in time domain when the cor. coef and p-values are plotted. 
 # 6. widthwin:  rolwincor_1win's output named "widthwin", which indicates
 # 	     	 window's size to compute the rolling window correlations. 
 # 7. KCASE:     please use KCASE="BIVAR" for the bi-variate case and 
 #	         KCASE="MULVAR" for the multi-variate case.
 # 8. varX:  name of the first variable for the bi-variate case, e.g. "X" 
 #           (please note that "X" is a vector if KCASE="BIVAR" and a matrix 
 #           if KCASE="MULVAR". For the multi-variate case the names for 
 #           "X" (the independent variables) will be defined as: 
 # 	     varX=paste("X1", "X2",..., sep=", "). 
 # 9. varY:  name of the second variable for the bi-variate case, e.g. 
 #           "Y". For the multi-variate case "Y" is the dependent variable. 
 # 10. coltsX: the color to be used to plot the first variable (e.g. "X")
 #             for the bi-variate case. For the multi-variate case the 
 #	       colors will be defined as: coltsX=c("red","blue",...), 
 #             by the default the color is "black". 	
 # 11. coltsY: the color to be used to plot the second or dependent 
 #  	       variable "Y", by default the color is "blue".  
 #########################################################################

 # 12. rmltrd:  remove (by default is TRUE; please use "FALSE" 
 #              otherwise) linear trends in the data analysed. 
 # 13. Scale:   scale (by default is "TRUE"; please use "FALSE" 
 #   	        otherwise) or "normalize" or "standardize" the data 
 #	        analysed. 
 # 14. HeigWin1: proportion of window to plot the variables (>R?layout).
 # 15. HeigWin2: proportion of window to plot cor. coef. and p-values. 
 # 16. colCOEF: the color to be used to plot the correlation coefficients, 
 #		by default the color is "black". 
 # 17. colPVAL: the color to be used to plot the p-values, by the default
 #	        the color is "gray". 
 # 18. CEXLAB:  the size (cex) for the labels, the default value is 1.15. 
 # 19. CEXAXIS: the size (cex.axis) for the axis, the default value is 1.05. 
 # 20. LWDtsX:  the line width/s (lwd) for the variable/s "X". Please don't 
 #              forget to use LWDtsX = c(1,1,2,...)) for the multi-variate 
 #		case, it has a default value of 1. 
 # 21. LWDtsY:  the line width (lwd) for variable Y, by default is 1. 
 # 22. LWDcoef: the line width (lwd) to be used to plot the correlation 
 #    	        coefficients, it has a default value of 1. 
 # 23. LWDpval: the line width (lwd) to be used to plot p-values, the 
 # 		default value is 1. 
 # 24. NUMLABX: this parameter is only used for the multi-variate case,
 #	        and indicates the number of labels for the X's axis
 #		(length.out), it has a default value of 5. 
 # 25. parcen:  this parameter is only used for the multi-variate case,
 # 	        and contains two parameters to control the position of 
 #		the title (by in seq, please look at R>?seq), fist value 
 #		(by default is 0.5) is used to center the title and the 
 #		other (by default is 25) to define the spaces between 
 # 		the names of the variables (please look at R>?mtext).
 #---------------------------------------------------------------------#

 #---------------------------------------------------------------------#
  # Check 1: removing linear trend - if rmltrd=TRUE
  if(isTRUE(rmltrd)) {
   dtrd_tmp  <- pracma::detrend(inputdata[,-1])
   inputdata <- cbind(inputdata[,1], dtrd_tmp) 
  } 

  # Check 2: scaling data: [X_i - mean(X_i)]/sd(X-i), mean=0 & sd=1
  if(isTRUE(Scale)) { 
   inputdata <- cbind(inputdata[,1], apply(inputdata[,-1], 2, scale))
  }

 # Auxiliary code to be used in the plots! 
  # Transforming input data to time series objects! 
  #Deltat <- diff(inputdata[,1]) 
  NL     <- dim(inputdata)[1]
  NP     <- dim(inputdata)[2]
  freq   <- length((inputdata[,1]))/length(unique(inputdata[,1]))
  ts1    <- ts(inputdata[,2], start=c(inputdata[1,1],1), 
            end=c(inputdata[NL,1],freq), deltat=1/freq) # ts1 = Y
  ts2    <- ts(inputdata[,3:NP], start=c(inputdata[1,1],1), 
	    end=c(inputdata[NL,1],freq), deltat=1/freq) #ts2 = X (X1, X2, ..., Xp) 
  Z      <- cbind(ts1, ts2) # To be used in the running cor./fun. (rollapply)
  #
  time.runcor <- time(ts1)
  Nrun        <- length(time.runcor) 
  #
  rc.ts1_ts2         <- corcoefs[,2]
  CORTD_pval_ts1_ts2 <- pvalues[,2] # p-values corrected or not corrected 

  # Setting up the graphical parameters
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar=c(3.90, 4.05, 2.2, 3.3) + 0.1) 
  layout(matrix(c(1,2), 2, 1, byrow=FALSE), heights=c(HeigWin1, HeigWin2))
 #---------------------------------------------------------------------#

 #---------------------------------------------------------------------#
 # Bi-variate case: rolwincor_1win
 if (KCASE == "BIVAR") { 
 ####################################################################
 # Setting up the graphical parameters
 ####################################################################
 #### pending to use open par and close, it's very important for CRAN
 ####################################################################
 # Plot data 
 plot(ts1, t="l", col=coltsX, las=1, xlab="", ylab="", xaxt="n", 
  yaxt="n", xaxs="i", yaxs="i", xlim=c(time(ts1)[1], time(ts1)[NL]), 
  main=paste(varX, " vs. ", varY, sep=""), lwd=LWDtsX)
 points(ts2, t="l", col=coltsY, las=1, xlab="", ylab="", lwd=LWDtsY)
 axis(1, at=pretty(time.runcor), labels=pretty(time.runcor), 
  cex.axis=CEXAXIS)
 axis(2, at=pretty(ts1), labels=pretty(ts1), col=coltsX, las=1, 
  cex.axis=CEXAXIS)
 axis(4, at=pretty(ts2), labels=pretty(ts2), col=coltsY, las=1, 
  cex.axis=CEXAXIS, cex.lab=CEXLAB)
 mtext(1, text="Time", line=2.5, cex=CEXLAB)
 mtext(2, text=varX, col=coltsX, line=2, cex=CEXLAB, las=3)
 mtext(4, text=varY, col=coltsY, line=2, cex=CEXLAB, las=3)
 # Plot coef. cor. values 
 plot(time.runcor[left_win:(Nrun-righ_win)], rc.ts1_ts2, t="l", yaxt="n", 
  xaxs="i", yaxs="i", xlim=c(range(time.runcor)[1], range(time.runcor)[2]),  
  ylim=c(-1,1), cex.lab=CEXLAB, cex.axis=CEXAXIS, lwd=LWDcoef,  
  xlab="Time", ylab="Dynamic correlation coefficient", las=1, 
   col="black", main=paste(varX, " vs. ", 
   varY, " (", "window size = ", widthwin, ")", sep="")) 
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
 } 

 # Multi-variate case: rolwinmulcor_1win
 if (KCASE == "MULVAR") { 
 ####################################################################
 # Setting up the graphical parameters
 ####################################################################
 #### pending to use open par and close, it's very important for CRAN
 ####################################################################
# Plot data (plot 1) 
 # To be used in "plot" and "mtext"
 miny <- min(apply(Z, FUN=min, 2))
 maxy <- max(apply(Z, FUN=max, 2))
 strspl_vnam <- unlist(strsplit(varX, ","))
 compos_vnam <- paste(strspl_vnam, rep("(", length(strspl_vnam)), 
                      coltsX, rep(")", length(strspl_vnam)), sep="", 
		      collapse = ", ") 
 # Plot the time series 
 plot(ts1, t="l", col=coltsY, las=1, xlab="", ylab="", xaxt="n", 
  yaxt="n", xaxs="i", yaxs="i", xlim=c(time(ts1)[1], time(ts1)[NL]), 
  ylim=c(1.10*miny, 1.10*maxy), lwd=LWDtsY)

  labnam   <- c(paste(varY, " <-   ", sep=""), strspl_vnam)
  colrs    <- c(coltsY, coltsX)

  firtX  <- nchar(varX[1])
  bytim  <- cumsum(nchar(c(paste(varY, " <-  ", sep=""), varX)))
  schar  <- sum(nchar(c(paste(varY, " <-  ", sep=""), varX)))
  at_tmp <- round((time(ts1)[Nrun] - time(ts1)[1])/round(schar/(length(labnam)))/parcen[1])
  at_times <- seq(time(ts1)[at_tmp], length.out=dim(Z)[2], by=parcen[2]) 

 for(i in 1:dim(Z)[2]) { 
  if(i == 1) 
   mtext(at=at_times[i], labnam[i], col=colrs[i], adj=0.5, font=2) 
  if (i > 1 & i < dim(Z)[2])
   mtext(at=at_times[i], paste(labnam[i], ",", sep=""), col=colrs[i], adj=0.5, font=2) 
  if (i == dim(Z)[2]) 
    mtext(at=at_times[i], labnam[i], col=colrs[i], adj=0.5, font=2)
  }

 for(j in 1:(NP-2)){
 points(ts2[,j], t="l", col=coltsX[j], las=1, xlab="", ylab="", 
  lwd=LWDtsX[j])
 }
 axis(1, at=floor(seq(time.runcor[1], time.runcor[Nrun], length.out=NUMLABX)),
  labels=floor(seq(time.runcor[1], time.runcor[Nrun], length.out=NUMLABX)),
  cex.axis=CEXAXIS)
 axis(2, at=pretty(ts1), labels=pretty(ts1), col=coltsY, las=1, 
  cex.axis=CEXAXIS)
 mtext(1, text="Time", line=2.5, cex=CEXLAB)
 # Plot coef. cor. values (plot 2) 
 plot(time.runcor[left_win:(Nrun-righ_win)], rc.ts1_ts2, t="l", xaxt="n", yaxt="n", 
  xaxs="i", yaxs="i", xlim=c(range(time.runcor)[1], range(time.runcor)[2]),  
  ylim=c(-1,1), cex.lab=CEXLAB, cex.axis=CEXAXIS, lwd=LWDcoef,  
  xlab="Time", ylab="Dynamic correlation coefficient", las=1, 
   col="black", main=paste(varY, " <- ", varX, " (", 
   "window size = ", widthwin, ")", sep="")) 
 axis(2, seq(-1, 1, 0.20), col=colPVAL, las=1, cex.axis=CEXAXIS)
 abline(h=0, col=1, lty=1)
 par(new=T)
 plot(time.runcor[left_win:(Nrun-righ_win)], CORTD_pval_ts1_ts2, t="l", 
  col=colPVAL, xlab="", ylab="", xaxs="i", yaxs="i", lwd=LWDpval,  
  xlim=c(range(time.runcor)[1], range(time.runcor)[2]), ylim=c(0,1), axes=FALSE) 
 axis(1, at=floor(seq(time.runcor[1], time.runcor[Nrun], length.out=NUMLABX)),
  labels=floor(seq(time.runcor[1], time.runcor[Nrun], length.out=NUMLABX)),
  cex.axis=CEXAXIS)
 axis(4, seq(0,1,0.1), col=colPVAL, las=1, cex.axis=CEXAXIS)
 mtext(4, text="P-value", col=colPVAL, line=2.25, las=3, cex=CEXLAB)
 abline(h=0.05, col=colPVAL, lty=1)
 abline(h=0.10, col=colPVAL, lty=2)
 }
 
 } # End of the function plot_1win
