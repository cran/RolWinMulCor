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

 plot_heatmap <- function(inputdata, corcoefs, pvalues, left_win, 
	          righ_win, Rwidthwin, KCASE="", typewidthwin="", 
	          widthwin_1=3, widthwin_N=dim(inputdata)[1], varX="", 
	 	  varY="", rmltrd=TRUE, Scale=TRUE, coltsX=c("black"), 
	          coltsY="blue", CEXLAB=1.15, CEXAXIS=1.05, LWDtsX=1, 
                  LWDtsY=1, NUMLABX=5, parcen=c(0.5,25))  
 { 

 #---------------------------------------------------------------------#
 # plot_heatmps uses the outputs of the functions rolwincor_heatmap and 
 # rolwinmulcor_heatmap to plot as a heat maps the correlation 
 # coefficients and p-values (corrected or not corrected) for the 
 # bi-variate and multi-variate cases. 
 #---------------------------------------------------------------------#
 #########################################################################
 # Parameters (from 1 to 12) that MUST be provided: 
 # 1. inputdata: the same data that were used in the rolwincor_1win function. 
 # 2. corcoefs:  rolwincor_1win's output named "Correlation_coefficients"
 # 3. pvalues:   rolwincor_1win's output named "P_values_corrected" or 
 #               P_values_not_corrected. 
 # 4. left_win:  rolwincor_1win's output named "left_win", which is used 
 #     	  	 in time domain when the cor. coef and p-values are plotted. 
 # 5. righ_win:  rolwincor_1win's output named "righ_win", which is used   
 #     	  	 in time domain when the cor. coef and p-values are plotted. 
 # 6. Rwidthwin: rolwincor_1win's output named "widthwin", which indicates
 # 	     	 window's size to compute the rolling window correlations. 
 # 7. KCASE:     please use KCASE="BIVAR" for the bi-variate case and 
 #	         KCASE="MULVAR" for the multi-variate case.
 # 8. typewidthwin: to plot all the possible window-lengths ("FULL") or 
 # 		    only a band of window-lengths ("PARTIAL"). 
 # 9. widthwin_1: first for the size of the windows (by default is 3) 
 # 	 	  when the option typewidthwin="PARTIAL" is selected. 
 # 10.widthwin_N: last value for the size of the windows (by default is 
 #                dim(inputdata)[1], i.e. (number of datums in inputdata) 
 # 		  when the option typewidthwin="PARTIAL" is selected.
 #########################################################################
 # 11. varX: name of the first variable for the bi-variate case, e.g. "X" 
 #           (please note that "X" is a vector if KCASE="BIVAR" and a matrix 
 #           if KCASE="MULVAR". For the multi-variate case the names for 
 #           "X" (the independent variables) will be defined as: 
 # 	     varX=paste("X1", "X2",..., sep=", "). 
 # 12. varY: name of the second variable for the bi-variate case, e.g. 
 #           "Y". For the multi-variate case "Y" is the dependent variable. 
 #########################################################################
 # 13. rmltrd:  remove (by default is TRUE; please use "FALSE" 
 #              otherwise) linear trends in the data analysed. 
 # 14. Scale:   scale (by default is "TRUE"; please use "FALSE" 
 #   	        otherwise) or "normalize" or "standardize" the data 
 #	        analysed. 
 # 15. coltsX:  the color to be used to plot the first variable (e.g. "X")
 #              for the bi-variate case. For the multi-variate case the 
 #	        colors will be defined as: coltsX=c("red","blue",...), 
 #              by the default the color is "black". 	
 # 16. coltsY:  the color to be used to plot the second or dependent 
 #  	        variable "Y", by default the color is "blue".  
 # 17. CEXLAB:  the size (cex) for the labels, the default value is 1.15. 
 # 18. CEXAXIS: the size (cex.axis) for the axis, the default value is 1.05. 
 # 19. LWDtsX:  the line width/s (lwd) for the variable/s "X". Please don't 
 #              forget to use LWDtsX = c(1,1,2,...)) for the multi-variate 
 #		case, it has a default value of 1.  
 # 20. LWDtsY:  the line width (lwd) for variable Y, by default is 1. 
 # 21. NUMLABX: this parameter is only used for the multi-variate case,
 #	        and indicates the number of labels for the X's axis
 #		(length.out), it has a default value of 5. 
 # 22. parcen:  this parameter is only used for the multi-variate case,
 # 	        and contains two parameters to control the position of 
 #		the title (by in seq, please look at R>?seq), fist value 
 #		(by default is 0.5) is used to center the title and the 
 #		other (by default is 25) to define the spaces between 
 # 		the names of the variables (please look at R>?mtext).
 #---------------------------------------------------------------------#

 #---------------------------------------------------------------------#
 # Auxiliary code to be used in the plots! 
 # Some set ups for plotting! 
  if (typewidthwin == "FULL") { 
   rango     <- range(Rwidthwin)
   NoYaxis   <- floor(length(Rwidthwin)/5)  
   } 
  if (typewidthwin == "PARTIAL") { 
   nwin      <- length(seq(widthwin_1, widthwin_N, 2)) 
   Rwidthwin <- seq(widthwin_1, widthwin_N, 2) # length(Rwidthwin) = nwin 
  if (nwin <= 10) { 
   rango     <- range(Rwidthwin)
   NoYaxis   <- 2 #floor(length(Rwidthwin)/3)
   }
  if (nwin > 10) { 
   rango     <- range(Rwidthwin)
   NoYaxis   <- floor(length(Rwidthwin)/5) 
   }  
  }   

  # Check 1: removing linear trend - if rmltrd=TRUE
# if(isTRUE(rmltrd)) {
#  ts1.tmp   <- cbind(inputdata[ ,1], c(pracma::detrend(inputdata[ ,2])))
#  ts2.tmp   <- cbind(inputdata[ ,2], c(pracma::detrend(inputdata[ ,3])))
#  inputdata <- cbind(ts1.tmp[ ,1], ts1.tmp[,2], ts2.tmp[ ,2])  
# } 
 if(isTRUE(rmltrd)) {
  dtrd_tmp  <- pracma::detrend(inputdata[,-1])
  inputdata <- cbind(inputdata[,1], dtrd_tmp) 
  } 

 # Check 2: scaling data: [X_i - mean(X_i)]/sd(X-i), mean=0 & sd=1
  if(isTRUE(Scale)) { 
  inputdata <- cbind(inputdata[,1], apply(inputdata[,-1], 2, scale))
 }

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
  the_matrixCOR <- corcoefs
  the_matrixPVA <- pvalues # p-values corrected or not corrected 

  # Setting up the graphical parameters
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar=c(2.35, 4.95, 2.2, 3.5) + 0.1) 
  layout(matrix(c(1,2,3), 3, 1, byrow=FALSE), heights=c(2, 3.5, 0.65))
  # Palette!
  Ncol     <- length(Rwidthwin)
  # Number of colors MUST be length(Rwidthwin) or Ncol 
  Palette  <- colorspace::diverge_hcl(4*Ncol, c=c(100,0), l=c(50,90), power=1.3)
  # Colorbar! 
  rangev   <- seq(min(the_matrixCOR, na.rm=TRUE), max(the_matrixCOR, 
                  na.rm=TRUE), length.out=11)
  rangebar <- matrix(rangev, nrow=1, ncol=11, byrow=TRUE)
  # to make transparent a color (e.g. green), to be used for pvalues > 0.05 
  makeTransparent <- function(someColor, alpha=0.15) scales::alpha(someColor, alpha)
 #---------------------------------------------------------------------#

 #---------------------------------------------------------------------
 # Bi-variate case: rolwincor_1win
 if (KCASE == "BIVAR") { 
 # -------------------------------------------------------------
 # Plot data (plot 1)
 plot(ts1, t="l", col=coltsX, las=1, xlab="", ylab="", xaxt="n", 
  yaxt="n", xaxs="i", yaxs="i", xlim=c(time(ts1)[1], time(ts1)[NL]), 
  main=paste(varX, " vs. ", varY, sep=""), lwd=LWDtsX)
 points(ts2, t="l", col=coltsY, las=1, xlab="", ylab="", lwd=LWDtsY)
 axis(1, at=floor(seq(time.runcor[1], time.runcor[Nrun], length.out=NUMLABX)),
  labels=floor(seq(time.runcor[1], time.runcor[Nrun], length.out=NUMLABX)),
  cex.axis=CEXAXIS)
 axis(2, at=pretty(ts1), labels=pretty(ts1), col=coltsX, las=1, 
  cex.axis=CEXAXIS)
 axis(4, at=pretty(ts2), labels=pretty(ts2), col=coltsY, las=1, 
  cex.axis=CEXAXIS, cex.lab=CEXLAB)
 mtext(1, text="Time", line=2.75, cex=CEXLAB)
 mtext(2, text=varX, col=coltsX, line=2.5, cex=CEXLAB, las=3)
 mtext(4, text=varY, col=coltsY, line=2.5, cex=CEXLAB, las=3)
 # -------------------------------------------------------------
 # Multiscale Window correlation (plot 2)
 # To take into account the statistical significance in the image-plot! 
 new_the_matrixCOR <- the_matrixCOR
 id.xy <- which(the_matrixPVA <= 0.05, arr.ind=TRUE) 
 if (length(id.xy) > 0) {
 for (k in 1:dim(id.xy)[1]) { 
  new_the_matrixCOR[id.xy[k,1], id.xy[k,2]] <- NA
 } 
 } 
 if (length(id.xy) == 0) {
 cat("W A R N I N G: there is no correlation coefficients statistically
 significant. \n") 
 }
 image(t(the_matrixCOR), xlab="", ylab="", las=1, 
  col=Palette, xaxt="n", yaxt="n")
 image(t(new_the_matrixCOR), xlab="", ylab="", las=1, 
  col=makeTransparent("green"), xaxt="n", yaxt="n", add=TRUE)
 contour(t(the_matrixCOR), add=TRUE, drawlabels=TRUE)
 axis(1, at=seq(0, 1, length.out=NUMLABX), labels=floor(seq(time.runcor[1], 
  time.runcor[Nrun], length.out=NUMLABX)), cex.axis=CEXAXIS)
 # Set up Y axis 
 # To verify if NoYaxis is odd or even number! 
 if(NoYaxis %% 2 == 0) { 
 # Case 1: 
 NoYaxis_CASE1 <- NoYaxis  
 to_at_CASE1   <- seq(rango[1], rango[2], by=NoYaxis_CASE1) 
 length_CASE1  <- length(to_at_CASE1)
 } else {
 # Case 2: 
 NoYaxis_CASE2 <- NoYaxis - 1 
 to_at_CASE2   <- seq(rango[1], rango[2], by=NoYaxis_CASE2) 
 length_CASE2  <- length(to_at_CASE2)
 }
 # 
 if (NoYaxis %% 2 == 0) { 
 axis(2, at=seq(0, 1, length.out=length_CASE1), 
  labels=to_at_CASE1, cex.axis=0.95*CEXAXIS, las=1)
 } else { 
   axis(2, at=seq(0, 1, length.out=length_CASE2), 
    labels=to_at_CASE2, cex.axis=0.95*CEXAXIS, las=1)
  }
 mtext(1, text="Time", line=2.75, cex=CEXLAB)
 mtext(2, text="Time-scales", line=3.35, cex=CEXLAB)
 # -------------------------------------------------------------
 # Colorbar (plot 3) 
  image(z=t(rangebar), axes=FALSE, col=Palette, frame.plot=TRUE,
   yaxt="n", xaxt="n") 
  axis(1, at=round(seq(0,1,length.out=11),2), labels=round(rangebar, 
   digits=2), cex.axis=CEXAXIS, las=1)
} 

 # Multi-variate case: rolwinmulcor_1win
 if (KCASE == "MULVAR") { 
 # -------------------------------------------------------------
 # To be used in "plot" and "mtext" 
  # Plot time series (plot 1)
  miny <- min(apply(Z, FUN=min, 2))
  maxy <- max(apply(Z, FUN=max, 2))
  plot(ts1, t="l", col=coltsY, las=1, xlab="", ylab="", xaxt="n", 
   yaxt="n", xaxs="i", yaxs="i", xlim=c(time(ts1)[1], time(ts1)[NL]), 
   ylim=c(1.10*miny, 1.10*maxy), lwd=LWDtsY)
 # To plot the names of variables with their corresponding colors
  labnam   <- c(paste(varY, " <-   ", sep=""), varX)
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
 mtext(1, text="Time", line=2.75, cex=CEXLAB)
 # -------------------------------------------------------------
 # Multiscale Window correlation (plot 2)
 # To take into account the statistical significance in the image-plot! 
 new_the_matrixCOR <- the_matrixCOR
 id.xy <- which(the_matrixPVA <= 0.05, arr.ind=TRUE) 
 if (length(id.xy) > 0) {
 for (k in 1:dim(id.xy)[1]) { 
  new_the_matrixCOR[id.xy[k,1], id.xy[k,2]] <- NA
 } 
 } 
 if (length(id.xy) == 0) {
 cat("W A R N I N G: there is no correlation coefficients statistically
 significant. \n") 
 }
 image(t(the_matrixCOR), xlab="", ylab="", las=1, 
  col=Palette, xaxt="n", yaxt="n")
 image(t(new_the_matrixCOR), xlab="", ylab="", las=1, 
  col=makeTransparent("green"), xaxt="n", yaxt="n", add=TRUE)
 contour(t(the_matrixCOR), add=TRUE, drawlabels=TRUE)
 axis(1, at=seq(0, 1, length.out=NUMLABX), labels=floor(seq(time.runcor[1], 
  time.runcor[Nrun], length.out=NUMLABX)), cex.axis=CEXAXIS)
 # Set up Y axis 
 # To verify if NoYaxis is odd or even number! 
 if(NoYaxis %% 2 == 0) { 
 # Case 1: 
 NoYaxis_CASE1 <- NoYaxis  
 to_at_CASE1   <- seq(rango[1], rango[2], by=NoYaxis_CASE1) 
 length_CASE1  <- length(to_at_CASE1)
 } else {
 # Case 2: 
 NoYaxis_CASE2 <- NoYaxis - 1 
 to_at_CASE2   <- seq(rango[1], rango[2], by=NoYaxis_CASE2) 
 length_CASE2  <- length(to_at_CASE2)
 }
 # 
 if (NoYaxis %% 2 == 0) { 
 axis(2, at=seq(0, 1, length.out=length_CASE1), 
  labels=to_at_CASE1, cex.axis=0.95*CEXAXIS, las=1)
 } else { 
   axis(2, at=seq(0, 1, length.out=length_CASE2), 
    labels=to_at_CASE2, cex.axis=0.95*CEXAXIS, las=1)
  }
 mtext(1, text="Time", line=2.75, cex=CEXLAB)
 mtext(2, text="Time-scales", line=3.35, cex=CEXLAB)
 # -------------------------------------------------------------
 # Colorbar (plot 3) 
 image(z=t(rangebar), axes=FALSE, col=Palette, frame.plot=TRUE,
  yaxt="n", xaxt="n") 
 axis(1, at=round(seq(0,1,length.out=11),2), labels=round(rangebar, 
  digits=2), cex.axis=CEXAXIS, las=1)
 }

 } # End of the function plot_1win
