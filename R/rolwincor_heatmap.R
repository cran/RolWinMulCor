######################################################################
#:: rolwincor_heatmap function - R package RolWinMulCor              #
#:: bi-variate case!				 		     # 
#:: Programmed by Josué M. Polanco-Martinez a.k.a jomopo             #
#:: Email: josue.m.polanco@gmail.com                                 #
######################################################################
#   Copyright (C) 2020 by Josué M. Polanco-Martínez 	             #
#   This file/code is part of the R package RolWinCor                #
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

 rolwincor_heatmap <- function(inputdata, varnametsX="", 
		varnametsY="", units="", coltsX="black", coltsY="blue", 
		CEXLAB=1.15, CEXAXIS=1.5, LWDtsX=1, LWDtsY=1, 
 	        CorMethod="spearman", typewidthwin="FULL", widthwin_1=3, 
		widthwin_N=dim(inputdata)[1], 
		Align="center", pvalcorectmethod="BH", rmltrd="Y", 
	        Scale="Y", device="screen", Hfig=900, Wfig=900, 
	        resfig=150, Hpdf=7, Wpdf=7, NUMLABX=5, ofilename) {

 #---------------------------------------------------------------------#
 # The rolwincor_heatmap function estimates the rolling (running)
 # window correlation between TWO time series (bi-variate case) 
 # sampled on identical time points for all the possible (i.e. from 3 to 
 # the number of elements of the time series under analysis) window-lengths 
 # or for a band of window-lengths and plots the correlation coefficients 
 # and their respective p-values as a heat map.
 #---------------------------------------------------------------------#
 # Description of variables (INPUTS): 
 # 1.  inputdata:  matrix of 3 columns: time, variable 1 and variable 2 
 # 2.  varnametsX: name of variable 1 
 # 3.  varnametsY: name of variable 2
 # 4.  units:      time's unit (e.g. days, weeks, years, etc) 
 # 5.  coltsX:     the color to be used to plot variable 1 
 # 6.  coltsY:     the color to be used to plot variable 2 
 # 7.  CEXLAB:     size for labels
 # 8.  CEXAXIS:    size for axis 
 # 9.  LWDtsX:     the line-width when variable 1 is plotted 
 # 10. LWDtsY:     the line-width when variable 2 is plotted 
 # 11. CorMethod:  method used to estimate the correlations, by default 
 #                 is "spearman" but other options ("pearson" and "kendall")
 #                 can be used (please look at: R>?cor.test)
 # 12. widthwin_1: first value of the size of the windows (by default is 3) 
 # 	 	   when the option typewidthwin="PARTIAL" is selected. 
 # 13. widthwin_N: last value of the size of the windows (by default is 
 #                 dim(inputdata)[1], i.e. (number of datums in inputdata) 
 # 		   when the option typewidthwin="PARTIAL" is selected 
 # 14. typewidthwin: "FULL" is to estimate and to plot windows from 3, 5,..., 
 #                 to dim(inputdata)[1]), the other option is "PARTIAL" 
 # 15. Align:      to align the rolling object, RolWinMulCor ONLY uses 
 #                 the "center" option (please look at: R>?running) 
 # 16. pvalcorectmethod: p-value correction method (please look at: 
 #                 R>?p.adjust), by default we use the method of Benjamini 
 # 		   & Hochberg (BH), but other 6 methods can be used
 # 17. rmltrd:     remove (by default is "Y" or "y"; please use "N" or 
 #                 "n" otherwise) linear trends in the data analysed 
 # 18. Scale:      scale (by default is "Y" or "y"; please use "N" or "n" 
 #   	           otherwise) to "normalize" or "standardize" the data analysed
 # 19. device:     the kind of output plot (please look at: R?device), 
 #                 you can use: "png", "jpeg/jpg", "eps", "pdf" and "screen" 
 # 20. Hfig:       plot's height (device) in "png" and "jpg" format (>R?pdf)
 # 21. Wfig:       plot's width (device) in "png" and "jpg" format (>R?postscript)
 # 22. resfig:     resolution for the plot in "png" and "jpg" format (>R?jpeg)
 # 23. Hpdf:       plot's height (device) in "pdf" or "eps" format (>R?pdf)
 # 24. Wpdf:       plot's width (device) in "pdf" or "eps" format (>R?postscript)
 # 25. NUMLABX:    number of labels for ALL X's axis 
 # 26. ofilename:  output file name
 #------------------------------------------------------------------------#

 # To reset the par options modified within this function!
 oldpar <- par(no.readonly = TRUE)
 on.exit(par(oldpar))   

 #:: Devices options: png, jpg, eps & pdf! 
 path <- tempdir()
 if (device == "png" || device == "PNG") {
  fileout <- file.path(path, paste("heatmap_bivariate_", ofilename, 
	      ".png", sep=""))
  png(fileout, height=Hfig, width=Wfig, res=resfig) 
  par(mar=c(2.35, 4.95, 2.2, 3.5) + 0.1) 
  layout(matrix(c(1,2,3), 3, 1, byrow=FALSE), heights=c(2, 3.35, 0.85))
 }

 if (device == "jpeg" || device == "jpg" || device == "JPG") {
  fileout <- file.path(path, paste("heatmap_bivariate_", ofilename, 
      	      ".jpg", sep=""))
  jpeg(fileout, height=Hfig, width=Wfig, res=resfig)
  par(mar=c(2.35, 4.95, 2.2, 3.5) + 0.1) 
  layout(matrix(c(1,2,3), 3, 1, byrow=FALSE), heights=c(2, 3.35, 0.85))
 }

 if (device == "pdf" || device == "PDF") {
  fileout <- file.path(path, paste("heatmap_bivariate_", ofilename, 
	      ".pdf", sep=""))
  pdf(fileout, height=Hpdf, width=Wpdf)
  par(mar=c(2.35, 4.95, 2.2, 3.5) + 0.1) 
  layout(matrix(c(1,2,3), 3, 1, byrow=FALSE), heights=c(2, 3.5, 0.65))
 }

 if (device == "eps" || device == "EPS") {
  fileout <- file.path(path, paste("heatmap_bivariate_", ofilename, 
	      ".eps", sep=""))
  postscript(file=fileout, height=Hpdf, width=Wpdf)
  par(mar=c(2.35, 4.95, 2.2, 3.5) + 0.1) 
  layout(matrix(c(1,2,3), 3, 1, byrow=FALSE), heights=c(2, 3.5, 0.65))
 }

 # ----------------------------------------------------------------------- 
 # Check 1: number of columns, inputdata MUST contain three columns: 
 #          time, variable 1, and variable 2 
 if (dim(inputdata)[2] != 3) 
  stop("The input data MUST be an array or matrix with 3 columns 
   (first column the time and the second and third the variables 
   under analysis. Thank you for using our RolWinMulCor package.)")

 # Check 2: the time steps MUST be regular/evenly - no gaps! 
 Deltat <- diff(inputdata[,1])  # Deltat is the temporal resolution! 
 if (length(unique(Deltat)) != 1)
  stop("The input data have some gaps (it's irregular), please, consider 
   to face this before use RolWinCor. Thank you for using our 
   RolWinMulCor package.")

 # Check 3: removing linear trend - if rmltrd="Y" 
 if (rmltrd == "Y" || rmltrd == "y") {
  ts1.tmp   <- cbind(inputdata[ ,1], c(pracma::detrend(inputdata[ ,2])))
  ts2.tmp   <- cbind(inputdata[ ,2], c(pracma::detrend(inputdata[ ,3])))
  inputdata <- cbind(ts1.tmp[ ,1], ts1.tmp[,2], ts2.tmp[ ,2])  
 } 

 # Check 4: scaling data: [X_i - mean(X_i)]/sd(X-i), mean=0 & sd=1
 if (Scale == "Y" || Scale == "y") { 
  inputdata <- cbind(inputdata[,1], apply(inputdata[,-1], 2, scale))
 }

 if (typewidthwin == "PARTIAL") { 
  # Check 5: widthwin_1/N MUST be odd (of the form 2m + 1, where m is 1, 2,...)
  if (widthwin_1 %% 2 == 0 || widthwin_N %% 2 == 0) { 
   stop(paste("widthwin_1 or widthwin_N is/are EVEN and these (both) 
   MUST be ODD. Thank you for using our RolWinMulCor package."))
  }

  # Check 6: initial and final values for the window lengths! 
  if (widthwin_1 == widthwin_N) {
   stop("The initial and final value for the window-lengths are 
    the same. Please, modify these values. Thank you for using 
    our RolWinMulCor package.") 
   }
 }
 # ----------------------------------------------------------------------- 

 # Transforming input data to time series object 
 NL  <- dim(inputdata)[1]
 ts1 <- ts(inputdata[ ,2], start=inputdata[1,1], end=inputdata[NL,1], 
  	   deltat=unique(Deltat))
 ts2 <- ts(inputdata[ ,3], start=inputdata[1,1], end=inputdata[NL,1], 
	   deltat=unique(Deltat))

 time.runcor <- time(ts1)
 Nrun        <- length(time.runcor)   

 # ------------------------------------------------------------------
 # Procedure to estimate number of windows to compute the correlations
 # ------------------------------------------------------------------
 # typewidthwin indicates the way how the rolling window is computed:
 # "FULL" the window correlations are computed for all the window-lengths 
 # from 3 to NL (number of datums in inputdata). PARTIAL the window 
 # correlations are computed from widthwin_1 to widthwin_N. 
 # nwin is the maximun number of windows in the heatmap.
 # NL is the number of elements of the time series under study.
 # ------------------------------------------------------------------
 if (typewidthwin == "FULL") { 
  if (NL %%2 == 0) { 
   nwin <- NL/2 - 1 
  } else {
    if (NL %%2 != 0)   { 
    nwin <- floor(NL/2) 
    } 
  } 
   Rwidthwin <- seq(3, NL, 2) # length(Rwidthwin) = nwin 
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
 # ------------------------------------------------------------------

 # Array/matrix to save the cor. coef. and p-values 
 the_matrixCOR <- array(NA, c(nwin, NL-2))
 the_matrixPVA <- array(NA, c(nwin, NL-2))
 
 #############	        Auxiliary   function              ############
 #Function to estimate the correlation coefficients and p-values 
 pvalues.fun    <- function(ts1, ts2){cor.test(ts1, ts2, 
 		     method=CorMethod)$p.value}
 cor.def        <- function(ts1, ts2){cor(ts1, 
	 	     ts2, method=CorMethod)}

 #############	START:	The big-for 	#############
 for (w in 1:nwin) { 
  pvalrc.ts1_ts2 <- gtools::running(ts1, ts2, fun=pvalues.fun,
                     width=Rwidthwin[w], align=Align) 
  rc.ts1_ts2     <- gtools::running(ts1, ts2, fun=cor.def,
		     width=Rwidthwin[w], align=Align) 
  # 
  ncompa         <- length(rc.ts1_ts2) 
  CORTD_pval_ts1_ts2  <- round(p.adjust(pvalrc.ts1_ts2, 
                          method=pvalcorectmethod, n=ncompa), 4)
  left_win <- (Rwidthwin[w] - 1)/2 
  righ_win <- (Rwidthwin[w] - 1)/2 + 1
  the_matrixCOR[w,left_win:(Nrun-righ_win)] <- rc.ts1_ts2 
  the_matrixPVA[w,left_win:(Nrun-righ_win)] <- CORTD_pval_ts1_ts2
 }
 #############	END:	The big-for 	#############
 
 # To take into account the statistical significance in the image-plot! 
 id.xy <- which(the_matrixPVA > 0.05, arr.ind=TRUE) 
 for (k in 1:dim(id.xy)[1]) { 
   the_matrixCOR[id.xy[k,1], id.xy[k,2]] <- NA
 } 

 # Splitting up the plot -some pieces of code come from my W2CWM2C R pack. 
 # Palette!
 Ncol     <- length(Rwidthwin)
 # Number of colors MUST be length(Rwidthwin) or Ncol 
 Palette  <- colorspace::diverge_hcl(4*Ncol, c=c(100,0), l=c(50,90), power=1.3)
 # Colorbar! 
 rangev   <- seq(min(the_matrixCOR, na.rm=TRUE), max(the_matrixCOR, 
                 na.rm=TRUE), length.out=11)
 rangebar <- matrix(rangev, nrow=1, ncol=11, byrow=TRUE)
 
 # -------------------------------------------------------------
 # Setting up the graphical parameters
 if (device == "screen") { 
  par(mar=c(2.35, 4.95, 2.2, 3.5) + 0.1) 
  layout(matrix(c(1,2,3), 3, 1, byrow=FALSE), heights=c(2, 3.5, 0.65))
 }
 # -------------------------------------------------------------
 # Plot data (plot 1)
 plot(ts1, t="l", col=coltsX, las=1, xlab="", ylab="", xaxt="n", 
  yaxt="n", xaxs="i", yaxs="i", xlim=c(time(ts1)[1], time(ts1)[NL]), 
  main=paste(varnametsX, " vs. ", varnametsY, sep=""), lwd=LWDtsX)
 points(ts2, t="l", col=coltsY, las=1, xlab="", ylab="", lwd=LWDtsY)
 axis(1, at=floor(seq(time.runcor[1], time.runcor[Nrun], length.out=NUMLABX)),
  labels=floor(seq(time.runcor[1], time.runcor[Nrun], length.out=NUMLABX)),
  cex.axis=CEXAXIS)
 axis(2, at=pretty(ts1), labels=pretty(ts1), col=coltsX, las=1, 
  cex.axis=CEXAXIS)
 axis(4, at=pretty(ts2), labels=pretty(ts2), col=coltsY, las=1, 
  cex.axis=CEXAXIS, cex.lab=CEXLAB)
 mtext(1, text="Time", line=2.75, cex=CEXLAB)
 mtext(2, text=varnametsX, col=coltsX, line=2.5, cex=CEXLAB, las=1)
 mtext(4, text=varnametsY, col=coltsY, line=2.5, cex=CEXLAB, las=1)
 # -------------------------------------------------------------
 # Multiscale Window correlation (plot 2)
 image(t(the_matrixCOR), xlab="", ylab="", las=1, 
  col=Palette, xaxt="n", yaxt="n")
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
 mtext(2, text="Scales", line=3.35, cex=CEXLAB)
 # -------------------------------------------------------------
 # Colorbar (plot 3) 
if(1){
 image(z=t(rangebar), axes=FALSE, col=Palette, frame.plot=TRUE,
  yaxt="n", xaxt="n") 
 axis(1, at=round(seq(0,1,length.out=11),2), labels=round(rangebar, 
  digits=2), cex.axis=CEXAXIS, las=1)
}

 if (device != "screen") 
 dev.off()

 # -------------------------------------------------------------
 # Outputs
 namesLS     <- c("matcor", "pvalscor", "NoWindows", "Windows")  
 LIST        <- list(the_matrixCOR, the_matrixPVA, nwin, Rwidthwin)
 names(LIST) <- namesLS

 return(LIST)  

} # End-Main-function 
 
