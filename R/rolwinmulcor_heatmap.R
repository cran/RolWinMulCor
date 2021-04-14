######################################################################
#:: rolwinmulcor_heatmap function - R package RolWinMulCor           #
#:: Multi-variate case!				 		     # 
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
#   along with RolWinMurCor If not, see <http://www.gnu.org/licenses/>.
#
#####################################################################

 rolwinmulcor_heatmap <- function(inputdata, varnametsY="", 
		varnametsX="", typewidthwin="FULL", widthwin_1=7, 
	        widthwin_N=dim(inputdata)[1], Align="center", 
	        pvalcorectmethod="BH", rmltrd=TRUE, Scale=TRUE) {

 #------------------------------------------------------------------------#
 # rolwinmulcor_heatmap function estimates the rolling (running) window 
 # correlation (coefficients and their respective p-values) among several 
 # time series (multi-variate case) sampled on identical time points for 
 # all the possible window-lengths (i.e. from three to the number of 
 # elements in the time series under analysis) or for a band of 
 # window-lengths. 
 #------------------------------------------------------------------------#
 # Description of variables (INPUTS): 
 # 1. inputdata:  matrix of P columns (time, dependent variable ("Y"), and
 #     		  independent variables ("X1", "X2", ... "Xp-1"). The 
 #		  number of independent variables MUST be at least 2. 
 # 2. varnametsY: name of the dependent variable ("Y"), please note that 
 #                you MUST be define the name of this variable. 
 # 3. varnametsX: name of the independent variables "X1", "X2", ... please
 #  		  note that you MUST define the names of these variables. 
 # 4. typewidthwin: "FULL" is to estimate and plot windows from 2, 4, ..., 
 #                to dim(inputdata)[1]) if Align is equal to "left" or 
 #                "right", or from 3, 5,..., to dim(inputdata)[1]) if  
 #                Align is "center". The other option is "PARTIAL", please 
 # 		  you should take into account that widthwin_1 and 
 #                widthwin_1 MUST be ODD if the option Align is "center".
 # 5. widthwin_1: first value for the size of the windows (by default is 5) 
 # 	 	  when the option typewidthwin="PARTIAL" is selected. 
 # 6. widthwin_N: last value for the size of the windows (by default is 
 #                dim(inputdata)[1], i.e. (number of datums in inputdata) 
 # 		  when the option typewidthwin="PARTIAL" is selected.
 # 7. Align:      to align the rolling object, RolWinMulCor uses three 	
 #                options: "left, "center", and "right" (please look at: 
 # 		  R>?running). However, there are some restrictions, 
 #                which have been described lines above. We recommend
 # 	          to use the "center" option to ensure that variations in 
 #                the correlations are aligned with the variations in the 
 # 	          relationships of the variables under study, rather than 
 #                being shifted to left or right (Polanco-Martínez 2019), 
 # 		  but this imply that the window-lengths MUST be ODD. 
 # 8. pvalcorectmethod: p-value correction method (please look at: 
 #                R>?p.adjust), by default we use the method of Benjamini 
 # 		  & Hochberg (BH) (1995), but other six methods can be used. 
 # 9. rmltrd:     remove (by default is "TRUE"; please use "FALSE" 
 #                otherwise) linear trends in data analysed. 
 # 10. Scale:     scale (by default is "TRUE"; please use "FALSE"  
 #   	          otherwise) or "normalize" or "standardize" the 
 # 		  data analysed. 
 #------------------------------------------------------------------------#

 # ----------------------------------------------------------------------- 
 # Check 1: number of columns, inputdata MUST contain at least four
 #          columns: time, dependent variable Y, and two (or more)
 # 	    independent variables X1, X2,...  
 if (dim(inputdata)[2] < 4) 
  stop("\n The input data MUST be an array or matrix with at least 4 
   columns (first column the time, second column the dependent variable 
   (Y) and the others columns the independent variables (X1, X2,...). 
   Thank you for using our RolWinMulCor package. \n")

 # Check 2: the time steps MUST be regular - no gaps! 
 #Deltat <- diff(inputdata[,1])  # Deltat is the temporal resolution! 
 #if (length(unique(Deltat)) != 1)
  cat("\n W A R N I N G: The input data must be regular (evenly spaced 
   time). Otherwise, please, consider to address this drawback before 
   using NonParRolCor. We recommend you our BINCOR package and method 
   (also in CRAN; https://cran.r-project.org/package=BINCOR), but other 
   packages and methods can be used. Thank you so much for using our 
   NonParRolCor package. \n")

############################################################################
 # Check 3.1: if Align="center" only estimate windows of the form 2p + 1
 # That is, widthwin_1/N MUST be odd
 if (typewidthwin == "PARTIAL" & Align == "center") { 
  if (widthwin_1 %% 2 == 0 || widthwin_N %% 2 == 0) { 
   stop("\n widthwin_1 or widthwin_N is/are EVEN and these (both) 
   MUST be ODD. Thank you for using our RolWinMulCor package. \n")
  }

 # Check 3.2: initial and final values for the window lengths! 
 if (widthwin_1 == widthwin_N) {
   stop("\n The initial and final value for the window-lengths are 
    the same. Please, modify these values. Thank you for using 
    our RolWinMulCor package. \n") 
  }
 }
############################################################################

 # Check 4: removing linear trend - if rmltrd="Y" 
 if(isTRUE(rmltrd)) {
  dtrd_tmp  <- pracma::detrend(inputdata[,-1])
  inputdata <- cbind(inputdata[,1], dtrd_tmp) 
  } 

 # Check 5: scaling data: [X_i - mean(X_i)]/sd(X-i), mean=0 & sd=1
 if(isTRUE(Scale)) { 
  inputdata <- cbind(inputdata[,1], apply(inputdata[,-1], 2, scale))
 }

 # ----------------------------------------------------------------------- 

 # Transforming input data to time series object 
 NL   <- dim(inputdata)[1]
 NP   <- dim(inputdata)[2]
 freq <- length((inputdata[,1]))/length(unique(inputdata[,1]))
 ts1  <- ts(inputdata[,2], start=c(inputdata[1,1],1), 
            end=c(inputdata[NL,1],freq), deltat=1/freq) # ts1 = Y
 ts2  <- ts(inputdata[,3:NP], start=c(inputdata[1,1],1), 
	    end=c(inputdata[NL,1],freq), deltat=1/freq) #ts2 = X (X1, X2,..., Xp) 
 Z   <- cbind(ts1, ts2) # To be used in the running cor./fun. (rollapply)

 time.runcor <- time(ts1)
 Nrun        <- length(time.runcor)   

 # ------------------------------------------------------------------
 # Procedure to estimate the widows and the number of windows to 
 # compute the rolling window correlations
 # ------------------------------------------------------------------
 # typewidthwin indicates how the rolling window are computed:
 # "FULL" the window correlations are computed for all the window-lengths 
 # from 2 or 3 to NL (number of datums in inputdata). PARTIAL the window 
 # correlations are computed from widthwin_1 to widthwin_N. 
 # nwin is the maximun number of windows in the heatmap. 
 # NL is the number of elements of the time series under study.
 # ------------------------------------------------------------------
 #
 if (Align == "left" || Align == "right")  { 
  if (typewidthwin == "FULL") { 
   if (NL %%2 == 0) { 
    nwin <- NL/2 - 1 
   } else {
    if (NL %%2 != 0)   { 
     nwin <- floor(NL/2) 
    } 
   } 
   Rwidthwin <- seq(4, NL, 2) # length(Rwidthwin) = nwin 
   #Rwidthwin <- seq(8, NL, 2) # length(Rwidthwin) = nwin 
  }
  if (typewidthwin == "PARTIAL") { 
   nwin      <- length(seq(widthwin_1, widthwin_N, 2)) 
   Rwidthwin <- seq(widthwin_1, widthwin_N, 2) # length(Rwidthwin) = nwin 
  }  
 }

#
 if(Align == "center") { 
  if (typewidthwin == "FULL") { 
   if (NL %%2 == 0) { 
    nwin <- NL/2 - 1 
   } else {
    if (NL %%2 != 0)   { 
     nwin <- floor(NL/2) 
    } 
   } 
   Rwidthwin <- seq(3, NL, 2) # length(Rwidthwin) = nwin 
   #Rwidthwin <- seq(7, NL, 2) # length(Rwidthwin) = nwin 
  }
  if (typewidthwin == "PARTIAL") { 
   nwin      <- length(seq(widthwin_1, widthwin_N, 2)) 
   Rwidthwin <- seq(widthwin_1, widthwin_N, 2) # length(Rwidthwin) = nwin 
  } 
 }
 
 # ------------------------------------------------------------------
 # Computing the rolling window (running) correlation and p-values  

 # Array/matrix to save the cor. coef. and p-values 
 the_matrixCOR 	   <- array(NA, c(nwin, NL-2))
 the_matrixPVA     <- array(NA, c(nwin, NL-2))
 the_matrixPVA_NOC <- array(NA, c(nwin, NL-2)) # p-values not corrected
 
 #############	BEGIN:	The big-for 	#############
 #if (typewidthwin == "FULL") { 
 for (w in 1:nwin) { 
  rc.ts1_ts2_tmp <- t( zoo::rollapply(Z, width=Rwidthwin[w], 
    #############        Auxiliary   function             ############
    FUN=function(Z) { 
    P <- dim(Z)[2]
    # where Z[,1] is the dependent and Z[,2:P] the independent variable/s
    lm_estimate <- lm(Z[,1] ~ Z[,2:P], data=as.data.frame(Z))
    summary_lm  <- summary(lm_estimate) 
    adjRsqu     <- summary_lm$adj.r.squared
     #cat("w, adjRsqu:", w, adjRsqu) 
     #if (adjRsqu < -1.0) adjRsqu <- NA
    Fstat       <- summary_lm$fstatistic[1]
    pvalueF     <- pf(summary_lm$fstatistic[1], summary_lm$fstatistic[2], 
                   summary_lm$fstatistic[3], lower.tail=F)
    namesLS     <- c("Squared_multiple_cor_coef", "F-statistics", "P-value")
     if (Rwidthwin[w] < 7) {LIST <- cbind(NA, NA, NA)} 
      else LIST        <- cbind(adjRsqu, Fstat, pvalueF) 
    return(LIST)
    }, 
   by.column=FALSE, align=Align) ) 

  rc.ts1_ts2     <- rc.ts1_ts2_tmp[1,]
  pvalrc.ts1_ts2 <- rc.ts1_ts2_tmp[3,]
  ncompa         <- length(rc.ts1_ts2) 

  CORTD_pval_ts1_ts2  <- round(p.adjust(pvalrc.ts1_ts2, 
                          method=pvalcorectmethod, n=ncompa), 4)
  ###################################################################
    # if widthwin is even and left or right  
  if(Align == "left" || Align == "right") { 
   left_win <- Rwidthwin[w]/2 
   righ_win <- Rwidthwin[w]/2 
  }
  if(Align == "center") { 
   left_win <- (Rwidthwin[w] - 1)/2   
   righ_win <- (Rwidthwin[w] - 1)/2 + 1 
  }
  ###################################################################
  #left_win <- (Rwidthwin[w] - 1)/2 
  #righ_win <- (Rwidthwin[w] - 1)/2 + 1
  the_matrixCOR[w,left_win:(Nrun-righ_win)] <- unlist(rc.ts1_ts2)
  the_matrixPVA[w,left_win:(Nrun-righ_win)] <- CORTD_pval_ts1_ts2
  the_matrixPVA_NOC[w,left_win:(Nrun-righ_win)] <- pvalrc.ts1_ts2 
 }
 #############	END:	The big-for 	#############
 
 # To take into account the statistical significance in the image-plot! 
# id.xy <- which(the_matrixPVA > 0.05, arr.ind=TRUE) 
# for (k in 1:dim(id.xy)[1]) { 
#   the_matrixCOR[id.xy[k,1], id.xy[k,2]] <- NA
# } 

 # -------------------------------------------------------------
 # Outputs
 namesLS     <- c("matcor", "pvalscor", "pvalNOTcor", "NoWindows", 
	          "Windows", "left_win", "righ_win")  
 LIST        <- list(the_matrixCOR, the_matrixPVA, the_matrixPVA_NOC,
                 nwin, Rwidthwin, left_win, righ_win)
 names(LIST) <- namesLS

 return(LIST)  

} # End-Main-function 
 
