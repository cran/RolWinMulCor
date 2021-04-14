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

 rolwincor_heatmap <- function(inputdata, varX="", varY="", 
 	        CorMethod="pearson", typewidthwin="FULL", widthwin_1=3, 
		widthwin_N=dim(inputdata)[1], Align="center", 
		pvalcorectmethod="BH", rmltrd=TRUE, Scale=TRUE) {

 #---------------------------------------------------------------------#
 # The rolwincor_heatmap function estimates the rolling (running)
 # window correlation coefficient and their respective p-values between 
 # TWO time series (bi-variate case) sampled on identical time points 
 # for all the possible (i.e. from 3 to the number of elements of 
 # the time series under analysis) window-lengths or for a band of 
 # window-lengths (from widthwin_1 to widthwin_N). 
 #---------------------------------------------------------------------#
 # Description of variables (INPUTS): 
 # 1. inputdata:  matrix of 3 columns: time, variable 1 (e.g. "X")
 #                and variable 2 (e.g. "Y").
 # 2. varX:       name of variable 1 (e.g. "X").
 # 3. varY:       name of variable 2 (e.g. "Y").
 # 4. CorMethod:  method used to estimate the correlations, by default is 
 #                "pearson" but other two options ("spearman" and "kendall")
 #                can be used (please look at: R>?cor.test).
 # 5. typewidthwin: "FULL" is to estimate the windows from 2, 4, ..., 
 #                to dim(inputdata)[1]) if Align is equal to "left" or 
 #                "right", or from 3, 5,..., to dim(inputdata)[1]) if  
 #                Align is "center". The other option is "PARTIAL", please 
 # 		  you should take into account that widthwin_1 and 
 #                widthwin_1 MUST be ODD if the option Align is "center".
 # 6. widthwin_1: first value for the size of the windows (by default is 3) 
 # 	 	  when the option typewidthwin="PARTIAL" is selected. 
 # 7. widthwin_N: Last value for the size of the windows (by default is 
 #                dim(inputdata)[1], i.e. (number of datums in inputdata) 
 # 		  when the option typewidthwin="PARTIAL" is selected.
 # 8. Align:      to align the rolling object, RolWinMulCor uses three 	
 #                options: "left, "center", and "right" (please look at: 
 # 		  R>?running). However, there are some restrictions, 
 #                which have been described lines above. We recommend
 # 	          to use the "center" option to ensure that variations in 
 #                the correlations are aligned with the variations in the 
 # 	          relationships of the variables under study, rather than 
 #                being shifted to left or right (Polanco-Martínez 2019), 
 # 		  but this imply that the window-lengths MUST be ODD. 
 # 9. pvalcorectmethod: p-value correction method (please look at: 
 #                R>?p.adjust), by default we use the method of Benjamini 
 # 		  & Hochberg (BH) (1995), but other six methods can be used. 
 # 10. rmltrd:    remove (by default is TRUE; please use "FALSE" 
 #                otherwise) linear trends in the data analysed. 
 # 11. Scale:     scale (by default is "TRUE"; please use "FALSE" 
 #   	          otherwise) or "normalize" or "standardize" the data 
 #		  analysed. 
 #------------------------------------------------------------------------#

 # ----------------------------------------------------------------------- 
 # Check 1: number of columns, inputdata MUST contain three columns: 
 #          time, variable 1 (e.g. "X"), and variable 2 (e.g. "Y")  
 if (dim(inputdata)[2] != 3) 
  stop("\n The input data MUST be an array or matrix with 3 columns 
   (first column the time and the second and third columns the 
   variables under analysis. Thank you so much for using our 
   RolWinMulCor package. \n")

 # Check 2: the time steps MUST be regular/evenly - no gaps! 
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

 # Check 4: removing linear trend - if rmltrd=TRUE
 if(isTRUE(rmltrd)) {
  ts1.tmp   <- cbind(inputdata[ ,1], c(pracma::detrend(inputdata[ ,2])))
  ts2.tmp   <- cbind(inputdata[ ,2], c(pracma::detrend(inputdata[ ,3])))
  inputdata <- cbind(ts1.tmp[ ,1], ts1.tmp[,2], ts2.tmp[ ,2])  
 } 

 # Check 5: scaling data: [X_i - mean(X_i)]/sd(X-i), mean=0 & sd=1
  if(isTRUE(Scale)) { 
  inputdata <- cbind(inputdata[,1], apply(inputdata[,-1], 2, scale))
 }

 # ----------------------------------------------------------------------- 

 # Transforming input data to time series object 
 NL   <- dim(inputdata)[1]
 freq <- length((inputdata[,1]))/length(unique(inputdata[,1]))
 ts1  <- ts(inputdata[,2], start=c(inputdata[1,1],1), 
            end=c(inputdata[NL,1],freq), deltat=1/freq)
 ts2  <- ts(inputdata[,3], start=c(inputdata[1,1],1), 
	    end=c(inputdata[NL,1],freq), deltat=1/freq)

 time.runcor <- time(ts1)
 Nrun        <- length(time.runcor)   

 # ------------------------------------------------------------------
 # Procedure to estimate the windows and the number of windows to 
 # compute the rolling window correlations
 # ------------------------------------------------------------------
 # typewidthwin indicates how the rolling windows are computed:
 # "FULL" the window correlations are computed for all the window-lengths 
 # from 2 or 3 to NL (number of datums in inputdata). PARTIAL the window 
 # correlations are computed from widthwin_1 to widthwin_N. 
 # nwin is the maximum number of windows.
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
  }
  if (typewidthwin == "PARTIAL") { 
  nwin      <- length(seq(widthwin_1, widthwin_N, 2)) 
  Rwidthwin <- seq(widthwin_1, widthwin_N, 2) # 
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
  }
 
 if (typewidthwin == "PARTIAL") { 
  nwin      <- length(seq(widthwin_1, widthwin_N, 2)) 
  Rwidthwin <- seq(widthwin_1, widthwin_N, 2) # length(Rwidthwin) = nwin 
 } 
}
 # ------------------------------------------------------------------

 # Array/matrix to save the cor. coef. and p-values 
 the_matrixCOR     <- array(NA, c(nwin, NL-2)) # correlation coef. 
 the_matrixPVA     <- array(NA, c(nwin, NL-2)) # p-values corrected 
 the_matrixPVA_NOC <- array(NA, c(nwin, NL-2)) # p-values not corrected
 
 #############	        Auxiliary   function              ############
 #Function to estimate the correlation coefficients and p-values 
 cor_pval.fun <- function(ts1, ts2){
  res_estim   <- cor.test(ts1, ts2, method=CorMethod)
  c(correlation=res_estim$estimate, p.value=res_estim$p.value)
  }

 #############	START:	The big-for 	#############
 for (w in 1:nwin) { 
  cor_pval_run <- gtools::running(ts1, ts2, fun=cor_pval.fun,
                   width=Rwidthwin[w], align=Align)

  ncompa              <- length(cor_pval_run[1,]) 
  CORTD_pval_ts1_ts2  <- round(p.adjust(cor_pval_run[2,], 
                          method=pvalcorectmethod, n=ncompa), 4)
  #cat("ncompa", ncompa)
  ###################################################################
   # if widthwin is even and left or right  
  if(Align == "left" || Align == "right") { 
   left_win <- Rwidthwin[w]/2 
   righ_win <- Rwidthwin[w]/2 
  }
  #if(Align == "right") { 
  # left_win <- Rwidthwin[w]/2 
  # righ_win <- Rwidthwin[w]/2 
  #}
  # if widthwin is odd and center 
  if(Align == "center") { 
   left_win <- (Rwidthwin[w] - 1)/2   
   righ_win <- (Rwidthwin[w] - 1)/2 + 1 
  }
  ###################################################################
  #left_win <- (Rwidthwin[w] - 1)/2 
  #righ_win <- (Rwidthwin[w] - 1)/2 + 1
  the_matrixCOR[w,left_win:(Nrun-righ_win)]     <- cor_pval_run[1,] 
  the_matrixPVA[w,left_win:(Nrun-righ_win)]     <- CORTD_pval_ts1_ts2
  the_matrixPVA_NOC[w,left_win:(Nrun-righ_win)] <- cor_pval_run[2,] 
 }
 #############	END:	The big-for 	#############
 
 # To take into account the statistical significance in the image-plot! 
 #id.xy <- which(the_matrixPVA > 0.05, arr.ind=TRUE) 
 #for (k in 1:dim(id.xy)[1]) { 
 #  the_matrixCOR[id.xy[k,1], id.xy[k,2]] <- NA
 #} 

 # -------------------------------------------------------------
 # Outputs
 namesLS     <- c("matcor", "pvalscor", "pvalNOTcor", "NoWindows", 
	          "Windows", "CorMethod", "left_win", "righ_win")  
 LIST        <- list(the_matrixCOR, the_matrixPVA, the_matrixPVA_NOC, 
	         nwin, Rwidthwin, CorMethod, left_win, righ_win)
 names(LIST) <- namesLS

 return(LIST)  

} # End-Main-function 
 
