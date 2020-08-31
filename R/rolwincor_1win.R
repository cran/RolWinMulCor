######################################################################
#:: rolwincor_1win function - R package RolWinMulCor                 #
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

 rolwincor_1win <- function(inputdata, varX="", varY="", 
 		    CorMethod="pearson", widthwin=3, Align="center", 
	            pvalcorectmethod="BH", rmltrd=TRUE, Scale=TRUE)  { 

 #---------------------------------------------------------------------#
 # "rolwincor_1win" function estimates the rolling (running) window 
 # correlation (coefficients and their respective p-values) between TWO 
 # time series (bi-variate case) sampled on identical time points for 
 # ONLY ONE window-length.
 #---------------------------------------------------------------------#
 # Description of variables (INPUTS): 
 # 1.  inputdata:  matrix of 3 columns: time, variable 1 (e.g. "X"), 
 # 		   variable 2 (e.g. "Y").  
 # 2.  varX: 	   name of variable 1 (e.g. "X"). 
 # 3.  varY: 	   name of variable 2 (e.g. "Y").
 # 4.  CorMethod:  method used to estimate the correlations, by default 
 #                 is "pearson" but other options ("spearman" and "kendall")
 #                 are possible (please see at: R>?cor.test). 
 # 5. widthwin:    window's size to compute the rolling window correlations, 
 #   		   this MUST be >= 3. 
 # 6. Align:       to align the rolling object, RolWinMulCor uses 
 #                 the options: "left", "center" and "right." Please 
 #                 note that if "widthwin" is even it's not possible to 
 #   	    	   to use the option "center" (please look at: R>?running).
 # 7. pvalcorectmethod: p-value correction method (please look at: 
 #                 R>?p.adjust), by default we use the method of Benjamini 
 # 		   & Hochberg (BH) (1995), but other 6 methods can be used
 # 8. rmltrd:      remove (by default is "TRUE"; please use "FALSE"
 #                 otherwise) linear trend in the data analysed.   
 # 9. Scale:       scale (by default is "TRUE"; please use "FALSE" 
 # 		   otherwise): "normalized" or "standardized" the data 
 # 		   analysed. 
 #------------------------------------------------------------------------#

 #------------------------------------------------------------------------#
  # Check 1: number of columns  - MUST be three columns: 
  #          time, variable X and variable Y 
  if (dim(inputdata)[2] != 3) 
   stop("\n The input data MUST be an array/matrix with 3 columns (first  
    column the time and the second and third the variables under analysis.  
    Thank you for using our RolWinMulCor package. \n ")

  # Check 2: the times should be regular/evenly spaced - no gaps! 
  Deltat <- diff(inputdata[,1])  # Deltat is the temporal resolution! 
  if (length(unique(Deltat)) != 1)
   stop("\n The input data (time) have some gaps (it's irregular), 
    please consider to face this before using RolWinMulCor.  We 
    recommend you our BINCOR package and method (also in CRAN; 
    https://cran.r-project.org/package=BINCOR), but other packages 
    and methods can be used. Thank you for using our RolWinCor 
    package. \n")

############################################################################
  # Check 3: widthwin MUST be odd if Align="center" otherwise 
  # gtools::running will not work! 
  if (widthwin %% 2 == 0 & Align == "center" | widthwin < 3) { 
   stop(paste("\n 'widthwin' is EVEN and 'Align' has been defined as 
    `center' or 'widthwin' is < 3.  Thank you for using RolWinMulCor 
     package. \n"))
  }
############################################################################

  # Check 4: removing linear trend - if rmltrd=TRUE
  if(isTRUE(rmltrd)) {
   ts1.tmp   <- cbind(inputdata[,1], c(pracma::detrend(inputdata[,2])))
   ts2.tmp   <- cbind(inputdata[,2], c(pracma::detrend(inputdata[,3])))
   inputdata <- cbind(ts1.tmp[,1], ts1.tmp[,2], ts2.tmp[,2])  
  } 
 
  # Check 5: scaling data: [X_i - mean(X_i)]/sd(X-i), mean=0 & sd=1
  if(isTRUE(Scale)) { 
   inputdata <- cbind(inputdata[,1], apply(inputdata[,-1], 2, scale))
  }
  # ----------------------------------------------------------------------- 

  # Transforming input data to time series objects! 
  NL  <- dim(inputdata)[1]
  ts1 <- ts(inputdata[,2], start=inputdata[1,1], end=inputdata[NL,1], 
     	    deltat=unique(Deltat))
  ts2 <- ts(inputdata[,3], start=inputdata[1,1], end=inputdata[NL,1], 
	    deltat=unique(Deltat))

  # ------------------------------------------------------------------
  # Computing the rolling window (running) correlation and p-values  

  cor_pval.fun <- function(ts1, ts2){
   res_estim   <- cor.test(ts1, ts2, method=CorMethod)
   c(correlation=res_estim$estimate, p.value=res_estim$p.value)
  }
 
  cor_pval_run <- gtools::running(ts1, ts2, fun=cor_pval.fun,
                    width=widthwin, align=Align)

  time.runcor  <- time(ts1)
  Nrun         <- length(time.runcor) 
  
  # p-value correction (p.adjust)
  ncompa 	      <- length(cor_pval_run[1,]) 
  CORTD_pval_ts1_ts2  <- round(p.adjust(cor_pval_run[2,], 
                          method=pvalcorectmethod, n=ncompa), 4)

######################################################################
  # if widthwin is even  
  if(widthwin %% 2 == 0 &  Align == "left") { 
   left_win <- widthwin/2 
   righ_win <- widthwin/2 
  }
  if(widthwin %% 2 == 0 &  Align == "right") { 
   left_win <- widthwin/2 + 1  
   righ_win <- widthwin/2 - 1  
  }
  # if widthwin is odd 
  if(widthwin %% 2 != 0 &  Align == "left") { 
   left_win <- floor(widthwin/2)
   righ_win <- ceiling(widthwin/2)
  }
  if(widthwin %% 2 != 0 &  Align == "center") { 
   left_win <- (widthwin - 1)/2 + 1  
   righ_win <- (widthwin - 1)/2 
  }
  if(widthwin %% 2 != 0 &  Align == "right") { 
   left_win <- ceiling(widthwin/2) + 1 
   righ_win <- floor(widthwin/2) - 1 
  }
######################################################################

 # just for testing: 
 #cat("\n left_win", left_win, "\n")
 #cat("\n righ_win", righ_win, "\n")

 # Numerical output 
 namesLS <- c("Correlation_coefficients", "P_values_corrected", 
              "P_values_not_corrected", "CorMethod", "left_win", 
	      "righ_win", "widthwin") 
 LIST    <- list(cbind(time.runcor[left_win:(Nrun-righ_win)], 
             cor_pval_run[1,]), cbind(time.runcor[left_win:(Nrun-righ_win)], 
             cor_pval_run[2,]), cbind(time.runcor[left_win:(Nrun-righ_win)], 
             CORTD_pval_ts1_ts2), CorMethod, left_win, righ_win, widthwin)
 names(LIST) <- namesLS

 return(LIST)
 
} # End function 


