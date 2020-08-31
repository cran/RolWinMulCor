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
#   you can redistribute it and/or modify it under the terms of the GNU 
#   General Public License as published by the Free Software 
#   Foundation, either version 3 of the License, or (at your option) 
#   any later version.
#
#   RolWinMulCor is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with RolWinMulCor If not, see <http://www.gnu.org/licenses/>.
#
#####################################################################

 rolwinmulcor_1win <- function(inputdata, varnametsY="", varnametsX="", 
	               widthwin=5, Align="center", pvalcorectmethod="BH", 
	               rmltrd="TRUE", Scale="TRUE") { 

 #---------------------------------------------------------------------#
 # rolwinmulcor_1win function estimates the rolling (running) window 
 # correlation (coefficients and p-values) among multiple time series 
 # (multi-variate case) sampled on identical time points for ONLY ONE 
 # window-length.
 #---------------------------------------------------------------------#
 # Description of variables (INPUTS): 
 # 1. inputdata:  matrix of P columns (time, dependent variable ("Y"), 
 #                independent variables ("X1", "X2",... "Xp-1"). 
#####################################################################
 # 2. varnametsY: name of the dependent variable: "Y", please note that 
 #                you MUST be define the name of this variable.
 # 3. varnametsX: name of the independent variables "X1", "X2", ..., please
 #  		  note that you MUST define the names of these variables, 
 #                please use varnametsX=paste("X1", "X2",..., sep=", "). 
 # 4. widthwin:  window's size to compute the rolling window correlations, 
 #   		 this MUST be >= 5. 
 # 5. Align:     to align the rolling object, RolWinMulCor uses 
 #               the options: "left", "center" and "right." Please 
 #               note that if "widthwin" is even it's not possible to use 
 #   	    	 the Align option "center" (please look at: R>?running).
 # 6. pvalcorectmethod: p value correction method (please look at: 
 #               R>?p.adjust), by default we use the method of Benjamini 
 # 		 & Hochberg (BH) (1995), but other six methods can be used. 
 # 7. rmltrd:    remove (by default is TRUE; please use "FALSE" 
 # 	         otherwise) linear trend in the data analysed. 
 # 8. Scale:     scale (by default is "TRUE"; please use "FALSE" otherwise)
 # 		 or "normalized" or "standardized" the data analysed. 
 #------------------------------------------------------------------------#
 
 # ----------------------------------------------------------------------- 
  # Check 1: number of columns, inputdata MUST contain at least four
  #          columns: time, dependent variable Y, and two (or more)
  # 	     independent variables X1, X2,...  
  if (dim(inputdata)[2] < 4) 
   stop("\n The input data MUST be an array/matrix with at least 4 columns 
    (first column the time, second the dependent variable (Y) and 
    the others the independent variables (X1, X2,...). Thank you for 
    using our RolWinMulCor package. \n")

  # Check 2: the time steps MUST be regular/evenly - no gaps! 
  Deltat <- diff(inputdata[,1])  # Deltat is the temporal resolution! 
  if (length(unique(Deltat)) != 1)
   stop("The input data (time) have some gaps (it's irregular), please, 
    consider to face this drawback before using RolWinMulCor. We 
    recommend you our BINCOR package and method (also in CRAN; 
    https://cran.r-project.org/package=BINCOR), but other packages 
    and methods can be used. Thank you for using our RolWinCor 
    package. \n")

############################################################################
 # Check 3: widthwin MUST be odd if Align="center" otherwise 
 # gtools::running will not work! 
 if (widthwin %% 2 == 0 & Align == "center" | widthwin < 5) { 
  stop("\n 'widthwin' is EVEN and 'Align' has been defined as 
   `center' or 'widthwin' is < 5.  Thank you for using our 
   RolWinMulCor package. \n")
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

  # Transforming input data to time series object! 
  NL  <- dim(inputdata)[1]
  NP  <- dim(inputdata)[2]
  ts1 <- ts(inputdata[,2], start=inputdata[1,1], end=inputdata[NL,1], 
     	    deltat=unique(Deltat)) # ts1 = Y
  ts2 <- ts(inputdata[,3:NP], start=inputdata[1,1], end=inputdata[NL,1], 
	    deltat=unique(Deltat)) #ts2 = X (X1, X2, ..., Xp) 
  Z   <- cbind(ts1, ts2) # To be used in the running cor./fun. (rollapply)

  time.runcor <- time(ts1)
  Nrun        <- length(time.runcor) 

  # ------------------------------------------------------------------
  # Computing the rolling window (running) correlation and p-values  

  rc.ts1_ts2_tmp <- t( zoo::rollapply(Z, width=widthwin, 
    #############        Auxiliary   function             ############
    FUN=function(Z) { 
    P <- dim(Z)[2]
    # where Z[,1] is the dependent and Z[,2:P] the independent variable/s
    lm_estimate <- lm(Z[,1] ~ Z[,2:P], data=as.data.frame(Z))
    summary_lm  <- summary(lm_estimate) 
    adjRsqu     <- summary_lm$adj.r.squared
    Fstat       <- summary_lm$fstatistic[1]
    pvalueF     <- pf(summary_lm$fstatistic[1], summary_lm$fstatistic[2], 
                   summary_lm$fstatistic[3], lower.tail=F)
    namesLS     <- c("Squared_multiple_cor_coef", "F-statistics", "P-value")
    LIST        <- cbind(adjRsqu, Fstat, pvalueF) 
    return(LIST)
    }, 
   by.column=FALSE, align=Align) ) 

  rc.ts1_ts2     <- rc.ts1_ts2_tmp[1,]
  pvalrc.ts1_ts2 <- rc.ts1_ts2_tmp[3,]

  # p-value correction (p.adjust)
  ncompa 	      <- length(rc.ts1_ts2) 
  CORTD_pval_ts1_ts2  <- round(p.adjust(pvalrc.ts1_ts2, 
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

 # Numerical output 
 namesLS <- c("Correlation_coefficients", "P_values_corrected", 
              "P_values_not_corrected", "left_win", "righ_win", "widthwin") 
 LIST    <- list(cbind(time.runcor[left_win:(Nrun-righ_win)], 
             rc.ts1_ts2), cbind(time.runcor[left_win:(Nrun-righ_win)], 
             pvalrc.ts1_ts2), cbind(time.runcor[left_win:(Nrun-righ_win)], 
             CORTD_pval_ts1_ts2), left_win, righ_win, widthwin)
 names(LIST) <- namesLS

 return(LIST)
 
} # End Main function 


