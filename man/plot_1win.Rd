\name{plot_1win}
\alias{plot_1win}
\alias{plot_1win}
\title{Plot the outputs of \code{\link{rolwincor_1win}} and \code{\link{rolwinmulcor_1win}} as a single one window 
}
\description{
The \code{\link{plot_1win}} function plots the correlation coefficients and their respective p-values (corrected or not corrected) as only one selected window-length (time-scale) using the outputs of the functions \code{\link{rolwincor_1win}} (bi-variate case) and \code{\link{rolwinmulcor_1win}} (multi-variate case). The \code{\link{plot_1win}} function is highly flexible since this contains several parameters to control the plot output. We would highlight that only the first 11 parameters (and \emph{LWDtsX} and \emph{LWDtsY} in case the users need to customize the line widths of the time series under study) must be defined by the users, the other parameters are not strictly necessary since these have been defined by default. A list of parameters are described in the following lines. 
}
\usage{
plot_1win(inputdata, corcoefs, pvalues, left_win, righ_win, widthwin, 
          KCASE="", varX="", varY="", coltsX=c("black"), coltsY="blue", 
          rmltrd=TRUE, Scale=TRUE, HeigWin1=2.05, HeigWin2=2.75, 
          colCOEF="black", colPVAL="gray", CEXLAB=1.15, CEXAXIS=1.05, 
          LWDtsX=c(1), LWDtsY=1, LWDcoef=1, LWDpval=1, NUMLABX=5, 
          parcen=c(0.5,25)) 
}
\arguments{
 \item{inputdata}{
Input data used in the functions \code{\link{rolwincor_1win}} or \code{\link{rolwinmulcor_1win}}. 
} 
 \item{corcoefs, pvalues}{
Correlation coefficients obtained from the functions \code{\link{rolwincor_1win}} or \code{\link{rolwinmulcor_1win}} (named \emph{Correlation_coefficients}) and p-values obtained from the aforementioned functions (named \emph{P_values_corrected} or \emph{P_values_not_corrected}).
}
 \item{left_win, righ_win}{
These parameters are used to accommodate the times in the rolling window correlations and are obtained from the functions \code{\link{rolwincor_1win}} or \code{\link{rolwinmulcor_1win}}, which have the same names.
}
 \item{widthwin}{
Window size to compute the rolling window correlations. This value can be an even or odd number of at least three (the default value), and this parameter is the same as the one used in \code{\link{rolwincor_1win}} or \code{\link{rolwinmulcor_1win}}.
}
 \item{KCASE}{
This parameter is used to activate the cases: ``BIVAR'' for the bi-variate or ``MULVAR'' for the multi-variate, and this must be the same as the one used in \code{\link{rolwincor_1win}} or \code{\link{rolwinmulcor_1win}}. 
}
 \item{varX}{ 
Name of the  ``first'' or independent variable, e.g. ``X'' (please note that ``X'' is a vector of one element if KCASE=``BIVAR'' and a vector of several elements if KCASE=``MULVAR''. For the multi-variate case the names for ``X'' (the independent variables) will be defined as: varX=paste(``X1'', ``X2'',..., sep='', ''). 
}
 \item{varY}{
Name of the ``second'' (bi-variate case) or dependent variable ( multi-variate case), e.g. ``Y''. 
}
 \item{coltsX, coltsY}{
Colors to be used when the variables are plotted, for the bi-variate case by default are ``black'' for ``X'' and ``blue'' for ``Y'', but other colors can be used. For the multi-variate case, colors for the dependent (``Y'') and independent variables (``X'') MUST be provided (e.g. coltsX=c("red","blue",...), coltsY="black").    
}
 \item{rmltrd}{
Remove (by default is ``TRUE''; ``FALSE'' otherwise) the linear trend in the time series under analysis. 
}
 \item{Scale}{
Scale (by default is ``TRUE''; ``FALSE'' otherwise) is used to ``normalize'' or ``standardize'' the time series under analysis.
}
 \item{HeigWin1, HeigWin2}{
Proportion of window's size to plot the time series under analysis (HeigWin1) and the rolling window correlation coefficients and p-values (HeigWin2) (look at: \R>?layout to get more information about ``layout''). By default HeigWin1 and HeigWin2 have values of 2.05 and 2.75, but other values can be used. 
}
 \item{colCOEF, colPVAL}{
The colors to be used when the correlation coefficients and their corresponding p-values are plotted, by default the colors are ``black'' and ``gray,'' but other colors can be used. 
}
 \item{CEXLAB, CEXAXIS}{
These parameters are used to plot the sizes of the X-axis and Y-axis labels and X- and Y-axis, by default these parameters have values of 1.15 and 1.05, respectively, but it is possible to use other values. 
}
 \item{LWDtsX, LWDtsY}{
Line-widths for the first and the second variable when these are plotted, for the bi-variate case by default these have values of 1, but other values (widths) can be used. For the multi-variate case and for the independent variables the line-widths MUST be provided (e.g. LWDtsX = c(1,2,...)). 
}
 \item{LWDcoef, LWDpval}{
The line-widths to be used when the correlation coefficients and their respective p-values are plotted, by default these parameters have a value of 1, but it is possible to use other values.
}
 \item{NUMLABX}{
Number of labels for (all) the X's axis, by the default is 5, but it is possible to use other values. 
}
 \item{parcen}{
These parameters contain two values: the first one is to control the position of the title, by default it is 0.5, but you should try with other close values to obtain the title centered, e.g. 0.4 or 0.8 (please avoid to use large values); the second value is to define the spaces between the names of variables, by default is 25 spaces, but you could try other values to fit properly the names of variables in the title. We use ``mtext'' to produce the title (please loot at \R>?mtext for more information). 
}
}
\details{
The \code{\link{plot_1win}} function plots the correlation coefficients and their respective p-values (corrected or not corrected) as only one selected window-length using the outputs of the functions \code{\link{rolwincor_1win}} (bi-variate case) and \code{\link{rolwinmulcor_1win}} (multi-variate case).  
}
\value{
Output: a single plot (via screen) of the correlation coefficients and their respective (corrected or not corrected) p-values.  
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
DeustoTech - Deusto Institute of Technology, \cr
Faculty of Engineering, University of Deusto, \cr
Avda. Universidades, 24, Bilbao, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com}, \email{josue.polanco@deusto.es} \cr
}
\references{
Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57 (1), 289-300. <URL: \url{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1995.tb02031.x}>. \cr

Polanco-Martínez, J. M. (2019). Dynamic relationship analysis between NAFTA stock markets using nonlinear, nonparametric, non-stationary methods. Nonlinear Dynamics, 97(1), 369-389. <URL: \url{https://doi.org/10.1007/s11071-019-04974-y}>. \cr 

Polanco-Martínez, J. M. (2020). RolWinMulCor : an R package for estimating rolling window multiple correlation in ecological time series. Ecological Informatics (Ms. ECOINF-D-20-00263 accepted for publication, 19/08/2020).\cr

Telford, R.: Running correlations -- running into problems (2013). <URL: \cr
\url{https://quantpalaeo.wordpress.com/2013/01/04/running-correlations-running-into-problems/}>. \cr
}
\examples{
#######################################################################
# Testing the bi-variate case (1 window), synthetic data!
#######################################################################
test_fun1 <- rolwincor_1win(syntDATA, varX="X", varY="Y", 
                            CorMethod="spearman", widthwin=21, 
                            Align="center", pvalcorectmethod="BH")
# Plotting the bi-variate case (1 window)
plot_1win(syntDATA, test_fun1$Correlation_coefficients, test_fun1$P_values_corrected, 
         test_fun1$left_win, test_fun1$righ_win, widthwin=21, KCASE="BIVAR", 
         varX="X", varY="Y")
#######################################################################
# Testing the multi-variate case (1 window), real-life ecological data! 
#######################################################################
test_fun2 <- rolwinmulcor_1win(YX_ecological_data, widthwin=61,
                            Align="center", pvalcorectmethod="BH")
# Plotting the bi-variate case (1 window), real-life ecological data 
plot_1win(YX_ecological_data, test_fun2$Correlation_coefficients, test_fun2$P_values_corrected, 
         test_fun2$left_win, test_fun2$righ_win, widthwin=21, KCASE="MULVAR", varY="PC1", 
         varX=paste("SST", "TSI", sep=", "), coltsY="black", coltsX=c("red", "orange"), 
         CEXLAB=1.15, CEXAXIS=1.65, LWDtsX=rep(2,2), LWDtsY=2, parcen=c(0.45,15))
}
\keyword{only_one_simple_plot}
\keyword{only_one_simple_window_bi_and_multi_variate}

