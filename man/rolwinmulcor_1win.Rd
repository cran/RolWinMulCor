\name{rolwinmulcor_1win}
\alias{rolwinmulcor_1win}
\alias{rolwinmulcor_1win}
\title{Estimate and plot the Rolling Window Correlation for the multi-variate case 
}
\description{
The \code{\link{rolwinmulcor_1win}} function estimates the rolling (running) window correlation among multiple time series (multi-variate case) sampled on identical time points for ONLY ONE window-length and plots the correlation coefficients and their respective p-values. The \code{\link{rolwinmulcor_1win}} function is highly flexible since this contains a great number of parameters to control the estimation of correlation and the features of the plot output. For example, \code{\link{rolwinmulcor_1win}} function contains parameters to remove the (linear) trend contained in the time series under analysis, to choose different p-value correction methods (which are used to address the multiple comparison problem), or to personalise the plot output. A list of parameters are described in the following lines. 
}
\usage{
rolwinmulcor_1win(inputdata, varnametsY="", varnametsX="", units="", 
                  coltsY="black", coltsX="", colCOEF="black", 
                  colPVAL="gray", CEXLAB=1.15, CEXAXIS=1.05, LWDtsY=1, 
                  LWDtsX="", LWDcoef=1, LWDpval=1, CorMethod="spearman", 
                  widthwin=3, Align="center", pvalcorectmethod="BH", 
                  rmltrd="Y", Scale="Y", device="screen", Hfig=1200, 
                  Wfig=900, resfig=150, Hpdf=5.15, Wpdf=13.5, 
                  HeigWin1=2.05, HeigWin2=2.75, ofilename)  
}
\arguments{
 \item{inputdata}{
Matrix of P columns: time, dependent variable (``Y''), and independent variables (``X1'', ``X2',..., ``Xp-2'').
} 
 \item{varnametsY}{
Name of the dependent variable: ``Y''. Please note that the name of this variable MUST be defined. 
} 
 \item{varnametsX}{
Name of the independent variables: ``X1'', ``X2',... ``Xp-2''. Please note that the names of these variables MUST defined in this way: paste("X1", "X2",..., sep=", ").
}
 \item{units}{
Time's unit (e.g. days, weeks, years, etc.) for the time series under analysis. Please note that the units MUST de defined. 
}
 \item{coltsY, coltsX}{
The colors to be used when the dependent (``Y'') and independent (``X1'', ``X2',...) variables are plotted, by default for the dependent variable the color is ``black'' (although other colors can be used) and the colors for the independent variables MUST be provided (e.g. coltsX=c("red","blue",...)). 
}
 \item{colCOEF, colPVAL}{
The colors to be used when the correlation coefficients and their corresponding p-values are plotted, by default the colors are ``black'' and ``gray,'' but t is possible to use other colors. 
}
 \item{CEXLAB, CEXAXIS}{
The parameters are used to plot the sizes of the X-axis and Y-axis labels and X- and Y-axis, by default these parameters have values of 1.15 and 1.05, respectively, but it is possible to use other values. 
}
 \item{LWDtsY, LWDtsX}{
These parameters are used to define the line-widths when the dependent and independent variables are plotted, by default and for the dependent variable this has a value of 1 (although other values (widths) can be used), and for the independent variables the line-widths MUST be provided (e.g. LWDtsX = c(1,2,...)). 
}
 \item{LWDcoef, LWDpval}{
The line-widths to be used when the correlation coefficients and their respective p-values are plotted, by default these parameters have a value of 1, but other values (widths) can be used. 
}
 \item{CorMethod}{
The method used to estimate the correlations, by default is ``spearman'' but other options (``pearson'' and ``kendall'') are available (please look at: \R>?cor.test).
}
 \item{widthwin}{
Window's size to compute the rolling window correlations, this MUST be odd (i.e. of the form \bold{2m + 1}, where \bold{m} is a natural number), by default a value of 3 is used, but this information MUST be provided. 
}
 \item{Align}{
To align the rolling object, RolWinMulCor ONLY uses the ``center'' option (please look at: \R>?running) to ensure that variations in the correlation are aligned with the variations in the relationship of the time series under study rather than being shifted. 
}
 \item{pvalcorectmethod}{
The p-value correction method to be used, by default the method of Benjamini & Hochberg (BH) is used since this is less conservative and performs much better than Bonferroni, but other six methods (e.g. Holm, Bonferroni, etc.) are available (please look at: \R>?p.adjust).
}
 \item{rmltrd}{
Remove (by default is ``Y'' or ``y''; please use ``N'' or ``n'' otherwise) the linear trend in the time series under analysis. 
}
 \item{Scale}{
Scale (by default is ``Y'' or ``y''; please use ``N'' or ``n'' otherwise) is used to ``normalize'' or ``standardize'' the time series under analysis.
}
 \item{device}{
Kind of plot output (please look at: \R>?device), there are five options: ``png'', ``jpeg/jpg'', ``eps'', ``pdf'' and ``screen.'' By default device is ``screen.''
}
 \item{Hfig, Wfig}{
Plot's height and width (for the device) for the plot in ``png'' and ``jpg'' format (look at \R>?png or \R>?jpg), by default Hfig and Wfig have values of 1200 and 900, but other values can be used. 
}
  \item{resfig}{
Image resolution (in ``ppi'') for the plot in ``png'' and ``jpg'' format (look at \R>?png , \R>?jpg or \R>?jpeg), by default this has a value of 150, but other mage resolutions can be used. 
}
 \item{Hpdf, Wpdf}{
Plot's height and width (for the device) for ``pdf'' or ``eps'' format (>\R?pdf or >\R?postscript), by default Hpdf and Wpdf are equal to 5.15 and 13.5, but it is possible to use other values. 
}

 \item{HeigWin1, HeigWin2}{
Proportion of window's size to plot the time series under analysis (HeigWin1) and the rolling window correlation coefficients and p-values (HeigWin2) (look at: \R>?layout to get more information). By default HeigWin1 and HeigWin2 have value of 2.05 and 2.75, but it is possible to use other values.
}
 \item{ofilename}{
Output file name. 
}
}
\details{
The \code{\link{rolwinmulcor_1win}} function estimates the rolling window correlation among multiple time series (multi-variate case) sampled on identical time points for ONLY ONE window-length and plots the rolling correlation coefficients and their respective p-values. \code{\link{rolwinmulcor_1win}} uses the functions \emph{running} (package:gtools), the function \emph{cor.def} that we have created for our package \R \emph{RolWinMulCor}, and the native \R function \emph{p.adjust} (package:stats). 
}
\value{
Outputs:

Plot output: \emph{screen} or \emph{'plot_multivariate_' + 'ofilename + .png, .jpg, .eps} or \emph{.pdf'}. 

Numerical output: a list containing \emph{Correlation_coefficients}, \emph{P_values_corrected}, and \cr 
\emph{P_values_not_corrected}, which are self-explanatory. 
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
DeustoTech - Deusto Institute of Technology, \cr
Faculty of Engineering, University of Deusto, \cr
Avda. Universidades, 24, Bilbao, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{http://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com}, \email{josue.polanco@deusto.es} \cr
}
\references{
Abdi H. Multiple correlation coefficient, in Encyclopedia of Measurement and
Statistics, N. J. Salkind, Ed. Sage, Thousand Oaks, CA, USA, 2007; 648-651. \cr
<URL: \url{https://personal.utdallas.edu/~herve/Abdi-MCC2007-pretty.pdf}>. \cr

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57 (1), 289-300. <URL: \url{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1995.tb02031.x}>. \cr

Polanco-Martínez, J. M. (2019). Dynamic relationship analysis between NAFTA stock markets using nonlinear, nonparametric, non-stationary methods. Nonlinear Dynamics, 97(1), 369-389. <URL: \url{https://doi.org/10.1007/s11071-019-04974-y}>. \cr 

Telford, R.: Running correlations -- running into problems (2013). <URL: \cr
\url{https://quantpalaeo.wordpress.com/2013/01/04/running-correlations-running-into-problems/}>. \cr
}
\examples{
 # Loading packages
 library("RolWinMulCor") 
 library("gtools") 
 library("pracma") 
 # Loading data set  
 data(synthetic_data)
 # Testing the function 'rolwinmulcor_1win' (multi-variate case)
 # Window length = 51 and plot output in screen format 
 test1_rolwinmulcor_1win <- rolwinmulcor_1win(synthetic_data, varnametsY="Y", 
                             varnametsX=paste("X1", "X2", "X3", sep=", "), 
                             coltsX = c("red", "blue", "green"), units="NU", 
                             widthwin=51, device="screen", ofilename="test1")
 #  Window length = 75 and plot output in screen format 
 test2_rolwinmulcor_1win <- rolwinmulcor_1win(synthetic_data, varnametsY="Y", 
                             varnametsX=paste("X1", "X2", "X3", sep=", "), 
                             coltsX = c("red", "blue", "green"), units="NU", 
                             widthwin=75, device="screen", ofilename="test2")
}
\keyword{rolwincor_multivariate}
\keyword{rolwincor_multivariate_simply_red}
