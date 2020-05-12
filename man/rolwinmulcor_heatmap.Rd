\name{rolwinmulcor_heatmap}
\alias{rolwinmulcor_heatmap}
\alias{rolwinmulcor_heatmap}
\title{Estimate the Rolling Window Correlation for the multi-variate case and plot the results as a heat map. 
}
\description{
The \code{\link{rolwinmulcor_heatmap}} function estimates the rolling (running) window correlation among several time series (multi-variate case) sampled on identical time points for all the possible (i.e. from 3 to the number of elements of the time series under analysis) window-lengths or for a band of window-lengths and plots the correlation coefficients and their respective p-values (<= 0.05) as a heat map. The \code{\link{rolwinmulcor_heatmap}} function is highly flexible since this contains a great number of parameters to control the estimation of correlation and features of the plot output. A list of parameters are described in the following lines. 
}
\usage{
rolwinmulcor_heatmap(inputdata, varnametsY="Y", varnametsX="", units="", 
                     coltsY="black", coltsX="", CEXLAB=1.15, CEXAXIS=1.5, 
                     LWDtsY=1, LWDtsX="", CorMethod="spearman", typewidthwin="FULL", 
                     widthwin_1=3, widthwin_N=dim(inputdata)[1],  
                     Align="center", pvalcorectmethod="BH", rmltrd="Y", 
                     Scale="Y", device="screen", Hfig=900, Wfig=900, 
                     resfig=150, Hpdf=7, Wpdf=7, NUMLABX=5, ofilename)
}
\arguments{
 \item{inputdata}{
Matrix of P columns: time, dependent variable (``Y''), and independent variables (``X1'', ``X2',..., ``Xp-2'').
} 
 \item{varnametsY}{
Name of the dependent variable: ``Y''. Please note that the name of this variable MUST be defined. 
} 
 \item{varnametsX}{
Name of the independent variables: ``X1'', ``X2',... ``Xp-2''. Please note that the names of these variables MUST be defined in this way: paste("X1", "X2",..., sep=", ").
}
 \item{units}{
Time's unit (e.g. days, weeks, years, etc.) for the variables under analysis. Please note that the units MUST be defined. 
}
 \item{coltsY, coltsX}{
The colors to be used when the dependent (``Y'') and independent (``X1'', ``X2',...) variables are plotted, by default for the dependent variable the color is ``black'' (although other colors can be used) and the colors for the independent variables MUST be provided (e.g. coltsX=c("red","blue",...)). 
}
 \item{CEXLAB, CEXAXIS}{
These parameters are used to plot the sizes of the X-axis and Y-axis labels and X- and Y-axis, by default these parameters have values of 1.15 and 1.5, respectively, but other values can be used. 
}
 \item{LWDtsY, LWDtsX}{
These parameters are used to define the line-widths when the dependent and independent variables are plotted, by default and for the dependent variable this has a values of 1 (although other values (widths) can be used), please note that for the independent variables the line-widths MUST be provided (e.g. LWDtsX = c(1,2,...)). 
}
 \item{CorMethod}{
The method used to estimate the correlations, by default is ``spearman'' but other options (``pearson'' and ``kendall'') are available (please look at: \R>?cor.test).
}
 \item{typewidthwin}{
There are two options to estimate the rolling window correlations and to plot the heat map: (1) typewidthwin=``PARTIAL'' to plot a band of windows from widthwin_1 to widthwin_N (both parameters were defined previously), or (2) typewidthwin=``FULL'' to plot all the  possible window-lengths (from 3 to dim(inputdata)[1]), please note that typewidthwin by default is FULL. 
}
 \item{widthwin_1}{
First value for the size (length) of the windows when the option typewidthwin=``PARTIAL'' is selected, by default is 3, but other values can be used. 
}
 \item{widthwin_N}{
Last value for the size (length) of the windows when the option typewidthwin=``PARTIAL'' is selected, by default is dim(inputdata)[1], but other widthwin_N values can be used (please note that widthwin_1 < widthwin_N).
}
  \item{Align}{
To align the rolling object, RolWinMulCor ONLY uses the ``center'' option (please look at: \R>?running)  to ensure that variations in the correlation are aligned with the variations in the relationship of the time series under study rather than being shifted.
}
 \item{pvalcorectmethod}{
The p-value correction method to be used, by default the method of Benjamini & Hochberg (BH) is used since this is less conservative and performs much better than Bonferroni, but other six methods (e.g. Holm, Bonferroni, etc.) are available (please look at: \R>?p.adjust).
}
 \item{rmltrd}{
Remove (by default is ``Y'' or ``y''; please use ``N'' or ``n'' otherwise) the linear trend in the time series under analysis. 
}
 \item{Scale}{
Scale (by default is ``Y'' or ``y''; please use ``N'' or ``n'' otherwise) is used to ``normalize'' or ``standardize''the time series under analysis.
}
 \item{device}{
Kind of plot output (please look at: \R>?device), there are five options: ``png'', ``jpeg/jpg'', ``eps'', ``pdf'' and ``screen.'' By default device is ``screen.''
}
 \item{Hfig, Wfig}{
Plot's height and width (for the device) for ``png'' and ``jpg'' format (look at \R>?png or \R>?jpg), by default Hfig and Wfig have values of 900, but other values can be used. 
}
  \item{resfig}{
Image resolution (in ``ppi'') for the plot in ``png'' and ``jpg'' format (look at \R>?png , \R>?jpg or \R>?jpeg), by default this has a value of 150, but other mage resolutions can be used. 
}
 \item{Hpdf, Wpdf}{
Plot's height and width (for the device) for ``pdf'' or ``eps'' format (>\R?pdf or >\R?postscript), by default Hpdf and Wpdf are equal to 7, but it is possible to use other values. 
}
 \item{NUMLABX}{
Number of labels for (all) the X's axis, by the default is 5, but it is possible to use other values. 
}
 \item{ofilename}{
Output file name. 
}
}
\details{
The \code{\link{rolwinmulcor_heatmap}} function estimates the rolling window correlation between multiple time series (multi-variate case) sampled on identical time points for all the possible window-lengths or for a band of window-lengths and plots the rolling correlation coefficients and their respective p-values as a heat map. \code{\link{rolwinmulcor_heatmap}} uses the functions \emph{running} (package:gtools), the function \emph{cor.def} that we have created for our \emph{RolWinMulCor} \R package, and the native \R function \emph{p.adjust} (package:stats). 
}
\value{
Outputs:

Plot output (heat map): \emph{screen} or \emph{'heatmap_multivariate_' + 'ofilename + .png, .jpg, .eps} or \emph{.pdf'}. 

Numerical output: two lists \emph{matcor} and \emph{pvalscor} containing the correlation matrix and their corresponding corrected p-values, \emph{NoWindows} and \emph{Windows} that contains the number of windows and the window-lengths. 
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
DeustoTech - Deusto Institute of Technology, \cr
Faculty of Engineering, University of Deusto, \cr
Avda. Universidades, 24, Bilbao, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{http://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com}, \email{josue.polanco@deusto.es}\cr
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
 # Testing the function rolwinmulcor_heatmap
 # typewidthwin="PARTIAL," window lengths from 21 to 31 and plot output in screen format 
 test1_rolwinmulcor_heatmap <- rolwinmulcor_heatmap(synthetic_data, varnametsY="Y", 
                                varnametsX=paste("X1", "X2", "X3", sep=", "), 
                                coltsX = c("red", "blue", "green"), units="NU", 
                                typewidthwin="PARTIAL", widthwin_1=21, widthwin_N=31, 
                                device="screen", ofilename="test1")
 # This example could takes a long time! 
 # typewidthwin="FULL" and pot output in PDF format 
 \donttest{
 test2_rolwinmulcor_heatmap <- rolwinmulcor_heatmap(synthetic_data, varnametsY="Y", 
                                varnametsX=paste("X1", "X2", "X3", sep=", "), 
                                coltsX = c("red", "blue", "green"), units="NU", 
                                typewidthwin="FULL", device="pdf", ofilename="test2")

 } 
}
\keyword{rolwincor_multivariate}
\keyword{rolwincor_multivariate_heatmap}

