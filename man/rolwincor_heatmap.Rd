\name{rolwincor_heatmap}
\alias{rolwincor_heatmap}
\alias{rolwincor_heatmap}
\title{Estimate the Rolling Window Correlation for the bi-variate case and plot the results as a heat map 
}
\description{
The \code{\link{rolwincor_heatmap}} function estimates the rolling (running) window correlation between TWO time series (bi-variate case) sampled on identical time points for all the possible (i.e. from 3 to the number of elements of the time series under analysis) window-lengths or for a band of window-lengths and plots the correlation coefficients and their respective p-values (<= 0.05) as a heat map. The \code{\link{rolwincor_heatmap}} function is highly flexible since this contains a great number of parameters to control the estimation of correlation and features of the plot output. A list of parameters are described in the following lines. 
}
\usage{
rolwincor_heatmap(inputdata, varnametsX="", varnametsY="", units="", 
                  coltsX="black", coltsY="blue", CEXLAB=1.15, CEXAXIS=1.5, 
                  LWDtsX=1, LWDtsY=1, CorMethod="spearman", typewidthwin="FULL", 
                  widthwin_1=3, widthwin_N=dim(inputdata)[1],  
                  Align="center", pvalcorectmethod="BH", rmltrd="Y", 
                  Scale="Y", device="screen", Hfig=900, Wfig=900, 
                  resfig=150, Hpdf=7, Wpdf=7, NUMLABX=5, ofilename)
}
\arguments{
 \item{inputdata}{
Matrix of 3 columns: time, first variable (``X''), and second variable (``Y''). 
} 
 \item{varnametsX, varnametsY}{
Names of the first (e.g. ``X'') and second (e.g. ``Y'') variable. Please note that the names of these two variables MUST be defined.
} 
 \item{units}{
Time's unit (e.g. days, weeks, years, etc.) for the variables under analysis. Please note that the units MUST be defined. 
}
 \item{coltsX, coltsY}{
The colors to be used to plot the two variables under study, by default the colors are ``black'' and ``blue,'' respectively, but it is possible to use other colors.
}
 \item{CEXLAB, CEXAXIS}{
These parameters are used to plot the sizes of X-axis and Y-axis labels and X- and Y-axis, by default these parameters have values of 1.15 and 1.5, respectively, but it is possible to use other values. 
}
 \item{LWDtsX, LWDtsY}{
These parameters are used to define the line-widths when the variables are plotted, by default these have values of 1, but other values (widths) can be used. 
}
 \item{CorMethod}{
The method used to estimate the correlations, by default is ``spearman'' but other options (``pearson'' and ``kendall'') are available (please look at: \R>?cor.test).
}
 \item{typewidthwin}{
There are two options to estimate the rolling correlations and to plot the heat map: (1) typewidthwin=``PARTIAL'' to plot a band of windows from widthwin_1 to widthwin_N (both parameters were previously defined), or (2) typewidthwin=``FULL'' to plot all the  possible window-lengths (from 3 to  dim(inputdata)[1]), by default is FULL. 
}
 \item{widthwin_1}{
First value for the size (length) of the windows when the option typewidthwin=``PARTIAL'' is selected, by default is 3, but others values can be used. 
}
 \item{widthwin_N}{
Last value for the size (length) of the windows when the option typewidthwin=``PARTIAL'' is selected, by default is dim(inputdata)[1], but other values can be used (please note that widthwin_1 < widthwin_N).
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
Plot's height and width (for the device) for ``png'' and ``jpg'' format (look at \R>?png or \R>?jpg), by default Hfig and Wfig have values of 900, but other values can be used. 
}
  \item{resfig}{
Image resolution (in ``ppi'') for the plot in ``png'' and ``jpg'' format (look at \R>?png , \R>?jpg or \R>?jpeg), by default this has a value of 150, but other image resolutions can be used. 
}
 \item{Hpdf, Wpdf}{
Plot's height and width (for the device) for ``pdf'' or ``eps'' format (>\R?pdf or >\R?postscript), by default Hpdf and Wpdf are equal to 7, but other values can be used. 
}
 \item{NUMLABX}{
Number of labels for (all) the X's axis, by the default is 5, but other number of labels can be utilized. 
}
 \item{ofilename}{
Output file name. 
}
}
\details{
The \code{\link{rolwincor_heatmap}} function estimates the rolling window correlation between TWO time series (bi-variate case) sampled on identical time points  for all the possible (i.e. from 3 to the number of elements of the time series under analysis) window-lengths or for a band of window-lengths and plots the rolling correlation coefficients and their respective p-values as a heat map. \code{\link{rolwincor_heatmap}} uses the functions \emph{running} (package:gtools) and the native \R functions \emph{cor}, \emph{cor.test}, and \emph{p.adjust} (package:stats). 
}
\value{
Outputs:

Plot output (heat map): \emph{screen} or \emph{'heatmap_bivariate_' + 'ofilename + .png, .jpg, .eps} or \emph{.pdf'}. 

Numerical output: two lists \emph{matcor} and \emph{pvalscor} containing the correlation matrix and their corresponding corrected p-values, \emph{NoWindows} and \emph{Windows} that contains the number of windows and the window-lengths. 
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
 # Testing the function rolwincor_heatmap
 #  typewidthwin="PARTIAL," window lengths from 21 to 31 and plot output in screen format 
 test1_rolwincor_heatmap <- rolwincor_heatmap(synthetic_data[,1:3], varnametsX="X", 
                             varnametsY="Y", units="NU", device="screen", typewidthwin="PARTIAL", 
                             widthwin_1=21, widthwin_N=31, ofilename="test1")
 # This example could takes a long time! 
 # typewidthwin="FULL" and pot output in PDF format 
 \donttest{
 test2_rolwincor_heatmap <- rolwincor_heatmap(synthetic_data[,1:3], varnametsX="X", 
                             varnametsY="Y", units="NU", device="pdf", typewidthwin="FULL", 
                             ofilename="test2")

 }
}
\keyword{rolwincor_bivariate}
\keyword{rolwincor_bivariate_heatmap}
