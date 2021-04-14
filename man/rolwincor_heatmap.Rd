\name{rolwincor_heatmap}
\alias{rolwincor_heatmap}
\alias{rolwincor_heatmap}
\title{Estimate the Rolling Window Correlation for the bi-variate case to plot its outputs as a heat map 
}
\description{
The \code{\link{rolwincor_heatmap}} function estimates the rolling window correlation coefficients and their respective p-values between TWO time series (bi-variate case) sampled on identical time points for all the possible (i.e. from 3 to the number of elements of the time series under analysis) window-lengths (time-scales) or for a band of window-lengths to be plotted as a heat map. To carry out the computational implementation we extend the works of Telford (2013), Polanco-Martínez (2019), and Polanco-Martínez (2020). The \code{\link{rolwincor_heatmap}} function is highly flexible since this contains several parameters to control the estimation of correlation. A list of parameters are described in the following lines. 
}
\usage{
rolwincor_heatmap(inputdata, varX="", varY="", CorMethod="pearson", 
                  typewidthwin="FULL", widthwin_1=3, 
                  widthwin_N=dim(inputdata)[1], Align="center", 
                  pvalcorectmethod="BH", rmltrd=TRUE, Scale=TRUE)
}
\arguments{
 \item{inputdata}{
Matrix of 3 columns: time, first variable (e.g. \eqn{X}), and second variable (e.g. \eqn{Y}). Please note that time's resolution can be annual, seasonal, monthly, etc, but the format must be in the following form: e.g. for monthly data, year 1 must be repeated 12 times (thus, these times will correspond to 12 rows), year 2 must be also repeated 12 times, and so on. 'Times' must be regular (equally spaced in time and not missing data).
} 
 \item{varX, varY}{
Names of the first (e.g. \eqn{X}) and second (e.g. \eqn{Y}) variable. Please note that the names of these two variables MUST be defined.
} 
 \item{CorMethod}{
The method used to estimate the correlations, by default is ``pearson'' but other options (``spearman'' and ``kendall'') are available (please look at: \R>?cor.test).
}
 \item{typewidthwin}{``FULL'' is to estimate the windows from 2, 4, ..., to dim(inputdata)[1]) if \code{Align} is equal to ``left'' or ``right'', or from 3, 5,..., to dim(inputdata)[1]) if Align is ``center''. The other option is ``PARTIAL'', please you should take into account that \code{widthwin_1} and \code{widthwin_1} MUST be ODD if the \code{Align} option is ``center''. 
}
 \item{widthwin_1}{
First value for the size (length) of the windows when the option \code{typewidthwin}=\cr 
``PARTIAL'' is selected, the minimum value is 3 (the default value), but you must define this parameter (please note that \code{widthwin_1} < \code{widthwin_N}). 
}
 \item{widthwin_N}{
Last value for the size (length) of the windows when the option \code{typewidthwin}=\cr 
``PARTIAL'' is selected, by default is dim(inputdata)[1], but you must define this parameter (please note that \code{widthwin_1} < \code{widthwin_N}).
}
  \item{Align}{
 To align the rolling object, RolWinMulCor uses three options: ``left'', ``center'', and ``right'' (please look at: \R>?running). However, there are some restrictions, which have been described lines above. We recommend to use the ``center'' option to ensure that variations in the correlations are aligned with the variations in the relationships of the variables under study, rather than being shifted to left or right (Polanco-Martínez 2019, 2020), but this imply that the window-lengths MUST be ODD. 
}
 \item{pvalcorectmethod}{
The p-value correction method to be used, by default the method of Benjamini and Hochberg (BH) (1995) is used since this is less conservative and performs much better than Bonferroni, but other five methods (Holm, Hochberg, Bonferroni, Hommel, and Benjamini and Yekutieli) are available (please look at: \R>?p.adjust). Moreover, \code{pvalcorectmethod} admits a pass-through option \cr 
named ``none'' (p-values will not be corrected).
}
 \item{rmltrd}{
Remove (by default is ``TRUE''; ``FALSE'' otherwise) the linear trend in the time series under analysis. 
}
 \item{Scale}{
Scale (by default is ``TRUE''; ``FALSE''  otherwise) is used to ``normalize'' or ``standardize'' the time series under analysis.
}
}
\details{
The \code{\link{rolwincor_heatmap}} function estimates the rolling window correlation between TWO time series (bi-variate case) sampled on identical time points for all the possible (i.e. from 3 to the number of elements of the time series under analysis) window-lengths (time-scales) or for a band of window-lengths to be plotted the rolling correlation coefficients and their respective p-values as a heat map. \code{\link{rolwincor_heatmap}} uses the functions \emph{running} (package:gtools), the native \R functions \emph{cor}, \emph{cor.test}, and \emph{p.adjust} (package:stats), and some pieces of code written specifically to our \R RolWinMulCor package. 
}
\value{
Outputs:

Numerical output: three lists \emph{matcor}, \emph{pvalscor}, and \emph{pvalNOTcor} containing the correlation matrix and their corresponding corrected and not corrected p-values, as well as \emph{NoWindows} and \emph{Windows} that contains the number of windows and the window-lengths (time-scales), and \emph{CorMethod}, \emph{left_win}, and \emph{righ_win}, which have been previously described.
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
BC3 - Basque Centre for Climate Change, Bilbao, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue-Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com}, \email{josue.polanco@bc3research.org}.
}
\references{
Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57 (1), 289-300. <URL: \doi{10.1111/j.2517-6161.1995.tb02031.x}>. \cr

Polanco-Martínez, J. M. (2019). Dynamic relationship analysis between NAFTA stock markets using nonlinear, nonparametric, non-stationary methods. Nonlinear Dynamics, 97(1), 369-389. <URL: \doi{10.1007/s11071-019-04974-y}>. \cr 

Polanco-Martínez, J. M. (2020). RolWinMulCor : an R package for estimating rolling window multiple correlation in ecological time series. Ecological Informatics, 60, 101163. <URL: \doi{10.1016/j.ecoinf.2020.101163}>. \cr

Telford, R.: Running correlations -- running into problems (2013). <URL: \cr 
\url{https://quantpalaeo.wordpress.com/2013/01/04/}>. 
}
\examples{
\donttest{
# Testing the function rolwincor_heatmap
# typewidthwin="PARTIAL," window lengths from 11 to 101 and plot output in screen format 
# Code used to produce the outputs to generate the Fig. 4 in Polanco-Martínez (2020). 
test_rolwincor_heatmap <- rolwincor_heatmap(syntDATA, varX="X", varY="Y", 
                            CorMethod="spearman", typewidthwin="PARTIAL", widthwin_1=11, 
                            widthwin_N=101, Align="center", pvalcorectmethod="BH", 
                            rmltrd=TRUE, Scale=TRUE)
# This example could takes a long time since typewidthwin="FULL" 
test2_rolwincor_heatmap <- rolwincor_heatmap(syntDATA, varX="X", varY="Y",
                             CorMethod="spearman", typewidthwin="FULL", Align="center", 
                             pvalcorectmethod="BH", rmltrd=TRUE, Scale=TRUE)
 }
}
\keyword{rolwincor_bivariate}
\keyword{rolwincor_bivariate_heatmap}
