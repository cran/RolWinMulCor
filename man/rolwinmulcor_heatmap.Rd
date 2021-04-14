\name{rolwinmulcor_heatmap}
\alias{rolwinmulcor_heatmap}
\alias{rolwinmulcor_heatmap}
\title{Estimate the Rolling Window Correlation for the multi-variate case to plot the results as a heat map 
}
\description{
The \code{\link{rolwinmulcor_heatmap}} function estimates the rolling window correlation coefficients and their respective p-values among several time series (multi-variate case) sampled on identical time points for all the possible (i.e. from 3 to the number of elements of the time series under analysis) window-lengths (time-scales) or for a band of window-lengths, and the outputs are used to produce a heat map. The multivariate case is based on the concept of multiple regression and generalizes the standard coefficient of correlation (the squared multiple correlation coefficient, or more appropriate, the adjusted coefficient of determination, \eqn{R^2}). We follow (and extend) the work of Abdi (2007) and Polanco-Martínez (2020) to implement computationally this technique. The \code{\link{rolwinmulcor_heatmap}} function is highly flexible since this contains several parameters to control the estimation of correlation and features of the plot output. A list of parameters are described in the following lines. 
}
\usage{
rolwinmulcor_heatmap(inputdata, varnametsY="", varnametsX="", 
                     typewidthwin="FULL", widthwin_1=7, 
                     widthwin_N=dim(inputdata)[1], Align="center", 
                     pvalcorectmethod="BH", rmltrd=TRUE, Scale=TRUE)
}
\arguments{
 \item{inputdata}{
Matrix of P columns: time, dependent variable (\eqn{Y}), and independent variables (\eqn{X_1, X_2,..., X_{P-2}}). Please note that time's resolution can be annual, seasonal, monthly, etc, but the format must be in the following form: e.g. for monthly data, year 1 must be repeated 12 times (thus, these times will correspond to 12 rows), year 2 must be also repeated 12 times, and so on. 'Times' must be regular (equally spaced in time and not missing data).
} 
 \item{varnametsY}{
Name of the dependent variable: \eqn{Y}. Please note that the name of this variable MUST be defined. 
} 
 \item{varnametsX}{
Name of the independent variables: \eqn{X_1, X_2,..., X_{p-2}}. Please note that the names of these variables MUST be defined in this way: \code{varnametsX}=c("X1", "X2",..., sep=", ").
}
\item{rmltrd}{
Remove (by default is ``TRUE''; ``FALSE'' otherwise) the linear trend in the time series under analysis. 
}
 \item{Scale}{
Scale (by default is ``TRUE''; ``FALSE'' otherwise) is used to ``normalize'' or ``standardize''the time series under analysis.
}
 \item{typewidthwin}{
`FULL'' is to estimate the windows from 2, 4, ..., to dim(inputdata)[1]) if \code{Align} is equal to ``left'' or ``right'', or from 3, 5,..., to dim(inputdata)[1]) if Align is ``center''. The other option is ``PARTIAL'', please you should take into account that \code{widthwin_1} and \code{widthwin_1} MUST be ODD if the \code{Align} option is ``center''.
}
 \item{widthwin_1}{
First value for the size (length) of the windows when the option \emph{typewidthwin}=``PARTIAL'' is selected, the minimum value is 3 (the default value), but you must define this parameter (please note that \code{widthwin_1} < \code{widthwin_N}).
}
 \item{widthwin_N}{
Last value for the size (length) of the windows when the option \emph{typewidthwin}=``PARTIAL'' is selected, by default is dim(inputdata)[1], but you must define this parameter (please note that \code{widthwin_1} < \code{widthwin_N}).
}
 \item{Align}{
 To align the rolling object, RolWinMulCor uses three options: ``left'', ``center'', and ``right'' (please look at: \R>?running). However, there are some restrictions, which have been described lines above. We recommend to use the ``center'' option to ensure that variations in the correlations are aligned with the variations in the relationships of the variables under study, rather than being shifted to left or right (Polanco-Martínez 2019, 2020), but this imply that the window-lengths (time-scales) MUST be ODD. 
}
 \item{pvalcorectmethod}{
The p-value correction method to be used, by default the method of Benjamini and Hochberg (BH) (1995) is used since this is less conservative and performs much better than Bonferroni, but other five methods (Holm, Hochberg, Bonferroni,  and Benjamini and Yekutieli) are available (please look at: \R>?p.adjust). Moreover, \code{pvalcorectmethod} admits a pass-through option named ``none'' (p-values will not be corrected).
}
}
\details{
The \code{\link{rolwinmulcor_heatmap}} function estimates the rolling window correlation coefficients and their respective p-values between multiple time series (multi-variate case) sampled on identical time points for all the possible window-lengths (time-scales) or for a band of window-lengths. \code{\link{rolwinmulcor_heatmap}} uses the functions \emph{rollapply} (package:zoo) that is able to tackle matrices, the native \R function \emph{p.adjust} (package:stats), and some pieces of code and an auxiliary function that we have created specifically for our function \code{\link{rolwinmulcor_heatmap}} and \R RolWinMulCor package. 
}
\value{
Outputs:

Numerical output: three lists \emph{matcor}, \emph{pvalscor}, and \emph{pvalNOTcor} containing the correlation coefficients and their corresponding corrected and not corrected p-values, as well as \emph{NoWindows}, \emph{Windows} that contains the number of windows and the window-lengths (time-scales), and \emph{CorMethod}, \emph{left_win}, and \emph{righ_win} that have been previously described.
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
BC3 - Basque Centre for Climate Change, Bilbao, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue-Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com}, \email{josue.polanco@bc3research.org}.
}
\references{
Abdi H. Multiple correlation coefficient, in Encyclopedia of Measurement and
Statistics, N. J. Salkind, Ed. Sage, Thousand Oaks, CA, USA, 2007; 648-651. \cr
<URL: \url{https://personal.utdallas.edu/~herve/Abdi-MCC2007-pretty.pdf}>. \cr

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57 (1), 289-300. <URL: \doi{10.1111/j.2517-6161.1995.tb02031.x}>. \cr

Polanco-Martínez, J. M. (2019). Dynamic relationship analysis between NAFTA stock markets using nonlinear, nonparametric, non-stationary methods. Nonlinear Dynamics, 97(1), 369-389. <URL: \doi{10.1007/s11071-019-04974-y}>. \cr 

Polanco-Martínez, J. M. (2020). RolWinMulCor: an R package for estimating rolling window multiple correlation in ecological time series. Ecological Informatics, 60, 101163. <URL: \doi{10.1016/j.ecoinf.2020.101163}>. 
}
\examples{
\donttest{
# Testing the function rolwinmulcor_heatmap
test_rolwinmulcor_heatmap <- rolwinmulcor_heatmap(YX_ecological_data, 
                               typewidthwin="PARTIAL", widthwin_1=11, widthwin_N=101, 
                               Align="center", pvalcorectmethod="BH")
}
}
\keyword{rolwincor_multivariate}
\keyword{rolwincor_multivariate_heatmap}

