\name{rolwincor_1win}
\alias{rolwincor_1win}
\alias{rolwincor_1win}
\title{Estimate the Rolling Window Correlation for the bi-variate case to be plotted its outputs as a single one window 
}
\description{
The \code{\link{rolwincor_1win}} function estimates (correlation coefficients and their respective p-values) the rolling (running) window correlation between TWO time series (bi-variate case) sampled on identical time points for ONLY ONE window-length (time-scale). To carry out the computational implementation we follow to Telford (2013), Polanco-Martínez (2019) and Polanco-Martínez (2020). The \code{\link{rolwincor_1win}} function is highly flexible since this contain several parameters to control the estimation of correlation. For example, \code{\link{rolwincor_1win}} function contain parameters to remove the (linear) trend contained in the time series under analysis or to choose different p-value correction methods (which are used to address the multiple comparison problem). A list of parameters are described in the following lines. 
}
\usage{
rolwincor_1win(inputdata, varX="", varY="", CorMethod="pearson", widthwin=3,
               Align="center", pvalcorectmethod="BH", rmltrd=TRUE, Scale=TRUE)
}
\arguments{
 \item{inputdata}{
Matrix of 3 columns: time, first variable (e.g. \eqn{X}), and second variable (e.g. \eqn{Y}). 
} 
 \item{varX, varY}{
Names of the first (e.g. \eqn{X}) and second (e.g. \eqn{Y}) variable. Please note that the names of these variables MUST be defined.
} 
 \item{CorMethod}{
The method used to estimate the correlations, by default is ``pearson,'' but other options (``spearman'' and ``kendall'') are available (please look at: \R>?cor.test).
}
 \item{widthwin}{
Window size to compute the rolling window correlations. This value can be an even or odd number of at least three (the default value), and this parameter MUST be provided. 
}
 \item{Align}{
To align the rolling object, RolWinMulCor ONLY uses the ``center'' option by default (please look at: \R>?running) to ensure that variations in the correlation are aligned with the variations in the relationship of the time series under study rather than being shifted (Polanco-Martínez, 2019; 2020), but the ``left'' and ``right'' options can be used, but if \emph{widthwin} is an even number it will not be possible to use the ``center'' option (please look at: \R>?running).
}
 \item{pvalcorectmethod}{
The p-value correction method to be used, by default the method of Benjamini and Hochberg (BH) (1995) is used since this is less conservative and performs much better than Bonferroni, but other five methods (Holm, Hochberg, Bonferroni, Hommel, and Benjamini and Yekutieli) are available (please look at: \R>?p.adjust). Moreover, \emph{pvalcorectmethod} admits a pass-through option named ``none'' (p-values will not be corrected).
}
 \item{rmltrd}{
Remove (by default is ``TRUE''; ``FALSE'' otherwise) the linear trend in the two time series under analysis. 
}
 \item{Scale}{
Scale (by default is ``TRUE''; ``FALSE'' otherwise) is used to ``normalize'' or ``standardize'' the time series under analysis.
}
}
\details{
The \code{\link{rolwincor_1win}} function estimates the rolling window correlation between TWO time series (bi-variate case) sampled on identical time points for ONLY ONE window-length (time-scale) and plots the rolling correlation coefficients and their respective p-values. \code{\link{rolwincor_1win}} uses the functions \emph{running} (package:gtools), the native \R functions \emph{cor}, \emph{cor.test}, and \emph{p.adjust} (package:stats), and some pieces of code written specifically to our \R RolWinMulCor package. 
}
\value{
Outputs:

Numerical output: a list containing \emph{Correlation_coefficients}, \emph{P_values_corrected}, and \emph{P_values_not_corrected}, which are self-explanatory, as well as \emph{CorMethod}, \emph{left_win}, \emph{righ_win}, and \emph{widthwin}, which indicate the method used to estimate the correlations, first and last time element of the rolling correlation matrix, and the window-length (time-scale).
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
# Testing the function rolwincor_1win (bi-variate case)
# Window length = 21 and plot output in screen format 
test_rolwincor_1win <- rolwincor_1win(syntDATA, varX="X", varY="Y", CorMethod="spearman", 
                                      widthwin=21, Align="center", pvalcorectmethod="BH")
}
\keyword{rolwincor_bivariate}
\keyword{rolwincor_bivariate_simply_plot}
