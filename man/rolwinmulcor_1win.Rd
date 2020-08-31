\name{rolwinmulcor_1win}
\alias{rolwinmulcor_1win}
\alias{rolwinmulcor_1win}
\title{Estimate the Rolling Window Correlation for the multi-variate case to be plotted the outputs as a single one window 
}
\description{
The \code{\link{rolwinmulcor_1win}} function estimates the rolling window correlation coefficients and their respective p-values among multiple time series (multi-variate case) sampled on identical time points for ONLY ONE window-length (time-scale). The multi-variate case is based on the concept of multiple regression and generalizes the standard coefficient of correlation (the squared multiple correlation coefficient, or more appropriate, the adjusted coefficient of determination, \eqn{R^2}). We follow to Abdi (2007) to implement computationally this technique. The \code{\link{rolwinmulcor_1win}} function is highly flexible since this contains several parameters to control the estimation of correlation. For example, \code{\link{rolwinmulcor_1win}} function contains parameters to remove the (linear) trend contained in the time series under analysis, to choose different p-value correction methods (which are used to address the multiple comparison problem). A list of parameters are described in the following lines. 
}
\usage{
rolwinmulcor_1win(inputdata, varnametsY="", varnametsX="", widthwin=5, 
                  Align="center", pvalcorectmethod="BH", 
                  rmltrd="TRUE", Scale="TRUE")
}
\arguments{
 \item{inputdata}{
Matrix of P columns: time, dependent variable (\eqn{Y}), and independent variables (\eqn{X_1, X_2,..., X_{P-2}}).
} 
 \item{varnametsY}{
Name of the dependent variable: \eqn{Y}. Please note that the name of this variable MUST be defined. 
} 
 \item{varnametsX}{
Name of the independent variables: \eqn{X_1, X_2,..., X_{p-2}}. Please note that the names of these variables MUST defined in this way: varnametsX=paste("X1", "X2",..., sep=", ").
}
 \item{rmltrd}{
Remove (by default is ``TRUE''; ``FALSE'' otherwise) the linear trend in the time series under analysis. 
}
 \item{Scale}{
Scale (by default is ``TRUE''; ``FALSE'' otherwise) is used to ``normalize'' or ``standardize'' the time series under analysis.
}
 \item{widthwin}{
Window size to compute the rolling window correlations. This value can be an even or odd number of at least three (the default value), and this parameter MUST be provided. 
}
 \item{Align}{
To align the rolling object, RolWinMulCor ONLY uses the ``center'' option by default (please look at: \R>?running) to ensure that variations in the correlation are aligned with the variations in the relationship of the time series under study rather than being shifted (Polanco-Martínez, 2019; 2020), but the ``left'' and ``right'' options can be used, but if widthwin is an even number it will not be possible to use the ``center'' option (please look at: \R>?running).
}
 \item{pvalcorectmethod}{
The p-value correction method to be used, by default the method of Benjamini and Hochberg (BH) (1995) is used since this is less conservative and performs much better than Bonferroni, but other five methods (Holm, Hochberg, Bonferroni,Hommel, and Benjamini and Yekutieli) are available (please look at: \R>?p.adjust). Moreover, \emph{pvalcorectmethod} admits a pass-through option named ``none'' (p-values will not be corrected).
}
}
\details{
The \code{\link{rolwinmulcor_1win}} function estimates the rolling window correlation coefficients and their respective p-values among multiple time series (multi-variate case) sampled on identical time points for ONLY ONE window-length. \code{\link{rolwinmulcor_1win}} uses the functions \emph{rollapply} (package:zoo) that is able to tackle matrices, the native \R function \emph{p.adjust} (package:stats), and some pieces of code and an auxiliary function that we have created specifically for our function \code{\link{rolwinmulcor_1win}} and RolWinMulCor \R package. 
}
\value{
Outputs:

Numerical output: three list containing \emph{Correlation_coefficients}, \emph{P_values_corrected}, and \emph{P_values_not_corrected}, which are self-explanatory, as well as \emph{left_win}, \emph{righ_win}, and \emph{widthwin}, which indicate the method used to estimate the correlations, first and last time element of the rolling correlation matrix, and the window-length (time-scale).
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
Abdi H. Multiple correlation coefficient, in Encyclopedia of Measurement and
Statistics, N. J. Salkind, Ed. Sage, Thousand Oaks, CA, USA, 2007; 648-651. \cr
<URL: \url{https://personal.utdallas.edu/~herve/Abdi-MCC2007-pretty.pdf}>. \cr

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57 (1), 289-300. <URL: \url{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1995.tb02031.x}>. \cr

Polanco-Martínez, J. M. (2019). Dynamic relationship analysis between NAFTA stock markets using nonlinear, nonparametric, non-stationary methods. Nonlinear Dynamics, 97(1), 369-389. <URL: \url{https://doi.org/10.1007/s11071-019-04974-y}>. \cr 

Polanco-Martínez, J. M. (2020). RolWinMulCor : an R package for estimating rolling window multiple correlation in ecological time series. Ecological Informatics (Ms. ECOINF-D-20-00263 accepted for publication, 19/08/2020).\cr
}
\examples{
 # Testing the function rolwinmulcor_1win (multi-variate case)
 test_rolwinmulcor_1win <- rolwinmulcor_1win(YX_ecological_data, widthwin=61,
                            Align="center", pvalcorectmethod="BH")
}
\keyword{rolwincor_multivariate}
\keyword{rolwincor_multivariate_simply_red}
