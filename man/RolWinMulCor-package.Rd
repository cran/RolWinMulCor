\name{RolWinMulCor-package}
\alias{RolWinMulCor-package}
\alias{RolWinMulCor}
\docType{package}
\title{Estimate the Rolling Window Multiple Correlation 
}
\description{
'RolWinMulCor' estimates the rolling (running) window correlation for 
the bi- and multi-variate cases between regular (sampled on identical 
time points) time series, with especial emphasis to ecological data 
although this can be applied to other kinds of data sets. 'RolWinMulCor' 
is based on the concept of rolling, running, or sliding window correlation 
and is useful to evaluate the evolution of correlation through time and 
time-scales. 'RolWinMulCor' contains six (four for estimations and two 
for plots) functions. The first two functions focus on the bi-variate 
case: (1) rolwincor_1win and (2) rolwincor_heatmap, estimate the
correlation coefficients and their respective p-values for only one 
window-length (time-scale) and considering all possible window-lengths 
or a band of window-lengths, respectively. The second two functions: (3) 
rolwinmulcor_1win and (4) rolwinmulcor_heatmap, are designed to analyze 
the multi-variate case, following the bi-variate case to visually display 
the results, but these two approaches are methodologically different: 
the multi-variate case estimate the adjusted coefficients of determination 
instead of the correlation coefficients. The last two functions: (5) 
plot_1win and (6) plot_heatmap, are used to represent graphically the 
outputs of the four aforementioned functions as simple plots or as heat 
maps. The six functions contained in 'RolWinMulCor' are highly flexible 
since this contains several parameters to control the estimation of 
correlation and the features of the plot output, e.g. to remove the linear 
trend contained in the time series under analysis, to choose different 
p-value correction methods (which are used to address the multiple 
comparison problem) or to personalise the plot output. The 'RolWinMulCor' 
package also provides examples with synthetic and real-life ecological 
time series to exemplify its use. 
}
\details{
\tabular{ll}{
Package: \tab RolWinMulCor\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2020-09-05\cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes\cr
}
RolWinMulCor package contains six functions: (1) \code{\link{rolwincor_1win}} estimates the rolling window correlation coefficients and their respective p-values for the bi-variate case for only one window-length or time-scale for the time series under study, (2) \code{\link{rolwincor_heatmap}} estimates the correlation coefficients and their corresponding p-values taking into account all the possible window-lengths that are determined by the number of elements of the time series under analysis or a band of window-lengths, (3) \code{\link{rolwinmulcor_1win}} estimates the rolling window correlation coefficients and their p-values for the multi-variate case for only one window-length or time-scale for the time series under study, (4) \code{\link{rolwinmulcor_heatmap}} estimates the correlation coefficients and their corresponding p-values for the multi-variate case taking into account all the possible window-lengths or a band of window-lengths, (5) \code{\link{plot_1win}} plots the correlation coefficients and their respective p-values (corrected or not corrected) as only one selected window-length using the outputs of the functions \code{\link{rolwincor_1win}} (bi-variate case) and \code{\link{rolwinmulcor_1win}} (multi-variate case), and (6) \code{\link{plot_heatmap}} plots the heat maps for the correlation coefficients and their respective p-values (corrected or not corrected) for all possible window-lengths (i.e., from five to the number of elements in the time series under analysis) or for a band of window-lengths using the outputs of the functions \code{\link{rolwincor_heatmap}} (bi-variate case) and \code{\link{rolwinmulcor_heatmap}} (multi-variate case). The bi-variate case follow from a methodological point of view to Telford (2013), Polanco-Martínez (2019), and Polanco-Martínez (2020) whereas the multi-variate case follow to Abdi (2007) and  Polanco-Martínez (2020). 
}
\note{
 Dependencies: \emph{stat}, \emph{gtools}, \emph{zoo}, \emph{pracma} and \emph{colorspace}. 
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
DeustoTech - Deusto Institute of Technology, \cr
Faculty of Engineering, University of Deusto, \cr
Avda. Universidades, 24, Bilbao, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com}, \email{josue.polanco@deusto.es} \cr\cr
Acknowledgement:\cr 
The author acknowledges to the SEPE (Spanish Public Service of Employment) for its funding support. Special thanks to the CRAN team (in particular to Martina Schmirl and Jelena Saf), to Ana-Maria Hereş and Jorge Curie for their helpful comments on the package, and for the three reviewers (this package is described in the paper recently accepted for publication that is cited in the References), in particular Reviewer #2, that provided some very useful suggestions to improve RolWinMulCor.  
}
\references{
Abdi H. Multiple correlation coefficient, in Encyclopedia of Measurement and
Statistics, N. J. Salkind, Ed. Sage, Thousand Oaks, CA, USA, 2007; 648-651.\cr
<URL: \url{https://personal.utdallas.edu/~herve/Abdi-MCC2007-pretty.pdf}>. \cr

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57 (1), 289-300. <URL: \url{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1995.tb02031.x}>. \cr

Polanco-Martínez, J. M. (2019). Dynamic relationship analysis between NAFTA stock markets using nonlinear, nonparametric, non-stationary methods. Nonlinear Dynamics, 97(1), 369-389. <URL: \url{https://doi.org/10.1007/s11071-019-04974-y}>. \cr 

Polanco-Martínez, J. M. (2020). RolWinMulCor : an R package for estimating rolling window multiple correlation in ecological time series. Ecological Informatics (Ms. ECOINF-D-20-00263 accepted for publication, 19/08/2020).\cr

Telford, R.: Running correlations -- running into problems (2013). <URL: \cr 
\url{https://quantpalaeo.wordpress.com/2013/01/04/running-correlations-running-into-problems/}>. \cr
}
