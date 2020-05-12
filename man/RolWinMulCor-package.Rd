\name{RolWinMulCor-package}
\alias{RolWinMulCor-package}
\alias{RolWinMulCor}
\docType{package}
\title{Estimate the Rolling Window Multiple Correlation 
}
\description{
'RolWinMulCor' estimates the rolling (running) window correlation for 
the bi- and multi-variate cases between regular (sampled on identical 
time points) time series, with especial emphasis to environmental data 
(although this can be applied to other kinds of data sets). 
'RolWinMulCor' is based on the concept of rolling or running window 
correlation and is useful to evaluate the evolution of correlation 
through time and time-scales. 'RolWinMulCor' contains four functions: 
(1) the first two are focused on the bi-variate case, one of them 
produces a simple plot of correlation coefficients and p-values (<=0.05) 
for only one window-length (time-scale) and the other function produces a 
heat map for the statistically significant (p-values <=0.05) correlation 
coefficients taking into account all the possible window-lengths (which 
are determined by the number of elements of the time series under 
analysis) or for a band of window-lengths; (2) the second two functions 
are designed to analyse the multi-variate case and follows the bi-variate 
case to display visually the results. The four functions contained in 
'RolWinMulCor' are highly flexible since this contains a great number of 
parameters to control the estimation of correlation and the features of 
the plot output, e.g. to remove the linear trend contained in the time 
series under analysis, to choose different p-value correction methods 
(which are used to address the multiple comparison problem) or to 
personalise the plot output (e.g. this can be displayed in the screen or 
can be saved as PNG, JPG, EPS or PDF formats). We would like to highlight 
that, to the best of our knowledge, there are few R packages (probably
the one) on CRAN that estimate rolling window correlation and produce a 
heat map for the bi-variate case and, especially, for the multi-variate 
case. The 'RolWinMulCor' package also provides examples with 
synthetic and real environmental time series to exemplify its use. 
}
\details{
\tabular{ll}{
Package: \tab RolWinMulCor\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2020-05-09\cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes\cr
}
RolWinMulCor package contains four functions: (1) \code{\link{rolwincor_1win}} (estimates and plots the rolling window correlation for the bi-variate case for only one window-length or time-scale for the time series under study), (2) \code{\link{rolwincor_heatmap}} (estimates and plots as a heat map the statistically significant (p-values <=0.05) correlation coefficients taking into account all the possible window-lengths that are determined by the number of elements of the time series under or or a band of window-lengths and plot the the correlation coefficients and their respective p-values as a heat map), (3) \code{\link{rolwinmulcor_1win}} (estimates and plots the rolling window correlation for the multi-variate case for only one window-length or time-scale for the time series under study) and, (4) \code{\link{rolwinmulcor_heatmap}} (estimates and plots the heat map for the multi-variate case).
}
\note{
 Dependencies: \emph{stat}, \emph{gtools}, \emph{pracma} and \emph{colorspace}. 
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
DeustoTech - Deusto Institute of Technology, \cr
Faculty of Engineering, University of Deusto, \cr
Avda. Universidades, 24, Bilbao, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{http://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com}, \email{josue.polanco@deusto.es} \cr\cr
Acknowledgement:\cr 
The author acknowledges to the SEPE (Spanish Public Service of Employment) for its funding support. Special thanks to Martina Schmirl and Jelena Saf (from CRAN team) for the helpful comments on the package. 
}
\references{
Abdi H. Multiple correlation coefficient, in Encyclopedia of Measurement and
Statistics, N. J. Salkind, Ed. Sage, Thousand Oaks, CA, USA, 2007; 648-651.\cr
<URL: \url{https://personal.utdallas.edu/~herve/Abdi-MCC2007-pretty.pdf}>. \cr

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57 (1), 289-300. <URL: \url{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1995.tb02031.x}>. \cr

Polanco-Martínez, J. M. (2019). Dynamic relationship analysis between NAFTA stock markets using nonlinear, nonparametric, non-stationary methods. Nonlinear Dynamics, 97(1), 369-389. <URL: \url{https://doi.org/10.1007/s11071-019-04974-y}>. \cr 

Telford, R.: Running correlations -- running into problems (2013). <URL: \cr 
\url{https://quantpalaeo.wordpress.com/2013/01/04/running-correlations-running-into-problems/}>. \cr
}
