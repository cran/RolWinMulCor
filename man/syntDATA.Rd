\name{syntDATA}
\alias{syntDATA}
\title{Synthetic data set to test the functions of \emph{RolWinMulCor}
}
\description{
The data set \code{\link{syntDATA}} contains three columns, the first one is the time and the next three columns are sinusoidal time series that contains two periodical signals (at 11 and 21, with a phase of zero and amplitude of 1 --- please note that these quantities are dimensionless) contaminated by Gaussian noise (with mean of 0 and standard deviation of 0.25) for the intervals 1--100 (signal at 11) and 200--400 (signal at 21) and Gaussian noise (with mean of 0 and standard deviation of 1) otherwise (more information about \code{\link{syntDATA}} in Polanco-Martínez (2020)). 
}
\usage{data(syntDATA)
}
\format{One file in ASCII format and columns are separated by spaces. 
}
\source{
Author's own production (Josué M. Polanco-Martínez). 
}
\references{
Polanco-Martínez, J. M. (2020). RolWinMulCor : an R package for estimating rolling window multiple correlation in ecological time series. Ecological Informatics, 60, 101163. <URL: \doi{10.1016/j.ecoinf.2020.101163}>. \cr
}
\keyword{dataset1}
