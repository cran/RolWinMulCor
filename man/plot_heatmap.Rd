\name{plot_heatmap}
\alias{plot_heatmap}
\alias{plot_heatmap}
\title{Plot the outputs of \code{\link{rolwincor_heatmap}} and  \code{\link{rolwinmulcor_heatmap}} as a heat map
}
\description{
The \code{\link{plot_heatmap}} function plots the correlation coefficients and their respective p-values (corrected or not corrected) as heat maps for all possible window-lengths (i.e., from five to the number of elements in the time series under analysis) or for a band of window-lengths using the outputs of the functions \code{\link{rolwincor_heatmap}} (bi-variate case) and \code{\link{rolwinmulcor_heatmap}} (multi-variate case). The \code{\link{plot_heatmap}} function is highly flexible since this contains several parameters to control the plot output. We would highlight that only the first 12 parameters (and \emph{LWDtsX} and \emph{LWDtsY} must be defined by the users since the others parameters are defined by default. A list of parameters are described in the following lines.}
\usage{
plot_heatmap(inputdata, corcoefs, pvalues, left_win, righ_win, Rwidthwin, 
     KCASE="", typewidthwin="", widthwin_1=3, widthwin_N=dim(inputdata)[1], 
     varX="", varY="", rmltrd=TRUE, Scale=TRUE, coltsX=c("black"), coltsY="blue", 
     CEXLAB=1.15, CEXAXIS=1.05, LWDtsX=1, LWDtsY=1, NUMLABX=5, parcen=c(0.5,25))
}
\arguments{
 \item{inputdata}{
Input data used in the functions \code{\link{rolwincor_heatmap}} or \code{\link{rolwinmulcor_heatmap}}. 
} 
 \item{corcoefs, pvalues}{
Correlation coefficients obtained from the functions \code{\link{rolwincor_heatmap}} or \code{\link{rolwinmulcor_heatmap}} (named \emph{Correlation_coefficients}) and p-values obtained from the aforementioned functions (named \emph{P_values_corrected} or \emph{P_values_not_corrected}).
}
 \item{left_win, righ_win}{
These parameters are used to accommodate the times in the rolling window correlations and are obtained from the functions \code{\link{rolwincor_heatmap}} or \code{\link{rolwinmulcor_heatmap}}, which have the same names.
}
 \item{Rwidthwin}{
Contains the window-sizes where the rolling window correlations are estimated by the functions \code{\link{rolwincor_heatmap}} or by \code{\link{rolwinmulcor_heatmap}}. 
}
 \item{KCASE}{
This parameter is used to activate the cases: ``BIVAR'' for the bi-variate or ``MULVAR'' for the multi-variate, and this must be the same label as the one used in \code{\link{rolwincor_1win}} or \code{\link{rolwinmulcor_1win}}. 
}
 \item{typewidthwin}{``FULL'' is to estimate the windows from 2, 4, ..., to dim(inputdata)[1]) if \emph{Align} is equal to ``left'' or ``right'', or from 3, 5,..., to dim(inputdata)[1]) if Align is ``center''. The other option is ``PARTIAL'', please you should take into account that \emph{widthwin_1} and \emph{widthwin_1} MUST be ODD if the \emph{Align} option is ``center''. 
}
 \item{widthwin_1}{
First value for the size (length) of the windows when the option \emph{typewidthwin}=``PARTIAL'' is selected, the minimum value is 3 (the default value), but you must define this parameter (please note that \emph{widthwin_1} < \emph{widthwin_N}). 
}
 \item{widthwin_N}{
Last value for the size (length) of the windows when the option \emph{typewidthwin}=``PARTIAL'' is selected, by default is dim(inputdata)[1], but you must define this parameter (please note that \emph{widthwin_1} < \emph{widthwin_N}).
}
 \item{varX}{ 
Name of the  ``first'' or independent variable, e.g. ``X'' (please note that ``X'' is a vector of one element if KCASE=``BIVAR'' and a vector of several elements if KCASE=``MULVAR''. For the multi-variate case the names for ``X'' (the independent variables) will be defined as: varX=paste(``X1'', ``X2'',..., sep='', ''). 
}
 \item{varY}{
Name of the ``second'' (bi-variate case) or dependent variable (multi-variate case), e.g. ``Y''. 
}
 \item{rmltrd}{
Remove (by default is ``TRUE''; ``FALSE'' otherwise) the linear trend in the time series under analysis. 
}
 \item{Scale}{
Scale (by default is ``TRUE''; ``FALSE'' otherwise) is used to ``normalize'' or ``standardize'' the time series under analysis.
}
 \item{coltsX, coltsY}{
Colors to be used when the variables are plotted, for the bi-variate case by default are ``black'' for ``X'' and ``blue'' for ``Y'', but other colors can be used. For the multi-variate case, colors for the dependent (``Y'') and independent variables (``X'') MUST be provided (e.g. coltsX=c("red","blue",...), coltsY="black"). 
}
 \item{CEXLAB, CEXAXIS}{
These parameters are used to plot the sizes of the X-axis and Y-axis labels and X- and Y-axis, by default these parameters have values of 1.15 and 1.05, respectively, but it is possible to use other values. 
}
 \item{LWDtsX, LWDtsY}{
Line-widths for the first and the second variable when these are plotted, for the bi-variate case by default these have values of 1, but other values (widths) can be used. For the multi-variate case and for the independent variables the line-widths MUST be provided (e.g. LWDtsX = c(1,2,...)). 
}
 \item{NUMLABX}{
Number of labels for (all) the X's axis, by the default is 5, but it is possible to use other values. 
}
 \item{parcen}{
These parameters contain two values: the first one is to control the position of the title, by default it is 0.5, but you should try with other close values to obtain the title centered, e.g. 0.4 or 0.8 (please avoid to use large values); the second value is to define the spaces between the names of variables, by default is 25 spaces, but you could try other values to fit properly the names of variables in the title. We use ``mtext'' to produce the title (please loot at \R>?mtext for more information). 
}
}
\details{
The \code{\link{plot_heatmap}} function plots the heat maps for the correlation coefficients and their respective p-values (corrected or not corrected) for all possible window-lengths (i.e., from five to the number of elements in the time series under analysis) or for a band of window-lengths. \code{\link{plot_heatmap}} uses the outputs of the functions  \code{\link{rolwincor_heatmap}} (bi-variate case) and \code{\link{rolwinmulcor_heatmap}} (multi-variate case). 
}
\value{
Output: a heat map (via screen) of the correlation coefficients and their respective (corrected or not corrected) p-values.}
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
\donttest{
######################################################################
# Testing the bi-variate case (heat map). Example: synthetic data! 
test_fun2 <- rolwincor_heatmap(syntDATA, varX="X", varY="Y",
                       CorMethod="spearman", typewidthwin="PARTIAL",
                       widthwin_1=11, widthwin_N=101, Align="center")
# Plotting the bi-variate case (heat map). Example: synthetic data! 
plot_heatmap(syntDATA, test_fun2$matcor, test_fun2$pvalscor, test_fun2$left_win, 
             test_fun2$righ_win, test_fun2$Windows, KCASE="BIVAR", typewidthwin="PARTIAL", 
           varX="X", varY="Y", widthwin_1=11, widthwin_N=101)
######################################################################
# Testing the bi-variate case (heat map). Example: real-life ecological data
######################################################################
SST_PC1 <- rolwincor_heatmap(YX_ecological_data[,c(1,3,2)], varX="SST",
                 varY="PC1", CorMethod="spearman", typewidthwin="FULL", 
                 Align="center", pvalcorectmethod="BH")
# Plotting the bi-variate case (heat map). Example: real-life ecological data 
plot_heatmap(YX_ecological_data[,c(1,3,2)], SST_PC1$matcor, SST_PC1$pvalscor,
     SST_PC1$left_win, SST_PC1$righ_win, SST_PC1$Windows, KCASE="BIVAR",
     typewidthwin="FULL", varX="SST", varY="PC1", coltsX="red", CEXLAB=1.15,
     CEXAXIS=1.65, coltsY="black", LWDtsX=2, LWDtsY=2)
######################################################################
# Testing the multi-variate case (heat map). Example: real-life ecological data
######################################################################
SST_TSI_PC1 <- rolwinmulcor_heatmap(YX_ecological_data, typewidthwin="FULL",
                                    Align="center", pvalcorectmethod="BH")
# Plotting the multi-variate case (heat map). Example: real-life ecological data
plot_heatmap(YX_ecological_data, SST_TSI_PC1$matcor, SST_TSI_PC1$pvalscor,
      SST_TSI_PC1$left_win, SST_TSI_PC1$righ_win, Rwidthwin=SST_TSI_PC1$Windows,
      KCASE="MULVAR", typewidthwin="FULL", varY="PC1", varX=c("SST", "TSI"),
      coltsY="black", coltsX=c("red", "orange"), CEXLAB=1.15, CEXAXIS=1.65,
      LWDtsX=rep(2,2), LWDtsY=2, parcen=c(0.45,15))
}
}
\keyword{heat_map_plot}
\keyword{plot_heatmap_bi_and_multi_variate}


