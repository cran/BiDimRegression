\name{BiDimRegression-package}
\alias{BiDimRegression-package}
\docType{package}
\title{package to calculate the bidimensional regression between two 2D configurations}
\description{An S3 class with a method for calculates the bidimensional regression between two 2D configurations

following the approach by Tobler (1965). }

\details{

\tabular{ll}{
Package: \tab BiDimRegression\cr
Type: \tab Package\cr
Version: \tab 1.0-6\cr
Date: \tab 2014-03-22\cr
License: \tab What license is it under?\cr
}
~~ Overview of the most important functions ~~
BiDimRegression <- function (coord) \\
summary.BiDimRegression <- function(object, ...) \\
print.summary.BiDimRegression <- function(x, ...) \\
print.BiDimRegression <- function(x, ...)\\}

\author{Who wrote it: Claus-Christian Carbon\\
Maintainer: Claus-Christian Carbon <ccc@experimental-psychology.com>}

\references{*** The routine is described in detail in:\\
Carbon, C. C. (2013). BiDimRegression: Bidimensional Regression Modeling Using R. Journal of Statistical Software,\\
Code Snippets, 52(1), 1-11 (URL http://www.jstatsoft.org/v52/c01/)\\
\\
*** The routine is based on the following literature:\\
1) Tobler, W. [R.](1965). Computation of the corresponding of geographical patterns. Papers of the Regional Science Association, 15, 131-139.\\
2) Tobler, W. R. (1966). Medieval distortions: Projections of ancient maps. Annals of the Association of American Geographers, 56(2), 351-360.\\
3) Tobler, W. R. (1994). Bidimensional regression. Geographical Analysis, 26(3), 187-212.\\
4) Friedman, A., & Kohler, B. (2003). Bidimensional regression: Assessing the configural similarity and accuracy of cognitive maps and other two-dimensional data sets. Psychological Methods, 8(4), 468-491.\\
5) Nakaya, T. (1997). Statistical inferences in bidimensional regression models. Geographical Analysis, 29(2), 169-186.\\
6) Waterman, S., & Gordon, D. (1984). A quantitative-comparative approach to analysis of distortion in mental maps. Professional Geographer, 36(3),\\ 326-337.}

\keyword{ regression }
\keyword{ bidimensional }
\keyword{ Pearson }
\keyword{ correlation }
\keyword{ AIC }
\keyword{ similarity }
\keyword{ comparison }

\examples{
data(NakayaData)
resultsBiDimRegr <- BiDimRegression(NakayaData)
print(resultsBiDimRegr)
summary(resultsBiDimRegr)
}
