\name{BiDimRegression}
\alias{BiDimRegression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~
function to calculate the bidimensional regression between two 2D configurations
}
\description{
%%  ~~ 
An S3 class with a method for calculates the bidimensional regression between two 2D configurations

following the approach by Tobler (1965). ~~
}
\usage{
BiDimRegression(coord)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coord}{
%%     ~~Describe \code{coord} here~~
}
}
\references{
%% ~put references to the literature/web site here ~
   Tobler, W. [R.](1965). Computation of the corresponding of geographical patterns. Papers of the Regional Science Association, 15, 131-139.
   Tobler, W. R. (1966). Medieval distortions: Projections of ancient maps. Annals of the Association of American Geographers, 56(2), 351-360.
   Tobler, W. R. (1994). Bidimensional regression. Geographical Analysis, 26(3), 187-212.
   Friedman, A., & Kohler, B. (2003). Bidimensional regression: Assessing the configural similarity and accuracy of cognitive maps and other two-dimensional data sets. Psychological Methods, 8(4), 468-491.
   Nakaya, T. (1997). Statistical inferences in bidimensional regression models. Geographical Analysis, 29(2), 169-186.
   Waterman, S., & Gordon, D. (1984). A quantitative-comparative approach to analysis of distortion in mental maps. Professional Geographer, 36(3), 326-337.
}
\author{
%%  ~~
Claus-Christian Carbon
}
\note{
%%  ~~further notes~~
The author is grateful to Waldo R. Tobler, now Professor Emeritus at the Department of Geography, University of California, Santa Barbara, for providing his original publications and his helpful correspondence and to Dan Montello for calling the author's attention to Tobler's work many years ago. I would like to thank Alinda Friedman, Gregory Francis, Jan de Leeuw, Achim Zeileis, two anonymous reviewers and Arne Terkowski and Chris Nappert for valuable comments on an earlier version of this paper, and Andrea Lyman and Vera M. Hesslinger for proofreading the manuscript. Last but not least, I am very indepted to Tomoki Nakaya, who has developed the original inference statistics of the overall models and the referring parameters and who helped me with reanalyzing these statistics as well as with ensuring the reliability of the used methods. Thank you!
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (coord) 
{
    timerStart <- proc.time()
    timerStart <- proc.time()
    n <- dim(coord)[1]
    vecZero <- c(rep(0, n))
    vecOne <- c(rep(1, n))
    A <- coord$depV1
    B <- coord$depV2
    X <- coord$indepV1
    Y <- coord$indepV2
    Am <- mean(A)
    Bm <- mean(B)
    Xm <- mean(X)
    Ym <- mean(Y)
    X2 <- sum(X^2)
    Y2 <- sum(Y^2)
    sumX2Y2 <- sum(X^2 + Y^2)
    A2 <- sum(A^2)
    B2 <- sum(B^2)
    sumA2B2 <- sum(A^2 + B^2)
    varA <- (sum((A - Am) * (A - Am)))/(n)
    varB <- (sum((B - Bm) * (B - Bm)))/(n)
    varX <- (sum((X - Xm) * (X - Xm)))/(n)
    varY <- (sum((Y - Ym) * (Y - Ym)))/(n)
    covABXY <- sum(((A - Am) + (B - Bm)) * ((X - Xm) + (Y - Ym)))/n
    euc_par <- data.frame(ps = 4L, name = "Euclidean")
    euc_dataMatrix <- matrix(c(vecOne, vecZero, vecZero, vecOne, 
        X, Y, -Y, X), ncol = 4)
    euc_target <- matrix(c(A, B), ncol = 1)
    euc_data <- data.frame(y = (euc_target), x1 = (euc_dataMatrix[, 
        1]), x2 = (euc_dataMatrix[, 2]), x3 = (euc_dataMatrix[, 
        3]), x4 = (euc_dataMatrix[, 4]))
    euc_regression <- lm(y ~ 0 + x1 + x2 + x3 + x4, data = euc_data)
    euc_alpha <- data.frame(coeff = c(NA, NA), SE = c(NA, NA), 
        tValue = c(NA, NA), pValue = c(NA, NA))
    euc_beta <- data.frame(coeff = c(NA, NA), SE = c(NA, NA), 
        tValue = c(NA, NA), pValue = c(NA, NA))
    for (iAlpha in 1:2) {
        for (iPar in 1:4) {
            euc_alpha[iAlpha, iPar] <- summary(euc_regression)$coeff[iAlpha, 
                iPar]
        }
    }
    for (iBeta in 1:2) {
        for (iPar in 1:4) {
            euc_beta[iBeta, iPar] <- summary(euc_regression)$coeff[iBeta + 
                2, iPar]
        }
    }
    euc_scaleFactorX <- sqrt(euc_beta$coeff[1] * euc_beta$coeff[1] + 
        euc_beta$coeff[2] * euc_beta$coeff[2])
    euc_scaleFactorY <- euc_scaleFactorX
    euc_angleRAD <- atan(euc_beta$coeff[2]/euc_beta$coeff[1])
    euc_angleDEG <- euc_angleRAD * 180/pi
    euc_shear = 0L
    euc_Apred <- euc_alpha$coeff[1] + euc_beta$coeff[1] * X - 
        euc_beta$coeff[2] * Y
    euc_Bpred <- euc_alpha$coeff[2] + euc_beta$coeff[2] * X + 
        euc_beta$coeff[1] * Y
    euc_Xpred <- (euc_alpha$coeff[2] * euc_beta$coeff[2] - B * 
        euc_beta$coeff[2] - euc_alpha$coeff[1] * euc_beta$coeff[1] + 
        A * euc_beta$coeff[1])/(euc_beta$coeff[1] * euc_beta$coeff[1] - 
        euc_beta$coeff[2] * euc_beta$coeff[2])
    euc_Ypred <- (euc_beta$coeff[1] * B - euc_alpha$coeff[2] * 
        euc_beta$coeff[1] + euc_alpha$coeff[1] * euc_beta$coeff[2] - 
        A * euc_beta$coeff[2])/(euc_beta$coeff[1] * euc_beta$coeff[1] + 
        euc_beta$coeff[2] * euc_beta$coeff[2])
    euc_r <- sqrt(sum(((euc_Apred - Am) * (euc_Apred - Am)) + 
        ((euc_Bpred - Bm) * (euc_Bpred - Bm)))/sum(((A - Am) * 
        (A - Am)) + ((B - Bm) * (B - Bm))))
    euc_rsqr <- euc_r * euc_r
    euc_F <- ((2 * n - euc_par$ps)/2) * (euc_rsqr/(1 - euc_rsqr))
    euc_df1 <- 2L
    euc_df2 <- 2 * n - euc_par$ps
    euc_p <- pf(euc_F, euc_df1, euc_df2, lower.tail = FALSE, 
        log.p = FALSE)
    euc_dDistanceXY <- sqrt(sum((X - euc_Xpred) * (X - euc_Xpred)) + 
        sum((Y - euc_Ypred) * (Y - euc_Ypred)))
    euc_dDistanceXYSqr <- euc_dDistanceXY * euc_dDistanceXY
    euc_dMaxXYSqr <- sum((X - Xm) * (X - Xm) + ((Y - Ym) * (Y - 
        Ym)))
    euc_dMaxXY <- sqrt(euc_dMaxXYSqr)
    euc_dDistanceAB <- sqrt(sum((A - euc_Apred) * (A - euc_Apred)) + 
        sum((B - euc_Bpred) * (B - euc_Bpred)))
    euc_dDistanceABSqr <- euc_dDistanceAB * euc_dDistanceAB
    euc_dMaxABSqr <- sum((A - Am) * (A - Am) + ((B - Bm) * (B - 
        Bm)))
    euc_dMaxAB <- sqrt(euc_dMaxABSqr)
    euc_diABSqr <- euc_dDistanceABSqr/euc_dMaxABSqr
    euc_diAB <- sqrt(euc_diABSqr)
    euc_diXYSqr <- euc_dDistanceABSqr/euc_dMaxXYSqr
    euc_diXY <- sqrt(euc_diXYSqr)
    euc_dAICso <- 2L * n * log(1 - euc_rsqr) + 2L * (euc_par$ps - 
        2L)
    aff_par <- data.frame(ps = 6L, name = "Affine")
    aff_dataMatrix <- matrix(c(vecOne, vecZero, vecZero, vecOne, 
        X, vecZero, Y, vecZero, vecZero, X, vecZero, Y), ncol = 6)
    aff_target <- matrix(c(A, B), ncol = 1)
    aff_data <- data.frame(y = (aff_target), x1 = (aff_dataMatrix[, 
        1]), x2 = (aff_dataMatrix[, 2]), x3 = (aff_dataMatrix[, 
        3]), x4 = (aff_dataMatrix[, 4]), x5 = (aff_dataMatrix[, 
        5]), x6 = (aff_dataMatrix[, 6]))
    aff_regression <- lm(y ~ 0 + x1 + x2 + x3 + x4 + x5 + x6, 
        data = aff_data)
    aff_alpha <- data.frame(coeff = c(NA, NA), SE = c(NA, NA), 
        tValue = c(NA, NA), pValue = c(NA, NA))
    aff_beta <- data.frame(coeff = c(NA, NA), SE = c(NA, NA), 
        tValue = c(NA, NA), pValue = c(NA, NA))
    for (iAlpha in 1:2) {
        for (iPar in 1:4) {
            aff_alpha[iAlpha, iPar] <- summary(aff_regression)$coeff[iAlpha, 
                iPar]
        }
    }
    for (iBeta in 1:4) {
        for (iPar in 1:4) {
            aff_beta[iBeta, iPar] <- summary(aff_regression)$coeff[iBeta + 
                2, iPar]
        }
    }
    aff_angleRAD <- atan(aff_beta$coeff[3]/aff_beta$coeff[1])
    aff_angleDEG <- aff_angleRAD * 180/pi
    if (aff_beta$coeff[1] < 0) {
        aff_angleDEG <- aff_angleDEG + 180
    }
    aff_shear <- (((aff_beta$coeff[4]/aff_beta$coeff[2]) * sin(aff_angleRAD)) + 
        cos(aff_angleRAD))/(((aff_beta$coeff[4]/aff_beta$coeff[2]) * 
        cos(aff_angleRAD)) - sin(aff_angleRAD))
    aff_scaleFactorX <- sqrt(aff_beta$coeff[1] * aff_beta$coeff[1] + 
        aff_beta$coeff[3] * aff_beta$coeff[3])
    if (is.nan(aff_shear)) {
        aff_shear <- (aff_beta$coeff[1] - cos(aff_angleRAD) * 
            aff_scaleFactorX)/aff_beta$coeff[3]
    }
    if (is.nan(aff_shear)) {
        aff_shear <- (sin(aff_angleRAD) * aff_scaleFactorX + 
            aff_beta$coeff[2])/aff_beta$coeff[4]
    }
    aff_scaleFactorY <- aff_beta$coeff[2]/(aff_shear * cos(aff_angleRAD) - 
        sin(aff_angleRAD))
    if (is.nan(aff_scaleFactorY)) {
        aff_scaleFactorY <- aff_scaleFactorX
    }
    aff_Apred <- aff_alpha$coeff[1] + aff_beta$coeff[1] * X + 
        aff_beta$coeff[2] * Y
    aff_Bpred <- aff_alpha$coeff[2] + aff_beta$coeff[3] * X + 
        aff_beta$coeff[4] * Y
    aff_Xpred <- -(aff_Bpred * aff_beta$coeff[2] - aff_Apred * 
        aff_beta$coeff[4] - aff_alpha$coeff[2] * aff_beta$coeff[2] + 
        aff_alpha$coeff[1] * aff_beta$coeff[4])/(aff_beta$coeff[2] * 
        aff_beta$coeff[3] + aff_beta$coeff[1] * aff_beta$coeff[4])
    aff_Ypred <- -(aff_Apred * aff_beta$coeff[3] - aff_Bpred * 
        aff_beta$coeff[1] - aff_alpha$coeff[1] * aff_beta$coeff[3] + 
        aff_alpha$coeff[2] * aff_beta$coeff[1])/(aff_beta$coeff[2] * 
        aff_beta$coeff[3] + aff_beta$coeff[1] * aff_beta$coeff[4])
    aff_r <- sqrt(sum(((aff_Apred - Am) * (aff_Apred - Am)) + 
        ((aff_Bpred - Bm) * (aff_Bpred - Bm)))/sum(((A - Am) * 
        (A - Am)) + ((B - Bm) * (B - Bm))))
    aff_rsqr <- aff_r * aff_r
    aff_F <- ((2 * n - aff_par$ps)/4) * (aff_rsqr/(1 - aff_rsqr))
    aff_df1 <- 4L
    aff_df2 <- 2 * n - aff_par$ps
    aff_p <- pf(aff_F, aff_df1, aff_df2, lower.tail = FALSE, 
        log.p = FALSE)
    aff_dDistanceXY <- sqrt(sum((X - aff_Xpred) * (X - aff_Xpred)) + 
        sum((Y - aff_Ypred) * (Y - aff_Ypred)))
    aff_dDistanceXYSqr <- aff_dDistanceXY * aff_dDistanceXY
    aff_dMaxXYSqr <- sum((X - Xm) * (X - Xm) + ((Y - Ym) * (Y - 
        Ym)))
    aff_dMaxXY <- sqrt(aff_dMaxXYSqr)
    aff_dDistanceAB <- sqrt(sum((A - aff_Apred) * (A - aff_Apred)) + 
        sum((B - aff_Bpred) * (B - aff_Bpred)))
    aff_dDistanceABSqr <- aff_dDistanceAB * aff_dDistanceAB
    aff_dMaxABSqr <- sum((A - Am) * (A - Am) + ((B - Bm) * (B - 
        Bm)))
    aff_dMaxAB <- sqrt(aff_dMaxABSqr)
    aff_diABSqr <- aff_dDistanceABSqr/aff_dMaxABSqr
    aff_diAB <- sqrt(aff_diABSqr)
    aff_diXYSqr <- aff_dDistanceABSqr/aff_dMaxXYSqr
    aff_diXY <- sqrt(aff_diXYSqr)
    aff_dAICso <- 2 * n * log(1 - aff_rsqr) + 2 * (aff_par$ps - 
        2)
    dAICea <- 2 * n * log((1 - aff_rsqr)/(1 - euc_rsqr)) + 2 * 
        (aff_par$ps - euc_par$ps)
    f_ea <- ((2 * n - aff_par$ps)/(aff_par$ps - euc_par$ps)) * 
        ((aff_rsqr - euc_rsqr)/(1 - aff_rsqr))
    df1_ea <- as.integer(aff_par$ps - euc_par$ps)
    df2_ea <- as.integer(2 * n - aff_par$ps)
    p_ea <- pf(f_ea, aff_df1, aff_df2, lower.tail = FALSE, log.p = FALSE)
    tTestDF <- aff_regression$df[1]
    timerStop <- proc.time()
    timeDeltaMS <- timerStop[3] - timerStart[3]
    res_euc <- data.frame(r = euc_r, rsqr = euc_rsqr, diABSqr = euc_diABSqr, 
        dMaxABSqr = euc_dMaxABSqr, diXYSqr = euc_diXYSqr, dMaxXYSqr = euc_dMaxXYSqr, 
        scaleFactorX = euc_scaleFactorX, scaleFactorY = euc_scaleFactorY, 
        angleDEG = euc_angleDEG, shear = euc_shear, ttestDF = tTestDF, 
        alpha1 = euc_alpha[1, ], alpha2 = euc_alpha[2, ], beta1 = euc_beta[1, 
            ], beta2 = euc_beta[2, ], beta3 = euc_beta[3, ], 
        beta4 = euc_beta[4, ], fValue = euc_F, df1 = euc_df1, 
        df2 = euc_df2, pValue = euc_p, dAICso = euc_dAICso)
    res_aff <- data.frame(r = aff_r, rsqr = aff_rsqr, diABSqr = aff_diABSqr, 
        dMaxABSqr = aff_dMaxABSqr, diXYSqr = aff_diXYSqr, dMaxXYSqr = aff_dMaxXYSqr, 
        scaleFactorX = aff_scaleFactorX, scaleFactorY = aff_scaleFactorY, 
        angleDEG = aff_angleDEG, shear = aff_shear, ttestDF = tTestDF, 
        alpha1 = aff_alpha[1, ], alpha2 = aff_alpha[2, ], beta1 = aff_beta[1, 
            ], beta2 = aff_beta[2, ], beta3 = aff_beta[3, ], 
        beta4 = aff_beta[4, ], fValue = aff_F, df1 = aff_df1, 
        df2 = aff_df2, pValue = aff_p, dAICso = aff_dAICso)
    euclideanVSaffine <- data.frame(dAIC = dAICea, fValue = f_ea, 
        pValue = p_ea, df1 = df1_ea, df2 = df2_ea)
    results_sum <- data.frame(euclidean = res_euc, affine = res_aff, 
        eucVSaff = euclideanVSaffine)
    class(results_sum) <- "BiDimRegression"
    invisible(results_sum)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ bidimensional }% __ONLY ONE__ keyword per line
\keyword{ Pearson }% __ONLY ONE__ keyword per line
\keyword{ correlation }% __ONLY ONE__ keyword per line
\keyword{ AIC }% __ONLY ONE__ keyword per line


