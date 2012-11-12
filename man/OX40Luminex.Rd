\name{OX40Luminex}
\docType{data}
\alias{OX40Luminex}
\alias{luminex}
\alias{lum}
\alias{print.luminex}

\usage{
  data(OX40Luminex)
}

\title{A sample data of beads based microRNA microarray}
\description{
A sample data of Luminex beads based microRNA microarray.
}
\format{A data frame with 3640 rows and 8 columns (variables).
  \tabular{lll}{
    \code{Gene}  \tab character  \tab MicroRNA names. \cr
    \code{Treat0} \tab numeric \tab measures under treatment with code '0'\cr
    \code{Treat1} \tab numeric \tab measures under treatment with code '1'\cr
    \code{Treat2} \tab numeric \tab measures under treatment with code '2'\cr
    \code{Treat3} \tab numeric \tab measures under treatment with code '3'\cr
    \code{Background} \tab numeric \tab measures of background noises\cr
    \code{Sample}  \tab numeric  \tab sample id. \cr
    \code{Pool}   \tab numeric  \tab Pool id. \cr
  }

 An object \code{lum} of class \code{luminx}.

}

\references{

 Wang, B., Zhang, S-G., Wang, X-F., Tan, M. and Xi, Y. (2012) "Testing
 for differentially-expressed microRNAs with errors-in-variables
 nonparametric regression", PLoS ONE 7(5):
 e37537. doi:10.1371/journal.pone.0037537.

 Wang, B., Wang, X-F. and Xi, Y. (2011) "Normalizing bead-based
 microRNA expression data: a measurement error model-based approach",
 Bioinformatics, 27(11), 1506-1512.

 Wang, B., Howell, P., Bruheim, S. Ju, J, Owen, L.B., Fodstad, O. and
 Xi, Y. (2011). "Systematic Evaluation of Three microRNA Profiling
 Platforms: Microarray, Beads Array, and Quantitative Real-Time PCR
 Array", PLoS ONE 6(2): e17167.

 Wang, B. and Wertelecki, W. (2012) "Density Estimation for Data With
 Rounding Errors", Computational Statistics and Data Analysis, (in
 press). doi: 10.1016/j.csda.2012.02.016. link

 Wang, B., Wang, X-F. and Xi, Y. (2011) "Normalizing bead-based
 microRNA expression data: a measurement error model-based approach",
 Bioinformatics, 27(11), 1506-1512.

 Wang, X-F. and Wang, B. (2011) "Deconvolution Estimation in
 Measurement Error Models: The R Package decon", Journal of Statistical
 Software, 39(10), 1-24.

}

\examples{
  data(OX40Luminex)
}

\keyword{datasets}
