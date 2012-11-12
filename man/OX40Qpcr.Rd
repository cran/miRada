\name{OX40Qpcr}
\docType{data}
\alias{qrtpcr}
\alias{qpcr}
\alias{print.qrtpcr}

\usage{
  data(OX40Qpcr)
}

\title{qRT-PCR results}
\description{
The qRT-PCR microRNA microarray profiles for 10 samples under four different treatments.
}
\format{A data frame with 27000 rows and 6 columns (variables).
  \tabular{lll}{
    \code{Sample}  \tab numeric  \tab sample id. \cr
    \code{Treat}   \tab numeric  \tab codes of treatments: 0,1,2,3. \cr
    \code{Gene}  \tab character  \tab MicroRNA names. \cr
    \code{dCT}   \tab numeric  \tab delta Ct values. \cr
    \code{ddCt}   \tab numeric  \tab delta delta Ct values. \cr
    \code{RQ}   \tab numeric  \tab \code{RQ = 2^(-ddCt)}. \cr
  }

 An object \code{qpcr} of class \code{qrtpcr}.

}

\references{
Wang, B., Wang, X-F. and Xi, Y. (2011) "Normalizing bead-based microRNA expression data: a measurement error model-based approach", Bioinformatics, 27(11), 1506-1512.

Wang, B., Howell, P., Bruheim, S. Ju, J, Owen, L.B., Fodstad, O. and Xi, Y. (2011). Systematic Evaluation of Three microRNA Profiling Platforms: Microarray, Beads Array, and Quantitative Real-Time PCR Array, PLoS ONE 6(2): e17167.
}

\examples{

  data(OX40Qpcr)

}

\keyword{datasets}
