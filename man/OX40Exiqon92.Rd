\name{OX40Exiqon}
\docType{data}
\alias{OX40}
\alias{OX40Exiqon}
\alias{ox40ex92}
\alias{exiqon}
\alias{print.exiqon}
\alias{summary.exiqon}

\title{A sample data of Exiqon microRNA microarray}
\description{
 
 A sample microRNA expression profiling data using miRCURY LNA
 microRNA Array based on miRbase 9.2 (Exiqon Inc., Denmark) for 40
 human osteosarcoma xenograft specimens, including 10 samples for each
 chemotherapeutic treatment (Cisplatin, Doxorubicin, and Ifosfamide)
 plus 10 non-treated samples.

}

\usage{
  data(OX40Exiqon)
}

\format{A data frame with 168192 rows and 7 columns (variables).
  \tabular{lll}{
    \code{Gene}  \tab character  \tab MicroRNA names. \cr
    \code{Flag} \tab numeric \tab an indictor showing the probe signal quality.\cr
    \code{Repeat} \tab numeric \tab an indictor showing the technical replicates.\cr
    \code{Treat} \tab numeric \tab treatment type: coded.\cr
    \code{Sample}  \tab numeric  \tab sample id. \cr
    \code{Signal} \tab numeric \tab fluorescence intensity measures\cr
    \code{Background} \tab numeric \tab measures of background noises\cr
  }
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

 data(OX40Exiqon)
 ex = structure(list(x = ox40ex92), class='exiqon')
 print(ex)

}

\keyword{datasets}
