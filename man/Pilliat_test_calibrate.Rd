% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/changepoint_test_Pilliat.R
\name{Pilliat_test_calibrate}
\alias{Pilliat_test_calibrate}
\title{Generates detection thresholds for the Pilliat algorithm for testing for a single change-point using Monte Carlo simulation}
\usage{
Pilliat_test_calibrate(
  n,
  p,
  N = 100,
  tol = 1/100,
  threshold_bj_const = 6,
  bonferroni = TRUE,
  rescale_variance = TRUE,
  fast = FALSE,
  debug = FALSE
)
}
\arguments{
\item{n}{Number of observations}

\item{p}{Number time series}

\item{N}{Number of Monte Carlo samples used}

\item{tol}{False error probability tolerance}

\item{threshold_bj_const}{Leading constant for \eqn{p_0} for the Berk-Jones statistic}

\item{bonferroni}{If \code{TRUE}, a Bonferroni correction applied and the detection thresholds for each statistic is chosen by simulating the leading constant in the theoretical detection thresholds}

\item{rescale_variance}{If \code{TRUE}, each row of the data is rescaled by a MAD estimate}

\item{fast}{If \code{FALSE}, a change-point test is applied to each candidate change-point position in each interval. If FALSE, only the mid-point of each interval is considered}

\item{debug}{If \code{TRUE}, diagnostic prints are provided during execution}
}
\value{
A list containing 
  \item{thresholds_partial}{vector of thresholds for the Partial Sum statistic (respectively for \eqn{t=1,2,4,\ldots,2^{\lfloor\log_2(p)\rfloor}} number of terms in the partial sum)}
  \item{threshold_dense}{threshold for the dense statistic}
  \item{thresholds_bj}{vector of thresholds for the Berk-Jones static (respectively for \eqn{x=1,2,\ldots,x_0})}
}
\description{
R wrapper for function choosing detection thresholds for the Dense, Partial sum and Berk-Jones statistics in the multiple change-point detection algorithm of \insertCite{pilliat_optimal_2022;textual}{HDCD} for single change-point testing using Monte Carlo simulation. When \code{Bonferroni==TRUE}, the detection thresholds are chosen by simulating the leading constant in the theoretical detection thresholds given in \insertCite{pilliat_optimal_2022;textual}{HDCD}, similarly as described in Appendix B in \insertCite{moen2023efficient;textual}{HDCD} for ESAC. When \code{Bonferroni==TRUE}, the thresholds for the Berk-Jones statistic are theoretical and not chosen by Monte Carlo simulation.
}
\examples{
library(HDCD)
n = 50
p = 50

set.seed(100)
thresholds_test_emp = Pilliat_test_calibrate(n,p, bonferroni=TRUE,N=100, tol=1/100)
set.seed(100)
thresholds_test_emp_without_bonferroni = Pilliat_test_calibrate(n,p, 
                                         bonferroni=FALSE,N=100, tol=1/100)
thresholds_test_emp # thresholds with bonferroni correction
thresholds_test_emp_without_bonferroni # thresholds without bonferroni correction

# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
Y = matrix(rnorm(n*p), ncol = n, nrow=p)

# Adding a single sparse change-point to X (and not Y):
X[1:5, 25:50] = X[1:5, 25:50] +2
resX = Pilliat_test(X, 
                    threshold_dense=thresholds_test_emp$threshold_dense, 
                    thresholds_bj = thresholds_test_emp$thresholds_bj, 
                    thresholds_partial = thresholds_test_emp$thresholds_partial
)
resX
resY = Pilliat_test(Y, 
                    threshold_dense=thresholds_test_emp$threshold_dense, 
                    thresholds_bj = thresholds_test_emp$thresholds_bj, 
                    thresholds_partial = thresholds_test_emp$thresholds_partial
)
resY
}
