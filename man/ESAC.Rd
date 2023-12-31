% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ESAC.R
\name{ESAC}
\alias{ESAC}
\title{Efficient Sparsity Adaptive Change-point estimator}
\usage{
ESAC(
  X,
  threshold_d = 1.5,
  threshold_s = 1,
  alpha = 1.5,
  K = 5,
  debug = FALSE,
  empirical = FALSE,
  tol = 0.001,
  N = 1000,
  thresholds = NULL,
  thresholds_test = NULL,
  threshold_d_test = threshold_d,
  threshold_s_test = threshold_s,
  fast = FALSE,
  rescale_variance = TRUE,
  trim = FALSE,
  NOT = TRUE,
  midpoint = FALSE
)
}
\arguments{
\item{X}{Matrix of observations, where each row contains a time series}

\item{threshold_d}{Leading constant for \eqn{\lambda(t) \propto r(t)} for \eqn{t= p}. Only relevant when \code{thresholds=NULL}}

\item{threshold_s}{Leading constant for \eqn{\lambda(t) \propto r(t)} for \eqn{t\leq \sqrt{p\log n}}.  Only relevant when \code{thresholds=NULL}}

\item{alpha}{Parameter for generating seeded intervals}

\item{K}{Parameter for generating seeded intervals}

\item{debug}{If \code{TRUE}, diagnostic prints are provided during execution}

\item{empirical}{If \code{TRUE}, detection thresholds are based on Monte Carlo simulation using \code{\link{ESAC_calibrate}}}

\item{tol}{If \code{empirical=TRUE}, \code{tol} is the false error probability tolerance}

\item{N}{If \code{empirical=TRUE}, \code{N} is the number of Monte Carlo samples used}

\item{thresholds}{Vector of manually chosen values of \eqn{\lambda(t)} for \eqn{t \in \mathcal{T}}, decreasing in \eqn{t}}

\item{thresholds_test}{Vector of manually chosen values of \eqn{\gamma(t)} for \eqn{t \in \mathcal{T}}, decreasing in \eqn{t}}

\item{threshold_d_test}{Leading constant for \eqn{\gamma(t) \propto r(t)} for \eqn{t=p}.  Only relevant when \code{empirical=FALSE} and \code{thresholds_test=NULL}}

\item{threshold_s_test}{Leading constant for \eqn{\gamma(t) \propto r(t)} for \eqn{t\leq \sqrt{p\log n}}. Only relevant when \code{empirical=FALSE} and \code{thresholds_test=NULL}}

\item{fast}{If \code{TRUE}, ESAC only tests for a change-point at the midpoint of each seeded interval}

\item{rescale_variance}{If \code{TRUE}, each row of the data is re-scaled by a MAD estimate using \code{\link{rescale_variance}}}

\item{trim}{If \code{TRUE}, interval trimming is performed}

\item{NOT}{If \code{TRUE}, ESAC uses Narrowest-Over-Threshold selection of change-points}

\item{midpoint}{If \code{TRUE}, change-point positions are estimated by the mid-point of the seeded interval in which the penalized score is the largest}
}
\value{
A list containing 
  \item{changepoints}{vector of estimated change-points}
  \item{changepointnumber}{number of changepoints}
  \item{CUSUMval}{the penalized score at the corresponding change-point in \code{changepoints}}
  \item{coordinates}{a matrix of zeros and ones indicating which time series are affected by a change in mean, with each row corresponding to the change-point in \code{changepoints}}
  \item{scales}{vector of estimated noise level for each series}
  \item{startpoints}{start point of the seeded interval detecting the corresponding change-point in \code{changepoints}}
  \item{endpoints}{end point of the seeded interval detecting the corresponding change-point in \code{changepoints}}
  \item{thresholds}{vector of values of \eqn{\lambda(t)} for \eqn{t \in \mathcal{T}} in decreasing order}
  \item{thresholds_test}{vector of values of \eqn{\gamma(t)} for \eqn{t \in \mathcal{T}} in decreasing order}
}
\description{
R wrapper for C function implementing the full ESAC algorithm \insertCite{@see @moen2023efficient}{HDCD}.
}
\examples{
library(HDCD)
n = 50
p = 50
set.seed(100)
# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# Adding a single sparse change-point:
X[1:5, 26:n] = X[1:5, 26:n] +1

# Vanilla ESAC:
res = ESAC(X)
res$changepoints

# Manually setting leading constants for \lambda(t) and \gamma(t)
res = ESAC(X, 
           threshold_d = 2, threshold_s = 2, #leading constants for \lambda(t)
           threshold_d_test = 2, threshold_s_test = 2 #leading constants for \gamma(t)
)
res$changepoints #estimated change-point locations

# Empirical choice of thresholds:
res = ESAC(X, empirical = TRUE, N = 100, tol = 1/100)
res$changepoints

# Manual empirical choice of thresholds (equivalent to the above)
thresholds_emp = ESAC_calibrate(n,p, N=100, tol=1/100)
res = ESAC(X, thresholds_test = thresholds_emp[[1]])
res$changepoints
}
\references{
\insertAllCited{}
}
