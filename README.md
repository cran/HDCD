
<!-- README.md is generated from README.Rmd. Please edit that file -->

# High-Dimensional Change-point Detection

<!-- badges: start -->
<!-- badges: end -->

HDCD contains efficient implementations of several multiple change-point
detection algorithms, including Efficient Sparsity Adaptive Change-point
estimator (ESAC) and Informative sparse projection for estimating
change-points (Inspect).

## Installation

You can install the development version of HDCD from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("peraugustmoen/HDCD")
```

## Example

This is a basic example which shows you how to run ESAC:

``` r
library(HDCD)
n = 50
p = 50
set.seed(100)
# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# Adding a single sparse change-point (at location \eta = 25):
X[1:5, 26:n] = X[1:5, 26:n] +2

# Vanilla ESAC:
res = ESAC(X)
res$changepoints
#> [1] 25

# Manually setting leading constants for \lambda(t) and \gamma(t)
res = ESAC(X,
           threshold_d = 2, threshold_s = 2, #leading constants for \lambda(t)
           threshold_d_test = 2, threshold_s_test = 2 #leading constants for \gamma(t)
)
res$changepoints #estimated change-point locations
#> [1] 25

# Empirical choice of thresholds:
res = ESAC(X, empirical = TRUE, N = 100, tol = 1/100)
res$changepoints
#> [1] 25


# Manual empirical choice of thresholds (equivalent to the above)
thresholds_emp = ESAC_calibrate(n,p, N=100, tol=1/100)
res = ESAC(X, thresholds_test = thresholds_emp[[1]])
res$changepoints
#> [1] 25
```

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
