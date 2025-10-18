
# infotest

<!-- badges: start -->
<!-- badges: end -->

The goal of infotest is to implement the Cameron & Trivedi (1990) Information Matrix test for regression models, decomposing it into tests of heteroskedasticity, skewness, and kurtosis.
## Installation

You can install the development version of infotest like so:

``` r
install.packages("pak")
pak::pkg_install("guliyevh/infotest")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(infotest)

# Fit a simple linear model
model <- lm(Sepal.Length ~ Sepal.Width + Petal.Length, data = iris)

# Run the Cameron & Trivedi (1990) IM test decomposition
infotest(model)

White's (1980) heteroskedasticity test
H0: Homoskedasticity
Ha: Unrestricted heteroskedasticity

chi2(5) = 16.3799
Prob > chi2 = 0.0058

Cameron & Trivedi's (1990) decomposition of IM-test

----------------------------------------
Source                     chi2     df          p
----------------------------------------
Heteroskedasticity        16.38      5     0.0058
Skewness                   0.15      2     0.9273
Kurtosis                   0.49      1     0.4834
Joint Test                17.02      8     0.0299
----------------------------------------


# Without White's (1980) heteroskedasticity test
infotest(model, white = FALSE)

Cameron & Trivedi's (1990) decomposition of IM-test

----------------------------------------
Source                     chi2     df          p
----------------------------------------
Heteroskedasticity        16.38      5     0.0058
Skewness                   0.15      2     0.9273
Kurtosis                   0.49      1     0.4834
Joint Test                17.02      8     0.0299
----------------------------------------

```
## References

Cameron, A. C., & Trivedi, P. K. (1990). The information matrix test and its applied alternatives. Econometric Theory, 6(2), 179–195.

White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. Econometrica, 48(4), 817–838.


## Citation

If you use this package, please cite:

Hasraddin Guliyev (2025). infotest: Information Matrix Test for Regression Models. R package version 0.1.0.
