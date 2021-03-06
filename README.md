# agnostic

The goal of agnostic is to perform agnostic hypothesis tests in R. In order to be better user friendly, this package is designed to work similarly to R base.

## Installation

You can install agnostic directly into R from github with devtools:

``` r
if(!("devtools" %in% rownames(installed.packages())))
  install.packages('devtools')
devtools::install_github("vcoscrato/agnostic")
```
## Example

This example performs univariate and multivariate analysis of variance (ANOVA / MANOVA) under the agnostic perspective.

``` r
library(agnostic)

# Test data
obs <- c(sample(20:30, size = 24, replace = TRUE))
obs2 <- rnorm(24)
groups <- as.factor(rep(c("A", "B", "c", "D"), 6))

# Build the model (ANOVA)
mod1 <- agnostic.aov(obs ~ groups)

# Get the output under the agnostic perspective
summary(mod1)

# Build the model (MANOVA)
mod2 <- agnostic.aov(cbind(obs, obs2) ~ groups)

# Get the output under the agnostic perspective
summary(mod2)
```
