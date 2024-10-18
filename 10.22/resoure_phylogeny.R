library(ggplot2)
library(deSolve)
library(tibble)
library(tidyr)
library(dplyr)
library(reshape2)
library(purrr)

set.seed(1)

nspec = 10
nres = 10

P = matrix(0, nrow = nres, ncol = nres)
diag(P) = 1
P[lower.tri(P)] = runif(((nres)^2 - nres) / 2 , min = 0, max = 1)
P[upper.tri(P)] = t(P[lower.tri(P)])