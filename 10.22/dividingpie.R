library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(reshape2)
library(corrplot)
source("KmeansGap.r")

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

set.seed(1)

MacArthur = 
  \(time, state, parms){
    N = state[1:params$nspec]
    R = state[(params$nspec + 1):(params$nspec + params$nres)]
    dNdt = with(parms, alpha + N * ((C %*% R) - m))
    dRdt = with(parms, R * (r * (1 - R / K) - t(C) %*% N))
    return(list(c(dNdt, dRdt)))
  }

circ_dist = function(vec1, vec2) {
  y = abs(outer(vec1, vec2, '-'))
  idx = which(y > 0.5)
  y[idx] = 1 - y[idx]
  return(y)
}

nspec = 120
nres = 8

res_trait_1 = seq(0, (nres - 1) / nres, l = nres)
spec_trait_1 = seq(0, 1, l = nspec)
dists_1 = circ_dist(spec_trait_1, res_trait_1)
w_1 = .7
C = exp(-w_1 * dists_1^2)

res_trait_2 = sample(res_trait_1, size = nres, replace = F)
spec_trait_2 = sample(spec_trait_1, size = length(spec_trait_1), replace = F)
dists_2 = circ_dist(spec_trait_2, res_trait_2)
w_2 = .3
C_2ax = exp(-(w_1 * dists_1^2 + w_2 * dists_2^2))

cor_matrix_2ax = cor(C_2ax)
corrplot(cor_matrix_2ax, type = "upper",  
         method = "square",  
         addCoef.col = "black",  
         tl.col = "black", tl.srt = 45, title = "2 axes, sigmas equal")

res_trait_3 = sample(res_trait_1, size = nres, replace = F)
spec_trait_3 = sample(spec_trait_1, size = length(spec_trait_1), replace = F)
dists_3 = circ_dist(spec_trait_3, res_trait_3)
w_2 = .15
w_3 = .15
C_3ax = exp(-(w_1 * dists_1^2 + w_2 * dists_2^2 + w_3 * dists_3^2))

cor_matrix_3ax = cor(C_3ax)
corrplot(cor_matrix_3ax, type = "upper",  
         method = "square",  
         addCoef.col = "black",  
         tl.col = "black", tl.srt = 45, title = "3 axes, all sigmas equal")