library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(reshape2)
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

circ_dist = \(vec1, vec2) {
  y = abs(outer(vec1, vec2, '-'))
  idx = which(y > 0.5)
  y[idx] = 1 - y[idx]
  return(y)
}

makeC = \(nspec, nres, ndims, weights) {
  res_traits = seq(0, (nres - 1) / nres, l = nres)
  spec_traits = matrix(data = runif(nspec * ndims, min = 0, max = 1), nrow = ndims, ncol = nspec)
  dists_list = list()
  for (i in 1:ndims) {
    dists_list[[i]] = circ_dist(spec_traits[i,], res_traits)
  }
  
  weights = rep(1 / ndims, ndims)
  
  weighted_sqr_dists_sum = 0
  for (i in 1:ndims) {
    weighted_sqr_dists_sum = weighted_sqr_dists_sum + weights[i] * dists_list[[i]]^2
  }
  
  C = exp(- weighted_sqr_dists_sum / .1)
  
  return(list(C = C, spec_traits = spec_traits))
}

nspec = 10
nres = 5

weights = rep(1 / ndims, ndims)

nsims = 5

ndim_vec = 1:2
richness_mat = matrix(0, nrow = length(ndim_vec), ncol = nsims)

for (j in 1:nsims) {
  for (i in 1:length(ndim_vec)) {
    ndims = ndim_vec[i]
    makeC_list = makeC(nspec, nres, ndims, weights)
    C = makeC_list$C
    C = C
    
    params = list(
      nspec = nspec,
      nres = nres,
      alpha = 0,
      r = 50,
      K = 1,
      m = .2,
      C = C
    )
    
    init_abuns = rep(5, params$nspec)
    init_res = rep(5, params$nres)
    init_state = c(init_abuns, init_res)
    
    sim = ode(y = init_state, times = seq(0, 30000), func = MacArthur, parms = params)
    sim.df = as.data.frame(sim)
    spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
    abuns.df = melt(spec.abuns, id.vars='time')
    eql = tail(sim, 1)[-1]
    eql_abuns = eql[0:nspec]
    num_coexist = length(eql_abuns[eql_abuns > 1])
    
    richness_mat[i,j] = num_coexist
  }
}