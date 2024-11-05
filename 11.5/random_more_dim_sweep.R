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

sim = \(params) {
  init_abuns = rep(5, params$nspec)
  init_res = rep(5, params$nres)
  init_state = c(init_abuns, init_res)
  
  sim = ode(y = init_state, times = seq(0, 5000), func = MacArthur, parms = params)
  sim.df = as.data.frame(sim)
  spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
  abuns.df = melt(spec.abuns, id.vars='time')
  eql = tail(sim, 1)[-1]
  eql_abuns = eql[0:nspec]
  
  
  return(eql_abuns)
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

ndim_vec = 1:5
nreps = 3
p_vals_more_dims = matrix(0, nrow = length(ndim_vec), ncol = nreps)
p_vals_2dims = matrix(0, nrow = length(ndim_vec), ncol = nreps)

for (i in 1:length(ndim_vec)) {
  nspec = 300
  nres = 5
  ndims = ndim_vec[i]
  weights = rep(1 / ndims, ndims)
  
  weights_2 = c(weights[i], 1 - weights[i])
  
  for (j in 1:nreps) {
    makeC_list = makeC(nspec, nres, ndims, weights)
    C = makeC_list$C
    makeC_list_2 = makeC(nspec, nres, 2, weights_2)
    C_2 = makeC_list_2$C
    
    params = list(
      nspec = nspec,
      nres = nres,
      alpha = 0.01,
      r = 1000,
      K = 1,
      m = .2,
      C = C
    )
    
    params_2 = list(
      nspec = nspec,
      nres = nres,
      alpha = 0.01,
      r = 1000,
      K = 1,
      m = .2,
      C = C_2
    )
    
    eql_abuns = sim(params)
    eql_abuns_2 = sim(params_2)
    
    plot(makeC_list_2$spec_traits[1,], eql_abuns, type = 'h', xlab = 'Species trait 1', ylab = 'Equilibrium Abundance',
         main = paste0('Number of dimensions = ', as.character(ndims)))
    
    kmg_data = as.data.frame(C)
    kmg_data$N = round(eql_abuns)
    kmg_gap = KmeansGap(dat = kmg_data, multiD = T, mink = 1, maxk = 10)
    
    kmg_data_2 = as.data.frame(C_2)
    kmg_data_2$N = round(eql_abuns_2)
    kmg_gap_2 = KmeansGap(dat = kmg_data_2, multiD = T, mink = 1, maxk = 10)   
    
    p_vals_more_dims[i,j] = kmg_gap$p.value
    p_vals_2dims[i,j] = kmg_gap_2$p.value
  }

}

p_vals_more_dims_vec = apply(p_vals_more_dims, 1, median)
p_vals_2dims_vec = apply(p_vals_2dims, 1, median)

plot(ndim_vec, p_vals_more_dims_vec, ylim = c(0,1), type = 'b', xlab = 'Number of other trait dimensions', ylab = 'p-value for clustering', col = 'red')
lines(ndim_vec, p_vals_2dims_vec, type = 'b', col = 'blue')