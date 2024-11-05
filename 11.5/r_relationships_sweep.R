library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(reshape2)
library(lsa)
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

nspec = 30
nres = 30

ntraits_vec = 1:5
nreps = 3
p_val = matrix(0, nrow = length(ntraits_vec), ncol = nreps)

for (i in 1:length(ntraits_vec)) {
  ntraits = ntraits_vec[i]
  
  for (j in 1:nreps) {
    res_traits = matrix(runif(nres * ntraits, min = 0, max = 1), nrow = nres, ncol = ntraits)
    res_dists = as.matrix(dist(res_traits, upper = T, diag = T))
    res_sim = exp(- res_dists^2 / .2)
    C = res_sim
    
    params = list(
      nspec = nspec,
      nres = nres,
      alpha = 0.05,
      r = 50,
      K = 10,
      m = .2,
      C = C
    )
    
    init_abuns = rep(5, params$nspec)+ rnorm(params$nspec, mean = 0, sd = 0.01)
    init_res = rep(5, params$nres)+ rnorm(params$nspec, mean = 0, sd = 0.01)
    init_state = c(init_abuns, init_res)
    
    sim = ode(y = init_state, times = seq(0, 1000, by = .01), func = MacArthur, parms = params, method = "rk4")
    sim.df = as.data.frame(sim)
    spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
    abuns.df = melt(spec.abuns, id.vars='time')
    eql = tail(sim, 1)[-1]
    eql_abuns = eql[0:nspec]
    
    p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')
    
    kmg_data = as.data.frame(C)
    kmg_data$N = round(eql_abuns)
    kmg_gap = KmeansGap(dat = kmg_data, multiD = T, mink = 1, maxk = 10)
    
    p_val[i,j] = kmg_gap$p.value
    
  }
}

p_vals_vec = apply(p_val, 1, median)

plot(ntraits_vec, p_vals_vec, ylim = c(0,1), type = 'b', xlab = 'Number of resource trait dimensions', ylab = 'p-value for clustering')
