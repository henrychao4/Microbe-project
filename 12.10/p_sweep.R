library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(reshape2)
library(randomForest)
source("KmeansGap.r")

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

set.seed(1)

cross_feed = \(p, uptake) {
  tot = matrix(0, nrow = nrow(uptake), ncol = ncol(uptake))
  for (i in 1:nrow(uptake)) {
    tot = tot + t(p[i,,]) %*% uptake
  }
  return(rowSums(tot))
}

MacArthur =
  \(time, state, params) {
    N = state[1:params$nspec]
    R = state[(params$nspec + 1):(params$nspec + params$nres)]
    
    uptake = t(t(C) * R) * N
    
    dNdt = params$alpha + ((1 - rowSums(rowSums(params$p, dims = 2))) * params$C %*% R - params$m) * N
    dRdt = R * (params$r * (1 - R / params$K) - t(params$C) %*% N) + cross_feed(params$p, uptake)
    
    return(list(c(dNdt, dRdt)))
  }

circ_dist = \(vec1, vec2) {
  y = abs(outer(vec1, vec2, '-'))
  idx = which(y > 0.5)
  y[idx] = 1 - y[idx]
  return(y)
}

find_num_alpha = \(current_state, model, params) {
  nspec = params$nspec
  current_abuns = current_state[0:nspec]
  dt = .1
  dNj = .1
  fwd_result = ode(y = current_state, times = seq(0,dt, by = dt/100), parms = params, func = MacArthur)
  fwd_state = as.numeric(tail(fwd_result, n = 1)[-1])
  fwd_abuns = current_state[0:nspec]
  num_alpha = matrix(0, nspec, nspec)
  for (j in 1:nspec) {
    pert_state = current_state
    pert_state[j] = pert_state[j] + dNj
    pert_abuns = pert_state[0:nspec]
    fwd_pert_result = ode(y = pert_state, times = seq(0, dt, by = dt/100), parms = params, func = MacArthur)
    fwd_pert_state = as.numeric(tail(fwd_pert_result, n = 1)[-1])
    fwd_pert_abuns = fwd_pert_state[0:nspec]
    num_alpha[,j] = (fwd_pert_abuns - pert_abuns) / (dt * dNj * pert_abuns) - (fwd_abuns - current_abuns) / (dt * dNj * current_abuns)
  }
  return(num_alpha)
}

nspec = 2
nres = 2

res_trait = seq(0, (nres - 1) / nres, l = nres)

spec_trait = seq(0, (nspec - 1) / nspec, l = nspec)
dists = circ_dist(spec_trait, res_trait)

C = matrix(c(1, 0, .9, 1), nrow = 2, ncol = 2)

p = array(rep(0, nspec * nres^2), dim = c(nspec, nres, nres))

p_vec = seq(.51, .52, by = .001)
alpha_vec = rep(0, length(p_vec))

for (i in 1:length(p_vec)) {
  
  p[1,1,2] = p_vec[i]
  
  r = c(10, 10)
  
  params = list(
    nspec = nspec,
    nres = nres,
    alpha = 0,
    r = r,
    K = 1,
    m = .2,
    C = C,
    p = p
  )
  
  init_abuns = rep(5, params$nspec)
  init_res = rep(5, params$nres)
  init_state = c(init_abuns, init_res)
  
  sim = ode(y = init_state, times = seq(0, 5000), func = MacArthur, parms = params)
  sim.df = as.data.frame(sim)
  spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
  abuns.df = melt(spec.abuns, id.vars='time')
  eql = tail(sim, 1)[-1]
  eql_abuns = eql[0:nspec]
  num_coexist = length(eql_abuns[eql_abuns > .1])
  
  num_alpha = find_num_alpha(eql, MacArthur, params)
  
  alpha_vec[i] = num_alpha[2,1]
}

plot(p_vec, alpha_vec, xlab = 'p_{112}', ylab = 'a_{12}', main = 'k = .4')
abline(a = 0, b = 0)

x_vec = seq(.1, .9, by = .1)
y_vec = c(.116, .205, .274, .31, .378, .419, .454, .485, .513)
plot(x_vec, y_vec, type = 'b', xlab = 'c_{12}', ylab = 'p_{112}', main = 'Parameter Combinations for alpha_{12} = 0')