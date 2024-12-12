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

find_num_alpha = \(current_state, model, params) {
  nspec = params$nspec
  current_abuns = current_state[0:nspec]
  dt = .01
  dNj = .01
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

nspec = 10
nres = 5

res_trait_1 = seq(0, (nres - 1) / nres, l = nres)
res_trait_2 = seq(0, (nres - 1) / nres, l = nres)

spec_trait_1 = runif(nspec, min = 0, max = 1)
spec_trait_2 = runif(nspec, min = 0, max = 1)
dists_1 = circ_dist(spec_trait_1, res_trait_1)
dists_2 = circ_dist(spec_trait_2, res_trait_2)
w_1 = 1
w_2 = 1 - w_1
C = exp(- ((w_1 * dists_1^2) + (w_2 * dists_2^2)) / .05)

params = list(
  nspec = nspec,
  nres = nres,
  alpha = 0.01,
  r = 100,
  K = 1,
  m = .2,
  C = C
)

init_abuns = rep(5, params$nspec)
init_res = rep(5, params$nres)
init_state = c(init_abuns, init_res)

sim = ode(y = init_state, times = seq(0, 50000), func = MacArthur, parms = params)
sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])

p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')
#print(p)

plot(spec_trait_1, eql_abuns, type = 'h', xlab = 'Species trait 1', ylab = 'Equilibrium Abundance')

num_alpha = find_num_alpha(eql, MacArthur, params)

analytical_alpha = matrix(0, nspec, nspec)
for (i in 1:nspec) {
  for (j in 1:nspec) {
    analytical_alpha[i,j] = - C[i,] %*% C[j,] * params$K / params$r
  }
}

trait_dists_1 = matrix(0, nrow = nspec, ncol = nspec)
trait_dists_2 = matrix(0, nrow = nspec, ncol = nspec)
C_dists = matrix(0, nrow = nspec, ncol = nspec)

for (i in 1:nspec) {
  for (j in 1:nspec) {
    trait_dists_1[i,j] = circ_dist(spec_trait_1[i], spec_trait_1[j])
    trait_dists_2[i,j] = circ_dist(spec_trait_2[i], spec_trait_2[j])
    
    C_dists[i,j] = sqrt(sum((C[i,] - C[j,])^2))
  }
}

trait_dists = sqrt(w_1 * (trait_dists_1)^2 + w_2 * (trait_dists_2)^2)

plot(trait_dists, num_alpha, ylim = c(-.001, 0), xlab = "Species - Species Trait Distances", ylab = "Numerical Alpha", main = "2 Dimensional Trait Space MacArthur Model")

plot(C_dists, num_alpha, ylim = c(-.001, 0), xlab = "C dists", ylab = "Numerical Alpha", main = "1 Dimensional Trait Space MacArthur Model")

C_data = data.frame(C_dists = as.vector(C_dists), num_alpha = as.vector(num_alpha))

plot(eql_abuns, diag(num_alpha), xlim = c(0, 10))

plot(rowSums(C), rowSums(num_alpha))

plot(rowSums(C), eql_abuns)

plot(analytical_alpha, num_alpha, xlab = "Raw Analytical Alpha", ylab = "Raw Numerical Alpha", main = paste0('r = ', as.character(params$r)))
abline(a = 0, b = 1)

plot(analytical_alpha / diag(analytical_alpha), num_alpha / diag(num_alpha), xlab = "Normalized Analytical Alpha", ylab = "Normalized Numerical Alpha")
abline(a = 0, b = 1)