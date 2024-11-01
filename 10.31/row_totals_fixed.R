library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(reshape2)
library(MEDseq)
source("KmeansGap.r")

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

MacArthur = \(time, state, parms){
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

set.seed(1)

nspec = 1000
nres = 3

nreps = 1
p_vals = rep(0, nreps)
opt_k = rep(0, nreps)

for (i in 1:nreps) {
  
  # res_trait_1 = seq(0, (nres - 1) / nres, l = nres)
  # spec_trait_1 = seq(0, (nspec - 1) / nspec, l = nspec)
  # dists_1 = circ_dist(spec_trait_1, res_trait_1)
  # C = exp(-dists_1^2 / .1)
  
  C = matrix(runif(nspec * nres, min = 0, max = 1), nrow = nspec, ncol = nres)
  
  row_totals = rowSums(C)
  min_row_tot = min(row_totals)
  
  for (j in 1:nrow(C)) {
    C[j,] = C[j,] * (min_row_tot / sum(C[j,]))
  }
  
  C = C / rowSums(C)
  
  params = list(
    nspec = nspec,
    nres = nres,
    alpha = 0,
    r = 1000,
    K = 1,
    m = .2,
    C = C
  )
  
  init_abuns = rep(5, params$nspec)
  init_res = rep(5, params$nres)
  init_state = c(init_abuns, init_res)
  
  sim = ode(y = init_state, times = seq(0, 15000), func = MacArthur, parms = params)
  sim.df = as.data.frame(sim)
  spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
  abuns.df = melt(spec.abuns, id.vars='time')
  eql = tail(sim, 1)[-1]
  eql_abuns = eql[0:nspec]
  num_coexist = length(eql_abuns[eql_abuns > .1])
  
  p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')
  #print(p)
  
  plot(eql_abuns, type = 'h')
  
  k_max = 10
  
  kmg_data = as.data.frame(C)
  kmg_data$N = round(eql_abuns)
  kmg_gap = KmeansGap(dat = kmg_data, multiD = T, mink = 1, maxk = k_max, numnulls = 100)
  
  p_vals[i] = kmg_gap$p.value
  opt_k[i] = kmg_gap$khat
}

print(p_vals)
print(opt_k)