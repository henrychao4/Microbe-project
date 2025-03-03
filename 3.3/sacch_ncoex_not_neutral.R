library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(reshape2)
library(scatterplot3d)
source("KmeansGap.r")

set.seed(3)

create_poly = \(p_length, nm) {
  # vec = sample(nm, p_length, replace = T, prob = c(1/3, 1/3, 1/3))
  vec = sample(nm, p_length, replace = T, prob = c(.325, .33, .335))
  values_to_count = 1:nm
  count_table = table(vec)
  counts = as.vector(count_table[as.character(values_to_count)])
  counts[is.na(counts)] = 0
  return(counts)
}

MacArthur = 
  \(time, state, parms){
    N = state[1:params$nspec]
    R = state[(params$nspec + 1):(params$nspec + params$np)]
    dNdt = with(parms, alpha + N * ((C %*% R) - m))
    dRdt = with(parms, R * (r * (1 - R / K) - t(C) %*% N))
    return(list(c(dNdt, dRdt)))
  }

nspec = 10
nm = 3
np = 5
p_length = 10

nsim = 20
ncoex_vec = rep(0, nsim)

for (k in 1:nsim) {
  polys = matrix(0, nrow = np, ncol = nm)
  
  for (i in 1:np) {
    polys[i,] = create_poly(p_length, nm)
  }
  
  cons_traits_mono = matrix(data = runif(nspec * nm), nrow = nspec, ncol = nm)
  cons_traits_mono = cons_traits_mono
  
  C = matrix(0, nrow = nspec, ncol = np)
  
  for (i in 1:nspec) {
    C[i,] = rowSums(t(t(polys) * cons_traits_mono[i,]))
  }
  
  C = (C * nm) / 50
  
  params = list(
    nspec = nspec,
    np = np,
    alpha = 0,
    r = 500,
    K = 1,
    m = .2,
    C = C
  )
  
  init_abuns = rep(5, params$nspec)
  init_res = rep(5, params$np)
  init_state = c(init_abuns, init_res)
  
  sim = ode(y = init_state, times = seq(0, 15000), func = MacArthur, parms = params)
  sim.df = as.data.frame(sim)
  spec.abuns = sim.df[-((nspec+2):(nspec + np + 2))]
  abuns.df = melt(spec.abuns, id.vars='time')
  eql = tail(sim, 1)[-1]
  eql_abuns = eql[0:nspec]
  num_coexist = length(eql_abuns[eql_abuns > .1])
  
  ncoex_vec[k] = num_coexist
}

hist(ncoex_vec)

# p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')
# 
# 
# 
# ll_trait_data = as.data.frame(cons_traits_mono)
# ll_trait_data$N = round(eql_abuns)
# colnames(ll_trait_data) = c('c1', 'c2', 'c3', 'N')
# print('Abundance of each monosaccharide:')
# print(colSums(polys))
# print(ll_trait_data[order(ll_trait_data$N, decreasing = T),])
# 
# plot(rowSums(cons_traits_mono), eql_abuns, xlab = 'Row Totals for Monosaccharide Preference', ylab = 'Equilibrium Abundance')