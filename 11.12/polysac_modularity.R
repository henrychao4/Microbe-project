library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(reshape2)
library(vegan)
library(wsyn)
source("KmeansGap.r")

set.seed(1)

create_poly = \(p_length, nm) {
  vec = sample(nm, p_length, replace = T)
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

nspec = 100
nm = 3
np = 10
p_length = 10

polys = matrix(0, nrow = np, ncol = nm)

for (i in 1:np) {
  polys[i,] = create_poly(p_length, nm)
}

adj = vegdist(polys, method = "bray", diag = T, upper = T)
#modularity(adj)


cons_traits_mono = matrix(data = runif(nspec * nm), nrow = nspec, ncol = nm)
#cons_traits_mono = cons_traits_mono / rowSums(cons_traits_mono)

C = matrix(0, nrow = nspec, ncol = np)

for (i in 1:nspec) {
  C[i,] = rowSums(t(t(polys) * cons_traits_mono[i,]))
}

C = C / p_length

params = list(
  nspec = nspec,
  np = np,
  alpha = 0.05,
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

p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')

kmg_data = as.data.frame(C)
kmg_data$N = round(eql_abuns)
kmg_gap = KmeansGap(dat = kmg_data, multiD = T, mink = 1, maxk = 5)

plot(kmg_gap$data$k, kmg_gap$data$gap, type = 'b')

print(kmg_gap)