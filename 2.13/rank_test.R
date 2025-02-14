library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(reshape2)
library(scatterplot3d)
library(Matrix)
source("KmeansGap.r")

set.seed(3)

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

nspec = 3
nm = 2
np = 3
p_length = 10

polys = matrix(0, nrow = np, ncol = nm)

for (i in 1:np) {
  polys[i,] = create_poly(p_length, nm)
}

cons_traits_mono = matrix(data = runif(nspec * nm), nrow = nspec, ncol = nm)
cons_traits_mono = cons_traits_mono / rowSums(cons_traits_mono)

scatterplot3d(cons_traits_mono, pch = 16, color = "blue", main = "3D Scatter Plot",
              xlab = "Affinity for Mono 1", ylab = "Affinity for Mono 2", zlab = "Affinity for Mono 3")

C = matrix(0, nrow = nspec, ncol = np)

for (i in 1:nspec) {
  C[i,] = rowSums(t(t(polys) * cons_traits_mono[i,]))
  #C[i,] = rowSums(t(t(polys) * cons_traits_mono[i,]^1.5))
}

C = C / (p_length)

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

plot(cons_traits_mono[,1], eql_abuns, type = 'h')

p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')

print(rankMatrix(C))