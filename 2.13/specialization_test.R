library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(reshape2)
library(scatterplot3d)
source("KmeansGap.r")

set.seed(1)

MacArthur = 
  \(time, state, parms){
    N = state[1:params$nspec]
    R = state[(params$nspec + 1):(params$nspec + params$nres)]
    dNdt = with(parms, alpha + N * ((C %*% R) - m))
    dRdt = with(parms, R * (r * (1 - R / K) - t(C) %*% N))
    return(list(c(dNdt, dRdt)))
  }

nspec = 30
nres = 3

specialization = .6

C = matrix(data = runif(nspec * nres), nrow = nspec, ncol = nres)

for (i in 1:nspec) {
  idx = floor(nres * ((i - 1)/nspec)) + 1
  C[i, idx] = 0
}
C = (C / rowSums(C)) * (1 - specialization)
for (i in 1:nspec) {
  idx = floor(nres * ((i - 1)/nspec)) + 1
  C[i, idx] = specialization
}

params = list(
  nspec = nspec,
  nres = nres,
  alpha = 0,
  r = c(40, 50, 60),
  K = 10,
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

p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle(paste0('Specializaion = ', as.character(specialization)))

print(p)