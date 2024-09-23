library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
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

nspec = 20
nres = 3

res_trait = seq(.2, .8, l = nres)
spec_trait = seq(0, .9, l = nspec) + rnorm(nspec, mean = 0, sd = .005)
dists = circ_dist(spec_trait, res_trait)
C = exp(-dists^2 / .1)

params = list(
  nspec = nspec,
  nres = nres,
  alpha = .005,
  r = 50,
  K = 1,
  m = .2,
  C = C
)

init_state = c(rep(1, params$nspec), rep(1, params$nres))

sim = ode(y = init_state, times = seq(0, 15000), func = MacArthur, parms = params)
sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])

p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')
#print(p)

plot(spec_trait, eql_abuns, type = 'h')

init_data = as.data.frame(C)
init_data$N = rep(1, nspec)
#init_data$N = sample(c(1,2,3), nspec, replace = T)

eql_data = as.data.frame(C)
eql_data$N = round(eql_abuns)

#because the method is permutation shuffling, gap curve is flat
init_gap = KmeansGap(dat = init_data, multiD = T, mink = 1, maxk = 9)

eql_gap = KmeansGap(dat = eql_data, multiD = T, mink = 1, maxk = 9)

plot(init_gap$data$k, init_gap$data$gap, type = 'b')

plot(eql_gap$data$k, eql_gap$data$gap, type = 'b')