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

nspec = 100
nres = 5

res_traits = seq(0, (nres - 1) / nres, l = nres)
spec_traits = seq(0, (nspec - 1) / nspec, l = nspec)
#spec_traits = runif(nspec, min = 0, max = 1)
dists = circ_dist(spec_traits, res_traits)
C = exp(-dists^2 / .05)

params = list(
  nspec = nspec,
  nres = nres,
  alpha = 0.01,
  r = 3 * nspec,
  K = 1,
  m = .2,
  C = C
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

p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')
#print(p)

plot(spec_traits, eql_abuns, type = 'h', xlab = 'Species trait 1', ylab = 'Equilibrium Abundance')

kmg_data = as.data.frame(C)
kmg_data$N = round(eql_abuns)
kmg_gap = KmeansGap(dat = kmg_data, multiD = T, mink = 1, maxk = 10)

plot(kmg_gap$data$k, kmg_gap$data$gap, type = 'b', xlab = 'k', ylab = 'Gap')

print(kmg_gap)