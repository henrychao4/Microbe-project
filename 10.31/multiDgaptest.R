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

circ_dist = function(vec1, vec2) {
  y = abs(outer(vec1, vec2, '-'))
  idx = which(y > 0.5)
  y[idx] = 1 - y[idx]
  return(y)
}

nspec = 81
nres = 3

res_trait_1 = seq(0, (nres - 1) / nres, l = nres)
res_trait_2 = seq(0, (nres - 1) / nres, l = nres)

spec_traits = seq(0, (sqrt(nspec) - 1) / sqrt(nspec), l = sqrt(nspec))
spec_trait_1 = rep(spec_traits, each = sqrt(nspec))
spec_trait_2 = rep(spec_traits, sqrt(nspec))

#spec_trait_2 = sample(spec_trait_2, size = length(spec_trait_2), replace = F)
dists_1 = circ_dist(spec_trait_1, res_trait_1)
dists_2 = circ_dist(spec_trait_2, res_trait_2)
w_1 = .5
w_2 = 1 - w_1
C = exp(- ((w_1 * dists_1^2) + (w_2 * dists_2^2)) / .1)

params = list(
  nspec = nspec,
  nres = nres,
  alpha = .05,
  r = 500,
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
# plot(spec_trait_1, eql_abuns, type = 'h')
# plot(spec_trait_2, eql_abuns, type = 'h')

kmg_data = as.data.frame(C)
kmg_data$N = round(eql_abuns)
kmg_gap = KmeansGap(dat = kmg_data, multiD = T, mink = 1, maxk = 10)

plot(kmg_gap$data$k, kmg_gap$data$gap, type = 'b', main = "wx = .5, wy = .5")

print(kmg_gap)