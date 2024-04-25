library(ggplot2)
library(tensor)
library(future)
library(parallel)
library(deSolve)
library(tibble)
library(tidyr)
library(dplyr)
library(reshape2)

model = 
  \(t, state, params){
    R = state[1:2]
    N = state[-(1:2)]
    dRdt = with(params, R * (r * (1 - R / K) - t(C) %*% N))
    dNdt = with(params, N * (C %*% R - d))
    return(list(c(dRdt,dNdt)))
  }

nresource = 2
r = c(1, 1)
K = c(1, 1)
res_trait = c(0, 1)
niche_width = c(.1, .1)
d = .2

parms = list(
  r = r,
  res_trait = res_trait,
  niche_width = niche_width
)

inv_growth_rate = 
  \(consumer_trait, res_eql, d, params){
    dists = outer(consumer_trait, params$res_trait, FUN = \(x, y) abs(x - y))
    C = exp(-(dists / params$niche_width) ^ 2)
    dNdt = with(params, (C %*% res_eql - d))
    return(-dNdt)
  }

#solving for an uninvadable community
opt_resident1_trait = optim(.1, inv_growth_rate, params = parms, res_eql = K, d = d, method = 'Brent', lower = 0, upper = 1)$par
opt_resident2_trait = 1 - opt_resident1_trait

dists = outer(c(opt_resident1_trait, opt_resident2_trait), res_trait, FUN = \(x, y) abs(x - y))
res_C = exp(-(dists / niche_width) ^ 2)
res_d = c(.2, .2)

init_state = c(rep(1, nresource), rep(1, nrow(res_C)))

sim = ode(y = init_state, times = seq(0, 1000, by = .1), func = model, parms = list(r = r, K = K, C = res_C, d = res_d))
resident_eql = tail(sim, 1)[-1]
R_star = resident_eql[1:2]

d_3 = 
  \(inv_trait, res_C, res_d, params){
    dists = outer(inv_trait, params$res_trait, FUN = \(x, y) abs(x - y))
    inv_C = exp(-(dists / params$niche_width) ^ 2)
    d = inv_C[1] * ((res_C[2,2] * res_d[1] - res_C[1,2] * res_d[2]) / (res_C[1,1] * res_C[2,2] - res_C[1,2] * res_C[2,1])) + inv_C[2] * ((res_C[1,1] * res_d[2] - res_C[2,1] * res_d[1]) / (res_C[1,1] * res_C[2,2] - res_C[1,2] * res_C[2,1]))
    return(d)
  }

test_traits = seq(0, 1, by = .01)
d = rep(0, length(test_traits))
for (i in 1:length(test_traits)) {
  d[i] = d_3(test_traits[i], res_C, res_d, params = parms)
}
# y = -inv_growth_rate(test_traits, R_star, params = parms, d = d)
plot(test_traits, d)