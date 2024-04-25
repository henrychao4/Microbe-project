library(ggplot2)
library(tensor)
library(future)
library(parallel)
library(deSolve)
library(tibble)
library(tidyr)
library(dplyr)
library(reshape2)

set.seed(1)

model = 
  \(t, state, params){
    R = state[1:2]
    N = state[-(1:2)]
    dRdt = with(params, R * (r * (1 - R / K) - t(C) %*% N))
    dNdt = with(params, N * (C %*% R - d))
    return(list(c(dRdt,dNdt)))
  }

nresident = 2
ninvader = 1
nresource = 2
res_c = matrix(c(.2, .1, .1, .2), nrow = nresident, ncol = nresource)
inv_c = .5 * matrix(runif(ninvader * nresource), nrow = ninvader, ncol = nresource)

res_d = c(.2, .2)
pert = .1 * runif(ninvader)
inv_d = inv_c[,1] * (res_c[2,2] * res_d[1] - res_c[1,2] * res_d[2]) / (res_c[1,1] * res_c[2,2] - res_c[1,2] * res_c[2,1]) + inv_c[,2] * (res_c[1,1] * res_d[2] - res_c[2,1] * res_d[1]) / (res_c[1,1] * res_c[2,2] - res_c[1,2] * res_c[2,1]) + pert
C = rbind(res_c, inv_c)
d = c(res_d, inv_d)
r = c(1,1)
K = c(10,10)

params = list(
  r = r,
  K = K,
  C = C,
  d = d
)

init_state = c(rep(10, nresource), rep(10, nresident), rep(10, ninvader))

sim = ode(y = init_state, times = seq(0, 1000, by = .001), func = model, parms = params) |>
  as.data.frame()

colnames(sim) = c('time', 'R1', 'R2', 'Nres1', 'Nres2', paste0('N', 1:ninvader))

t_extinct = rep(0, ninvader)

for (i in 1:length(t_extinct)) {
  spec = as.name(paste0('N', toString(i)))
  t_extinct[i] = min(filter(sim, !!spec  < .01)$time)
}

plot(pert, t_extinct)