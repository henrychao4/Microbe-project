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
res_c = matrix(c(.2, .1, .1, .6), nrow = nresident, ncol = nresource)
inv_c = c(.4, .3)

res_d = c(.2, .2)
inv_d = inv_c[1] * (res_c[2,2] * res_d[1] - res_c[1,2] * res_d[2]) / (res_c[1,1] * res_c[2,2] - res_c[1,2] * res_c[2,1]) + inv_c[2] * (res_c[1,1] * res_d[2] - res_c[2,1] * res_d[1]) / (res_c[1,1] * res_c[2,2] - res_c[1,2] * res_c[2,1])
inv_d = inv_d + .001
C = rbind(res_c, inv_c)
d = c(res_d, inv_d)
r = c(1, 1)
K = c(10, 10)

params = list(
  r = r,
  K = K,
  C = C,
  d = d
)

init_state = c(rep(10, length(r)), rep(10, length(res_d)))
initialize = ode(y = init_state, times = seq(0, 10000, by = .1), func = model, parms = list(r = r, K = K, C = res_c, d = res_d))
resident_eql = tail(initialize, 1)[-1]


invader_init = 1
new_init_state = c(resident_eql, invader_init)

sim = ode(y = new_init_state, times = seq(0, 1000, by = .1), func = model, parms = params) |>
  as.data.frame()

colnames(sim) = c('time', 'R1', 'R2', 'Nres1', 'Nres2', paste0('N', 1:ninvader))

#t_extinct = min(filter(sim, N1 < .01)$time)

sim.df = as.data.frame(sim)
spec.abuns = sim.df[-(2:3)]
abuns.df = melt(spec.abuns, id.vars='time')
p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
print(p)