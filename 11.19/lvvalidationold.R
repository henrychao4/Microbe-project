library(tidyverse)
library(deSolve)
library(reshape2)

lv = function(t, y, params) {
  dN = y * (params$r + params$alpha %*% y)
  return(list(c(dN)))
}

set.seed(1)
tmax = 3000
step = .1
n = 5
r = rep(.1, n)
alpha = matrix(-.25, n, n)
diag(alpha) = diag(alpha) - .1
alpha = alpha + matrix(rnorm(n^2,0,.1),nrow = n)
init_abuns = rep(1, n)

params = list(alpha = alpha, r = r, tmax = tmax, step = step)

sim = ode(y = init_abuns, times = seq(0,params$tmax, by = params$step), func = lv, parms = params)

sim.df = as.data.frame(sim)
abuns.df = melt(sim.df, id.vars='time')
p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
print(p)

current_state = tail(sim, n = 1)
current_abuns = as.numeric(current_state[1:n+1])
dt = .001
dNj = .001
fwd_state = ode(y = current_abuns, times = seq(0,dt, by = dt/10), func = lv, parms = params)
next_abuns = tail(fwd_state, n = 1)
next_abuns = as.numeric(next_abuns[1:n+1])

num_alpha = matrix(0, n, n)
for (j in 1:n) {
  pert_abuns = current_abuns
  pert_abuns[j] = pert_abuns[j] + dNj
  fwd_pert_state = ode(y = pert_abuns, times = seq(0,dt, by = dt/10), func = lv, parms = params)
  next_pert_abuns = tail(fwd_pert_state, n = 1)
  next_pert_abuns = as.numeric(next_pert_abuns[1:n+1])
  num_alpha[,j] = (next_pert_abuns - pert_abuns) / (dt * dNj * pert_abuns) - (next_abuns - current_abuns) / (dt * dNj * current_abuns)
}

plot(alpha, num_alpha)
abline(a = 0, b = 1)