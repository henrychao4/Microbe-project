library(tidyverse)
library(deSolve)
library(reshape2)

hoi_effects = \(beta, y) {
  ans = rep(0, length(y))
  for (i in 1:length(y)) {
    for (j in 1:length(y)) {
      for (k in 1:length(y)) {
        ans[i] = ans[i] + beta[i,j,k] * y[j] * y[k]
      }
    }
  }
  return(ans)
}

hoi_model = function(t, y, params) {
  dN = y * (params$r + params$alpha %*% y + hoi_effects(beta, y))
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

beta = array(rnorm(n^3, 0, 0), dim = c(n, n, n))


init_abuns = rep(1, n)

params = list(alpha = alpha, r = r, tmax = tmax, step = step)

sim = ode(y = init_abuns, times = seq(0,params$tmax, by = params$step), func = hoi_model, parms = params)

sim.df = as.data.frame(sim)
abuns.df = melt(sim.df, id.vars='time')
p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
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
  fwd_pert_state = ode(y = pert_abuns, times = seq(0,dt, by = dt/10), func = hoi_model, parms = params)
  next_pert_abuns = tail(fwd_pert_state, n = 1)
  next_pert_abuns = as.numeric(next_pert_abuns[1:n+1])
  num_alpha[,j] = (next_pert_abuns - pert_abuns) / (dt * dNj * pert_abuns) - (next_abuns - current_abuns) / (dt * dNj * current_abuns)
}

plot(alpha, num_alpha, xlab = expression(paste("True ", alpha, "s from LV Model")), ylab = expression(paste("Numerical ", alpha, "s")))
abline(a = 0, b = 1)