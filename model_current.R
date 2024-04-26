library(ggplot2)
library(tensor)
library(future)
library(parallel)
library(deSolve)
library(tibble)
library(tidyr)
library(dplyr)
library(reshape2)

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

set.seed(1)

nspec = 10
nres = 10

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(0, nspec),
  mu = rep(.05, nspec),
  rho = rep(.3, nres),
  delta = rep(.05, nres),
  epsilon = .9,
  beta = 1,
  h = 1,
  #I = diag(nspec),
  I = matrix(rbinom(nspec * nres, 1, 0.5), nrow = nspec, ncol = nres),
  p = array(rep(0, nspec * nres^2), dim=c(nspec, nres, nres))
)

init_state = c(rep(1, params$nspec), rep(1, params$nres), 0)

U = 
  \(params, R){
    ans = matrix(0, nrow = params$nspec, ncol = params$nres)
    for (j in 1:params$nres) {
      for (i in 1:params$nspec) {
        ans[i,j] = params$beta * R[j] / (params$h + R[j]) * params$I[i,j]
      }
    }
    return(ans)
  }

omega = 
  \(params, R){
    ans = matrix(0, nrow = params$nspec, ncol = params$nres)
    for (i in 1:params$nspec) {
      ans[i,] = U(params, R)[i,] %*% params$p[i,,]
    }
    return(ans)
  }

gamma = 
  \(params, R){
    ans = rep(0, params$nspec)
    for (i in 1:params$nspec) {
      ans[i] = params$epsilon * sum(U(params, R)[i,])
    }
    ans2 = params$epsilon * rowSums(U(params, R))
    return(c(ans, ans2))
  }

find_num_alpha =
  \(eql, eql_abuns, model, params){
    dt = .001
    dNj = .001
    fwd_result = ode(y = eql, times = seq(0, dt, by = dt/100), func = model, parms = params)
    fwd_state = as.numeric(tail(fwd_result, n = 1)[-1])
    fwd_abuns = fwd_state[0:nspec]
    num_alpha = matrix(0, params$nspec, params$nspec)
    for (j in 1:nspec) {
      pert_state = eql
      pert_state[j] = pert_state[j] + dNj
      pert_abuns = pert_state[0:nspec]
      fwd_pert_result = ode(y = pert_state, times = seq(0,dt, by = dt/100), func = model, parms = params)
      fwd_pert_state = as.numeric(tail(fwd_pert_result, n = 1)[-1])
      fwd_pert_abuns = fwd_pert_state[0:nspec]
      num_alpha[,j] = (fwd_pert_abuns - pert_abuns) / (dt * dNj * pert_abuns) - (fwd_abuns - eql_abuns) / (dt * dNj * eql_abuns)
    }
    return(num_alpha)
  }

model = 
  \(t, state, params){
    N = state[1:params$nspec]
    R = state[(params$nspec + 1):(params$nspec + params$nres)]
    W = tail(state, 1)
    dNdt = params$alpha + gamma(params, R) * N - params$mu * N
    dRdt = params$rho - params$delta * R + colSums(omega(params, R)) - colSums(U(params, R) * N)
    dWdt = (1 - params$epsilon) * sum(U(params, R) * N)
    
    return(list(c(dNdt, dRdt, dWdt)))
  }

sim = ode(y = init_state, times = seq(0, 5, by = 1), func = model, parms = params)

sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
print(p)

eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])
df = c(eql_abuns, num_coexist) |> t() |> as_tibble()
colnames(df) = c(paste0('N', seq(nspec)), 'num_coexist')

hamming_dist = \(x, y) {
  ans = sum(x != y)
  return(ans)
}
