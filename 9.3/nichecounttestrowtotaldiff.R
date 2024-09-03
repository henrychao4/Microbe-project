library(ggplot2)
library(tensor)
library(future)
library(parallel)
library(deSolve)
library(tibble)
library(tidyr)
library(dplyr)
library(reshape2)
library(purrr)
library(rootSolve)

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

set.seed(1)

nspec = 20
nres = 20
nruns = 10

makeI = \(nspec, nres, nconsume) {
  I = matrix(0, nrow = nspec, ncol = nres)
  cons_vec = c(rep(1, nconsume), rep(0, (nres - nconsume)))
  for (i in 1:nspec) {
    I[i,] = sample(cons_vec, nres, replace = T)
  }
  return(I)
}

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(0, nspec),
  mu = rep(.05, nspec),
  rho = rep(.3, nres),
  delta = rep(.05, nres),
  # mu = rep(0, nspec),
  # rho = rep(0, nres),
  # delta = rep(0, nres),
  epsilon = 1,
  beta = 1,
  h = 1,
  I = I,
  p = array(rep(0, nspec * nres^2), dim=c(nspec, nres, nres))
)

U = 
  \(params, R){
    ans = t(params$beta * R / (params$h + R) * t(params$I))
    return(ans)
  }

omega = 
  \(params, R){
    ans = rowSums(params$p, dims = 2) * U(params, R)
    return(ans)
  }

gamma = 
  \(params, R){
    ans = params$epsilon * rowSums((1 - rowSums(params$p, dims = 2)) * U(params, R))
    return(ans)
  }

model = 
  \(t, state, params){
    N = state[1:params$nspec]
    R = state[(params$nspec + 1):(params$nspec + params$nres)]
    W = tail(state, 1)
    dNdt = params$alpha + gamma(params, R) * N - params$mu * N
    dRdt = params$rho - params$delta * R + colSums(omega(params, R) * N) - colSums(U(params, R) * N)
    dWdt = sum((1 - params$epsilon) * (1 - rowSums(params$p, dims = 2)) * U(params, R) * N)
    
    return(list(c(dNdt, dRdt, dWdt)))
  }

nconsume_vec = seq(1, nres)
num_coexist_vec = rep(0, length(nconsume_vec))
lb = rep(0, length(nconsume_vec))
ub = rep(0, length(nconsume_vec))
init_state = c(rep(1, params$nspec), rep(1, params$nres), 0)

for (i in 1:length(nconsume_vec)) {
  nc = rep(0, nruns)
  for (j in 1:nruns) {
    nconsume = nconsume_vec[i]
    params$I = makeI(nspec, nres, nconsume)
    eql = runsteady(y = init_state, time = c(0, 500000), func = model, parms = params, stol = 1e-8)
    eql_abuns = eql$y[0:nspec]
    nc[j] = length(eql_abuns[eql_abuns > .1])
  }
  num_coexist_vec[i] = mean(nc)
  lb[i] = quantile(nc, probs = .025)
  ub[i] = quantile(nc, probs = .975)
}

df = data.frame(nconsume_vec = nconsume_vec, num_coexist_vec = num_coexist_vec, lb = lb, ub = ub)

p = ggplot(df, aes(x = nconsume_vec, y = num_coexist_vec)) + geom_line() + geom_point() +
  geom_ribbon(aes(ymin = lb, ymax = ub), fill = 'blue', color = 'black', linetype = "dotted", alpha = .1)
print(p)