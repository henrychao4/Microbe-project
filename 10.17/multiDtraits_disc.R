library(ggplot2)
library(tensor)
library(future)
library(parallel)
library(deSolve)
library(tibble)
library(tidyr)
library(dplyr)
library(reshape2)
library(cluster)
library(MEDseq)
library(purrr)
library(klaR)

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

set.seed(3)

nspec = 50
nres = 30

res_trait_1 = runif(nres, min = 0, max = (nres - 1) / nres)
res_trait_2 = runif(nres, min = 0, max = (nres - 1) / nres)
spec_trait_1 = seq(0, 1, l = nspec) + rnorm(nspec, mean = 0, sd = .005)
spec_trait_1 = sample(spec_trait_1, size = length(spec_trait_1), replace = F)
spec_trait_2 = seq(0, 1, l = nspec) + rnorm(nspec, mean = 0, sd = .005)
spec_trait_2 = sample(spec_trait_2, size = length(spec_trait_2), replace = F)
dists_1 = outer(spec_trait_1, res_trait_1, FUN = \(x, y) abs(x - y))
dists_2 = outer(spec_trait_2, res_trait_2, FUN = \(x, y) abs(x - y))
sigma_1 = .5
sigma_2 = .5
C = exp(-((dists_1 / sigma_1)^2 + (dists_2 / sigma_2)^2)^2)

I = (C > .3) * 1

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(0.05, nspec),
  mu = rep(.05, nspec),
  rho = rnorm(nres, mean = .6, sd = 0) + c(rep(0,9), .4, rep(0,9), .4, rep(0,9), .4), #3 good res
  #rho = rnorm(nres, mean = .6, sd = 0) + c(rep(0,4), .4, rep(0,4), .4, rep(0,4), .4, rep(0,4), .4, rep(0,4), .4, rep(0,4), .4), #6 good res
  delta = rep(.05, nres),
  epsilon = .9,
  beta = 1,
  h = 1,
  I = I,
  p = array(rep(0, nspec * nres^2), dim = c(nspec, nres, nres))
)

init_state = c(rep(1, params$nspec), rep(1, params$nres), 0)

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

best_wKModes = \(data, modes, weights, nruns) {
  df = list()
  tot.withindiff_vec = rep(0, nruns)
  for (i in 1:nruns) {
    wk_modes = wKModes(data, modes = modes, weights = weights)
    df[[i]] = wk_modes
    tot.withindiff_vec[i] = wk_modes$tot.withindiff
  }
  idx = which.min(tot.withindiff_vec)
  return(df[[idx]])
}

sim = ode(y = init_state, times = seq(0, 30000, by = 1), func = model, parms = params)

sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
#print(p)

eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])
df = c(eql_abuns, num_coexist) |> t() |> as_tibble()
colnames(df) = c(paste0('N', seq(nspec)), 'num_coexist')

plot(spec_trait_1, eql_abuns, type = 'h')
plot(spec_trait_2, eql_abuns, type = 'h')

k_max = 7
true_errs = rep(0, k_max)
null_errs = rep(0, k_max)
nboot = 30

for (k in 1:k_max) {
  kmode = best_wKModes(I, modes = k, weights = eql_abuns, nruns = 10)
  true_errs[k] = sum(kmode$withindiff)
}

replicate_err = matrix(0, nrow = nboot, ncol = k_max)
for (i in 1:nboot) {
  null_abuns = sample(eql_abuns, size = length(eql_abuns), replace = F)
  for (k in 1:k_max) {
    null_kmode = best_wKModes(I, modes = k, weights = null_abuns, nruns = 10)
    replicate_err[i,k] = sum(null_kmode$withindiff)
  }
  print(i)
}

null_errs = apply(replicate_err, 2, mean)

null_gaps = matrix(0, nrow = nboot, ncol = k_max)
for (k in 1:k_max) {
  null_gaps[,k] = null_errs[k] - replicate_err[,k]
}

null_max_gaps = apply(null_gaps, 1, max)
null_max_gap_q95 = quantile(null_max_gaps, .95)

gap = null_errs - true_errs

plot(1:k_max, gap, type = 'b', xlab = 'k', ylab = 'Gap')

print(paste0('True max gap for equilibrium abundances: ', as.character(max(gap)), '. 95th percentile under the null: ', as.character(null_max_gap_q95)))