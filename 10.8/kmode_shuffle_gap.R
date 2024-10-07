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

set.seed(1)

nspec = 30
nres = 30

rotate_vector = function(v, shift) {
  n = length(v)
  if (shift == 0) {
    return(v)
  }
  shift = shift %% n
  v_rotated = c(v[(n-shift+1):n], v[1:(n-shift)])
  return(v_rotated)
}

# makeI = \(nspec, nres) {
#   vec = c(rep(1, nspec/nres), rep(0, nres - nspec/nres))
#   traits = rep(0, nspec)
#   I = matrix(0, nrow = nspec, ncol = nres)
#   for (i in 1:nspec) {
#     traits[i] = sample(0:(nres-1), 1)
#     I[i,] = rotate_vector(vec, traits[i])
#   }
#   return(list(I = I, traits = traits))
# }

makeI = \(nspec, nres) {
  width = 7
  vec = c(rep(1, width), rep(0, nres - width))
  traits = rep(0, nspec)
  I = matrix(0, nrow = nspec, ncol = nres)
  for (i in 1:nspec) {
    traits[i] = (i + floor(width/2)) %% nres
    I[i,] = rotate_vector(vec, (i-1))
  }
  traits[traits == 0] = nres
  return(list(I = I, traits = traits))
}

obj = makeI(nspec, nres)
I = obj$I
traits = obj$trait

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(0.05, nspec),
  mu = rep(.05, nspec),
  rho = rnorm(nres, mean = .6, sd = 0) + c(rep(0,9), .4, rep(0,9), .4, rep(0,9), .4),
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

plot(traits, eql_abuns, type = 'h')

k_max = 7
true_errs = rep(0, k_max)
null_errs = rep(0, k_max)
nboot = 30
null_gaps = matrix(0, nrow = nboot, ncol = k_max)
for (k in 1:k_max) {
  kmode = best_wKModes(I, modes = k, weights = eql_abuns, nruns = 10)
  true_errs[k] = sum(kmode$withindiff)
  replicate_err = rep(0, nboot)
  for (i in 1:nboot) {
    null_abuns = sample(eql_abuns, size = length(eql_abuns), replace = F)
    null_kmode = best_wKModes(I, modes = k, weights = null_abuns, nruns = 10)
    replicate_err[i] = sum(null_kmode$withindiff)
  }
  null_errs[k] = mean(replicate_err)
  
  for (i in 1:nboot) {
    null_gaps[,k] = null_errs[k] - replicate_err
  }
  
  print(k)
}

null_max_gaps = apply(null_gaps, 1, max)
null_max_gap_q95 = quantile(null_max_gaps,.95)

gap = null_errs - true_errs

plot(1:k_max, gap, type = 'b')

best_wKModes(I, modes = 3, weights = eql_abuns, nruns = 10)

print(paste0('True max gap for equilibrium abundances: ', as.character(max(gap)), '. 95th percentile under the null: ', as.character(null_max_gap_q95)))