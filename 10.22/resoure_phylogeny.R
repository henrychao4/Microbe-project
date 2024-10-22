library(ggplot2)
library(deSolve)
library(tibble)
library(tidyr)
library(dplyr)
library(reshape2)
library(purrr)

set.seed(1)

nspec = 100
nres = 10
width = 5

P = matrix(0, nrow = nres, ncol = nres)
diag(P) = 9000
P[lower.tri(P)] = runif(((nres)^2 - nres) / 2, min = 0, max = 1)
#P[lower.tri(P)] = sample(x = c(.1, .9), size = ((nres)^2 - nres) / 2, replace = T, prob = c(.5, .5))
P[upper.tri(P)] = t(P)[upper.tri(P)]

spec_traits = sample.int(nres, nspec, replace = T)

C = matrix(0, nrow = nspec, ncol = nres)
I = matrix(0, nrow = nspec, ncol = nres)

for (i in 1:nspec) {
  C[i,] = P[spec_traits[i],]
  index = sample(1:nres, size = width, prob = C[i,], replace = FALSE)
  I[i,][index] = 1
}

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(0.01, nspec),
  mu = rep(.05, nspec),
  rho = rnorm(nres, mean = .6, sd = 0),
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
spec.abuns = sim.df[-((nspec + 2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
#print(p)

eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])
df = c(eql_abuns, num_coexist) |> t() |> as_tibble()
colnames(df) = c(paste0('N', seq(nspec)), 'num_coexist')

bootstrap_replicate = \(I = I, eql_abuns = eql_abuns, k_max){
  null_abuns = sample(eql_abuns, size = length(eql_abuns), replace = F)
  replicate_err = rep(0, k_max)
  for (k in 1:k_max) {
    null_kmode = best_wKModes(I, modes = k, weights = null_abuns, nruns = 10)
    replicate_err[k] = sum(null_kmode$withindiff)
  }
  return(replicate_err)
}

numCores = detectCores()
cl = makeCluster(numCores)

clusterExport(cl, varlist = c("bootstrap_replicate", "I", "eql_abuns", "best_wKModes", "wKModes"))

nboot = 50
k_max = 12
inputs = replicate(nboot, list(I = I, eql_abuns = eql_abuns, k_max = k_max), simplify = FALSE)

results = clusterApply(cl, inputs, function(input) {
  bootstrap_replicate(input$I, input$eql_abuns, input$k_max)
})

stopCluster(cl)

replicate_err = do.call(rbind, results)
null_errs = apply(replicate_err, 2, mean)

null_gaps = matrix(0, nrow = nboot, ncol = k_max)
for (k in 1:k_max) {
  null_gaps[,k] = null_errs[k] - replicate_err[,k]
}

true_errs = rep(0, k_max)
for (k in 1:k_max) {
  kmode = best_wKModes(I, modes = k, weights = eql_abuns, nruns = 10)
  true_errs[k] = sum(kmode$withindiff)
}

null_max_gaps = apply(null_gaps, 1, max)
null_k_opts = apply(null_gaps, 1, which.max)

null_max_gap_q95 = quantile(null_max_gaps, .95)

gap = null_errs - true_errs

plot(1:k_max, gap, type = 'b', xlab = 'k', ylab = 'Gap')

print(paste0('True max gap for equilibrium abundances: ', as.character(max(gap)), '. 95th percentile under the null: ', as.character(null_max_gap_q95)))

opt_k = which.max(gap)
null_opt_k_dist = null_max_gaps[null_k_opts == opt_k]
null_max_gap_opt_k_q95 = quantile(null_opt_k_dist, .95)

print(paste0('True max gap for equilibrium abundances: ', as.character(max(gap)), '. 95th percentile under the null for k = ', as.character(opt_k), ': ', as.character(null_max_gap_opt_k_q95)))

z_score = (max(gap) - mean(null_max_gaps)) / sd(null_max_gaps)
p_val = 1 - pnorm(z_score)

print(paste0('Z-score: ', as.character(z_score)))
print(paste0('P-value: ', as.character(p_val)))