library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(reshape2)
library(MEDseq)
source("KmeansGap.r")

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

set.seed(1)

MacArthur = \(time, state, parms){
  N = state[1:params$nspec]
  R = state[(params$nspec + 1):(params$nspec + params$nres)]
  dNdt = with(parms, alpha + N * ((C %*% R) - m))
  dRdt = with(parms, R * (r * (1 - R / K) - t(C) %*% N))
  return(list(c(dNdt, dRdt)))
}

circ_dist = function(vec1, vec2) {
  y = abs(outer(vec1, vec2, '-'))
  idx = which(y > 0.5)
  y[idx] = 1 - y[idx]
  return(y)
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

bootstrap_replicate = \(I = I, eql_abuns = eql_abuns, k_max){
  null_abuns = sample(eql_abuns, size = length(eql_abuns), replace = F)
  replicate_err = rep(0, k_max)
  for (k in 1:k_max) {
    null_kmode = best_wKModes(I, modes = k, weights = null_abuns, nruns = 10)
    replicate_err[k] = null_kmode$tot.withindiff
  }
  return(replicate_err)
}

discrete_gap = \(dat, k_max, nboot) {
  numCores = detectCores()
  cl = makeCluster(numCores)
  
  clusterExport(cl, varlist = c("bootstrap_replicate", "I", "eql_abuns", "best_wKModes", "wKModes"))
  inputs = replicate(nboot, list(I = I, eql_abuns = eql_abuns, k_max = k_max), simplify = FALSE)
  
  results = clusterApply(cl, inputs, function(input) {
    bootstrap_replicate(input$I, input$eql_abuns, input$k_max)
  })
  
  stopCluster(cl)
  
  replicate_err = do.call(rbind, results)
  replicate_err = log(replicate_err)
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
  true_errs = log(true_errs)
  
  null_max_gaps = apply(null_gaps, 1, max)
  null_k_opts = apply(null_gaps, 1, which.max)
  
  null_max_gap_q95 = quantile(null_max_gaps, .95)
  
  gap = null_errs - true_errs
  
  opt_k = which.max(gap)
  null_opt_k_dist = null_max_gaps[null_k_opts == opt_k]
  null_max_gap_opt_k_q95 = quantile(null_opt_k_dist, .95)
  
  z_score = (max(gap) - mean(null_max_gaps)) / sd(null_max_gaps)
  p_val = 1 - pnorm(z_score)
  
  return(list(k = 1:k_max, gap = gap, opt_k = opt_k, null_max_gap_q95 = null_max_gap_q95, z_score = z_score, p_val = p_val))
}

nspec = 240
nres = 6

res_trait_1 = seq(0, (nres - 1) / nres, l = nres)
res_trait_2 = seq(0, (nres - 1) / nres, l = nres)
spec_trait_1 = seq(0, (nspec - 1) / nspec, l = nspec)
#spec_trait_1 = sample(spec_trait_1, size = length(spec_trait_1), replace = F)
spec_trait_2 = seq(0, (nspec - 1) / nspec, l = nspec)
spec_trait_2 = sample(spec_trait_2, size = length(spec_trait_2), replace = F)
dists_1 = circ_dist(spec_trait_1, res_trait_1)
#dists_2 = 0
dists_2 = circ_dist(spec_trait_2, res_trait_2)
w_1 = .99
w_2 = 0.01
C = exp(- ((w_1 * dists_1^2) + (w_2 * dists_2^2)) / .01)
#C = exp(-pmax((dists_1^2 / sigma_1), (dists_2^2 / sigma_2)))

params = list(
  nspec = nspec,
  nres = nres,
  alpha = .05,
  r = 500,
  K = 1,
  m = .2,
  C = C
)

init_abuns = rep(5, params$nspec)
init_res = rep(5, params$nres)
init_state = c(init_abuns, init_res)

sim = ode(y = init_state, times = seq(0, 15000), func = MacArthur, parms = params)
sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])

p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')
#print(p)

plot(spec_trait_1, eql_abuns, type = 'h')

k_max = 8

kmg_data = as.data.frame(C)
kmg_data$N = round(eql_abuns)
kmg_gap = KmeansGap(dat = kmg_data, multiD = T, mink = 1, maxk = k_max)

print(kmg_gap)
plot(kmg_gap$data$k, kmg_gap$data$gap, type = 'b', xlab = 'k', ylab = 'Gap', main = 'Continuous Gap')

# I = matrix(0, nrow = nspec, ncol = nres)
# I[(C > 0) & (C < .2)] = 0
# I[(C > .2) & (C < .4)] = 1
# I[(C > .4) & (C < .6)] = 2
# I[(C > .6) & (C < .8)] = 3
# I[(C > .8) & (C < 1)] = 4

I = (C > .3) * 1

disc_gap = discrete_gap(dat = I, k_max = k_max, nboot = 50)

print(disc_gap)

plot(disc_gap$k, disc_gap$gap, type = 'b', xlab = 'k', ylab = 'Gap', main = 'Discrete Gap')