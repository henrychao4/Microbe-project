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

discrete_gap = \(dat, eql_abuns, k_max, nboot) {
  numCores = detectCores()
  cl = makeCluster(numCores)
  
  clusterExport(cl, varlist = c("bootstrap_replicate", "I", "eql_abuns", "best_wKModes", "wKModes"))
  inputs = replicate(nboot, list(I = dat, eql_abuns = eql_abuns, k_max = k_max), simplify = FALSE)
  
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

w_2_vec = seq(0, .1, by = .01)

cont_p_val_vec = rep(0, length(w_2_vec))
cont_opt_k_vec = rep(0, length(w_2_vec))

disc_p_val_vec = rep(0, length(w_2_vec))
disc_opt_k_vec = rep(0, length(w_2_vec))

three_bin_p_val = rep(0, length(w_2_vec))

for (w in 1:length(w_2_vec)) {
  res_trait_1 = seq(0, (nres - 1) / nres, l = nres)
  res_trait_2 = seq(0, (nres - 1) / nres, l = nres)
  spec_trait_1 = seq(0, (nspec - 1) / nspec, l = nspec)
  spec_trait_2 = seq(0, (nspec - 1) / nspec, l = nspec)
  spec_trait_2 = sample(spec_trait_2, size = length(spec_trait_2), replace = F)
  dists_1 = circ_dist(spec_trait_1, res_trait_1)
  dists_2 = circ_dist(spec_trait_2, res_trait_2)
  w_1 = 1 - w_2_vec[w]
  w_2 = w_2_vec[w]
  C = exp(-((w_1 * dists_1^2) + (w_2 * dists_2^2)) / .01)
  
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
  
  I = (C > .3) * 1
  disc_gap = discrete_gap(dat = I, eql_abuns = eql_abuns, k_max = k_max, nboot = 50)
  
  cont_p_val_vec[w] = kmg_gap$p.value
  cont_opt_k_vec[w] = kmg_gap$khat
  
  disc_p_val_vec[w] = disc_gap$p_val
  disc_opt_k_vec[w] = disc_gap$opt_k
  
}

plot(w_2_vec, cont_p_val_vec, ylim = c(0,1), type = 'b', xlab = 'Noise from other trait axis', ylab = 'p-value for clustering', col = 'red')
lines(w_2_vec, disc_p_val_vec, type = 'b', col = 'blue')
abline(h = .05)
legend("topleft", legend = c("Continuous", "Discrete"), col = c("red", "blue"), lty = 1)

plot(w_2_vec, cont_opt_k_vec, ylim = c(0, k_max), type = 'b', xlab = 'Noise from other trait axis', ylab = 'Optimal k from gap statistic', col = 'red')
lines(w_2_vec, disc_opt_k_vec, type = 'b', col = 'blue')
legend("topleft", legend = c("Continuous", "Discrete"), col = c("red", "blue"), lty = 1)