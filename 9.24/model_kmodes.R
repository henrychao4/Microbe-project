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

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

set.seed(1)

nspec = 30
nres = 10

# makeI = \(nspec, nres) {
#   I = matrix(0, nrow = nspec, ncol = nres)
#   cons_vec = c(rep(0, nres/2), rep(1, nres/2))
#   for (i in 1:nspec) {
#     I[i,] = sample(cons_vec, nres, replace = F)
#   }
#   return(I)
# }

makeI = \(nspec, nres) {
  I = matrix(rbernoulli(nspec * nres, p = .5), nrow = nspec, ncol = nres)
  return(I * 1)
}

library(purrr)

nclust = 5
n = 6
nspec = nclust * n
nres = 20

makeI = \(nclust, nspec, nres) {
  clust_traits = matrix(0, nrow = nclust, ncol = nres)
  nconsume = floor(nres/nclust)
  for (i in 1:nclust) {
    clust_traits[i, (nconsume * (i-1) + 1):(nconsume * i)] = rep(1, nconsume)
  }
  nspec = nclust * n
  I = matrix(0, nrow = nspec, nres)
  for (i in 1:nspec) {
    I[i,] = clust_traits[ceiling(i/n),]
    for (j in 1:nres) {
      flip = rbernoulli(1, .1)
      if (flip == 1) {
        if (I[i,j] == 0) {
          I[i,j] = 1
        }
        else {
          I[i,j] = 0
        }
      }
    }
  }
  return(I)
}

I = makeI(nclust, nspec, nres)

#I = makeI(nspec, nres)

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(0.01, nspec),
  mu = rep(.05, nspec),
  rho = rep(.3, nres),
  delta = rep(.05, nres),
  # mu = rep(0, nspec),
  # rho = rep(0, nres),
  # delta = rep(0, nres),
  epsilon = .9,
  beta = 1,
  h = 1,
  I = I,
  p = array(rep(0, nspec * nres^2), dim=c(nspec, nres, nres))
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

get_null_data = \(data) {
  n = nrow(data)
  m = ncol(data)
  null_data = matrix(0, nrow = n, ncol = m)
  
  for (j in 1:m) {
    col_j = data[,j]
    null_data[,j] = runif(n, min = min(col_j), max = max(col_j))
  }
  
  return(null_data)
}

best_wKModes = \(data, weights, modes, nruns) {
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

k_max = 10
pre_sim_gap_stat = clusGap(I, FUN = best_wKModes, K.max = k_max, B = 30, weights = rep(1, nspec), nruns = 10, verbose = F)
post_sim_gap_stat = clusGap(I, FUN = best_wKModes, K.max = k_max, B = 30, weights = eql_abuns, nruns = 10, verbose = F)
print(pre_sim_gap_stat)
print(post_sim_gap_stat)

plot(1:k_max, pre_sim_gap_stat$Tab[,1], type = 'b', col = 'red', xlab = 'k', ylab = 'logW', main = 'Initial abundances')
lines(1:k_max, pre_sim_gap_stat$Tab[,2], type = 'b', col = 'blue')
plot(1:k_max, pre_sim_gap_stat$Tab[,3], type = 'b', xlab = 'k', ylab = 'gap', main = 'Initial abundances')

plot(1:k_max, post_sim_gap_stat$Tab[,1], type = 'b', col = 'red', xlab = 'k', ylab = 'logW', main = 'At equilibrium')
lines(1:k_max, post_sim_gap_stat$Tab[,2], type = 'b', col = 'blue')
plot(1:k_max, post_sim_gap_stat$Tab[,3], type = 'b', xlab = 'k', ylab = 'gap', main = 'At equilibrium')

max_gap = max(post_sim_gap_stat$Tab[,3])

best_wKModes(I, weights = eql_abuns, modes = 5, nruns = 10)

dist(I, method = 'binary')

# n_bootstrap = 50
# boot_k = rep(0, n_bootstrap)
# boot_max_gaps = rep(0, n_bootstrap)
#
# for (i in 1:n_bootstrap) {
#   null_I = get_null_data(I)
#   clusgap = clusGap(null_I, FUN = wKModes, K.max = k_max, B = 20, weights = eql_abuns, verbose = F)
#   boot_max_gaps[i] = max(clusgap$Tab[,3])
# }
# 
# print(quantile(boot_max_gaps, .95))
# print(max_gap)
# 
# maxSE(post_sim_gap_stat$Tab[,3], post_sim_gap_stat$Tab[,4], method = 'Tibs2001SEmax')
# 
# best_wKModes(I, weights = eql_abuns, modes = 4, nruns = 10)
# clusGap(I, FUN = best_wKModes, K.max = k_max, B = 10, weights = eql_abuns, nruns = 10, verbose = F)