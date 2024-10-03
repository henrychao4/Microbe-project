library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
source("KmeansGap.r")

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

set.seed(1)

MacArthur = 
  \(time, state, parms){
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

best_kmeans = \(data, k, nruns) {
  df = list()
  tot.withinss_vec = rep(0, nruns)
  for (i in 1:nruns) {
    k_means = kmeans(data, k)
    df[[i]] = k_means
    tot.withinss_vec[i] = k_means$tot.withinss
  }
  idx = which.min(tot.withinss_vec)
  return(df[[idx]])
}

get_null_data = \(data) {
  n = nrow(data)
  m = ncol(data)
  null_data = matrix(0, nrow = n, ncol = m)
  
  for (j in 1:m) {
    null_data[,j] = sample(data[,j], n, replace = T)
    null_data[,j] = runif(n, min = min(data[,j]), max = max(data[,j]))
  }
  
  return(null_data)
}

get_null_shuffled_data = \(I, abuns) {
  shuffled_data = 0
  rounded_abuns = round(abuns)
  for (i in 1:nspec) {
    rep_spec = matrix(rep(I[i,], rounded_abuns[i]), nrow = rounded_abuns[i], ncol = ncol(I), byrow = T)
    shuffled_data = rbind(shuffled_data, rep_spec)
  }
  shuffled_data = shuffled_data[-1,]
  return(shuffled_data)
}

nspec = 20
nres = 3

res_trait = seq(.2, .8, l = nres)
spec_trait = seq(0, .9, l = nspec) + rnorm(nspec, mean = 0, sd = .005)
dists = circ_dist(spec_trait, res_trait)
C = exp(-dists^2 / .05)

params = list(
  nspec = nspec,
  nres = nres,
  alpha = .05,
  r = 50,
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

plot(spec_trait, eql_abuns, type = 'h')

eql_data = 0
rounded_eql_abuns = round(eql_abuns)
for (i in 1:nspec) {
  rep_spec = matrix(rep(C[i,], rounded_eql_abuns[i]), nrow = rounded_eql_abuns[i], ncol = ncol(C), byrow = T)
  eql_data = rbind(eql_data, rep_spec)
}
eql_data = eql_data[-1,]

k_max = 10
eql_clusgap = clusGap(eql_data, FUNcluster = kmodes, K.max = k_max, B = 5)

plot(1:k_max, eql_clusgap$Tab[,3], type = 'b', xlab = 'k', ylab = 'Gap', main = 'Equilibrium Abundances')

emax_gaps = rep(0, 5)

for (i in 1:5) {
  null_C = get_null_data(C)
  null_eql_data = get_null_shuffled_data(null_C, eql_abuns)
  ecg = clusGap(null_eql_data, FUNcluster = kmodes, K.max = k_max, B = 5, verbose = F)
  emax_gaps[i] = max(ecg$Tab[,3])
  print(i)
}

emax_gap = max(eql_clusgap$Tab[,3])

eq95 = quantile(emax_gaps, .95)

print(paste0('True max gap for equilibrium abundances: ', as.character(emax_gap), '. 95th percentile under the null: ', as.character(eq95)))

kmg_data = as.data.frame(C)
kmg_data$N = round(eql_abuns)
kmg_gap = KmeansGap(dat = kmg_data, multiD = T, mink = 1, maxk = 19)

plot(kmg_gap$data$k, kmg_gap$data$gap, type = 'b')