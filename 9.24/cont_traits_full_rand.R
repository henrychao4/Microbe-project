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

nspec = 20
nres = 3

C = matrix(rnorm(n = nspec * nres, mean = .5, sd = .1), nrow = nspec, ncol = nres)

params = list(
  nspec = nspec,
  nres = nres,
  alpha = .05,
  r = 50,
  K = 1,
  m = .2,
  C = C
)

init_state = c(rep(1, params$nspec), rep(1, params$nres))

sim = ode(y = init_state, times = seq(0, 15000), func = MacArthur, parms = params)
sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])

p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')
#print(p)

init_data = as.data.frame(C)
init_data$N = rep(1, nspec)
#init_data$N = sample(c(1,2,3), nspec, replace = T)

eql_data = as.data.frame(C)
eql_data$N = round(eql_abuns)

#because the method is permutation shuffling, gap curve is flat
init_gap = KmeansGap(dat = init_data, multiD = T, mink = 1, maxk = 9)

eql_gap = KmeansGap(dat = eql_data, multiD = T, mink = 1, maxk = 9)

plot(init_gap$data$k, init_gap$data$gap, type = 'b')

plot(eql_gap$data$k, eql_gap$data$gap, type = 'b')

dupe_data = 0
rounded_abuns = round(eql_abuns)
for (i in 1:nspec) {
  rep_spec = matrix(rep(C[i,], rounded_abuns[i]), nrow = rounded_abuns[i], ncol = ncol(C), byrow = T)
  dupe_data = rbind(dupe_data, rep_spec)
}
dupe_data = dupe_data[-1,]

init_clusgap = clusGap(C, FUNcluster = kmeans, K.max = 19, B = 20)
eql_clusgap = clusGap(dupe_data, FUNcluster = best_kmeans, K.max = 19, B = 20, nruns = 10)

plot(1:19, init_clusgap$Tab[,3], type = 'b')
plot(1:19, eql_clusgap$Tab[,3], type = 'b')
