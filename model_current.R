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

set.seed(2)

nspec = 20
nres = 20

makeI = \(nspec, nres) {
  I = matrix(0, nrow = nspec, ncol = nres)
  cons_vec = c(rep(0, nres/2), rep(1, nres/2))
  for (i in 1:nspec) {
    I[i,] = sample(cons_vec, nres, replace = F)
  }
  return(I)
}

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
  I = makeI(nspec, nres),
  p = array(rep(.01, nspec * nres^2), dim=c(nspec, nres, nres))
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

modularity = \(I) {
  
}

sim = ode(y = init_state, times = seq(0, 30000, by = 1), func = model, parms = params)

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

dists = matrix(0, nrow = nspec, ncol = nspec)
for (i in 1:nspec) {
  for (j in 1:nspec) {
    dists[i,j] = hamming_dist(params$I[i,], params$I[j,])
  }
}

num_alpha = find_num_alpha(eql, eql_abuns, model, params)
coexist_traits = params$I[eql_abuns > .1,]

overall_dists = matrix(0, nrow = nspec, ncol = nspec)
coex_dists = matrix(0, nrow = nrow(coexist_traits), ncol = nrow(coexist_traits))

for (i in 1:nspec) {
  for (j in 1:nspec) {
    overall_dists[i,j] = hamming_dist(params$I[i,], params$I[j,])
  }
}

for (i in 1:nrow(coexist_traits)) {
  for (j in 1:nrow(coexist_traits)) {
    coex_dists[i,j] = hamming_dist(coexist_traits[i,], coexist_traits[j,])
  }
}

overall_avg_ham_dist = mean(dists)
coexist_avg_ham_dist = mean(coex_dists)
