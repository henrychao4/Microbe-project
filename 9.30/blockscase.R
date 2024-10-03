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

nspec = 120
nres = 20

rotate_vector = function(v, shift) {
  n = length(v)
  if (shift == 0) {
    return(v)
  }
  shift = shift %% n
  v_rotated = c(v[(shift+1):n], v[1:shift])
  return(v_rotated)
}

#come up with a way to make the block matrix
makeI = \(nspec, nres) {
  vec = c(rep(1, nspec/nres), rep(0, nres - nspec/nres))
  traits = rep(0, nspec)
  I = matrix(0, nrow = nspec, ncol = nres)
  for (i in 1:nspec) {
    traits[i] = sample(0:(nres-1), 1)
    I[i,] = rotate_vector(vec, traits[i])
  }
  return(list(I = I, traits = traits))
}

obj = makeI(nspec, nres)
I = obj$I
traits = obj$trait

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(0, nspec),
  mu = rep(.05, nspec),
  rho = rep(.7, nres),
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

get_null_data = \(data) {
  n = nrow(data)
  m = ncol(data)
  null_data = matrix(0, nrow = n, ncol = m)
  
  for (j in 1:m) {
    null_data[,j] = sample(data[,j], n, replace = T)
  }
  
  return(null_data)
}

best_kmodes = \(data, modes, nruns) {
  df = list()
  tot.withindiff_vec = rep(0, nruns)
  for (i in 1:nruns) {
    k_modes = kmodes(data, modes = modes)
    df[[i]] = k_modes
    tot.withindiff_vec[i] = sum(k_modes$withindiff)
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

dupe_data = 0
rounded_abuns = round(eql_abuns)
for (i in 1:nspec) {
  rep_spec = matrix(rep(I[i,], rounded_abuns[i]), nrow = rounded_abuns[i], ncol = ncol(I), byrow = T)
  dupe_data = rbind(dupe_data, rep_spec)
}
dupe_data = dupe_data[-1,]

k_max = 7
init_clusgap = clusGap(I, FUNcluster = kmodes, K.max = k_max, B = 10, verbose = F)
eql_clusgap = clusGap(dupe_data, FUNcluster = kmodes, K.max = k_max, B = 10)

plot(1:k_max, init_clusgap$Tab[,3], type = 'b')
plot(1:k_max, eql_clusgap$Tab[,3], type = 'b')

imax_gaps = rep(0, 5)
emax_gaps = rep(0, 5)

for (i in 1:5) {
  null_data = get_null_data(I)
  null_dupe_data = get_null_data(dupe_data)
  icg = clusGap(null_data, FUNcluster = kmodes, K.max = 5, B = 10, verbose = F)
  ecg = clusGap(null_dupe_data, FUNcluster = kmodes, K.max = k_max, B = 10, verbose = F)
  imax_gaps[i] = max(icg$Tab[,3])
  emax_gaps[i] = max(ecg$Tab[,3])
}

imax_gap = max(init_clusgap$Tab[,3])
emax_gap = max(eql_clusgap$Tab[,3])

iq95 = quantile(imax_gaps, .95)
eq95 = quantile(emax_gaps, .95)

print(paste0('True max gap for initial abundances: ', as.character(imax_gap), '. 95th percentile under the null: ', as.character(iq95)))

print(paste0('True max gap for equilibrium abundances: ', as.character(emax_gap), '. 95th percentile under the null: ', as.character(eq95)))