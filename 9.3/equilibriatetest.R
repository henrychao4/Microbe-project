library(ggplot2)
library(tensor)
library(future)
library(parallel)
library(deSolve)
library(tibble)
library(tidyr)
library(dplyr)
library(reshape2)
library(rootSolve)

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

set.seed(1)

nspec = 20
nres = 10 #number of traits per species

makeOverlap = \(I) {
  overlap = matrix(0, nrow = nrow(I), ncol = ncol(I))
  for (i in 1:nrow(I)) {
    for (j in 1:ncol(I)) {
      if (I[i,j] == 1) {
        overlap[i,j] = sum(I[,j] == I[i,j]) - 1
      }
    }
  }
  return(overlap)
}

makeI = \(nspec, nres) {
  I = matrix(0, nrow = nspec, ncol = nres)
  cons_vec = c(rep(1, nres/2), rep(0, nres/2))
  for (i in 1:nspec) {
    I[i,] = sample(cons_vec, nres, replace = F)
  }
  while (sum(duplicated(I)) > 0) {
    for (i in 1:nspec) {
      I[i,] = sample(cons_vec, nres, replace = F)
    }
  }
  return(I)
}

I = makeI(nspec, nres)

overlap = makeOverlap(I)
total_overlap = rowSums(overlap)

min_rm_zero = \(x) {
  y = x[x != 0]
  return(min(y))
}

min_overlap = apply(overlap, 1, FUN = min_rm_zero)

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(0, nspec),
  mu = rep(.05, nspec),
  rho = rep(.7, nres),
  delta = rep(.05, nres),
  epsilon = 1,
  beta = 1,
  h = 1,
  I = I,
  p = array(rep(0, nspec * nres^2), dim=c(nspec, nres, nres))
)

init_state = c(rep(1, params$nspec), rep(1, params$nres))

U = 
  \(params, R){
    ans = t(params$beta * R * t(params$I))
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
    N[N < 0] = 0
    dNdt = params$alpha + gamma(params, R) * N - params$mu * N
    dRdt = params$rho - params$delta * R + colSums(omega(params, R) * N) - colSums(U(params, R) * N)
    
    if (19999 < t && t < 20000) {
      dNdt = dNdt + 5
    }
    
    return(list(c(dNdt, dRdt)))
  }

sim = ode(y = init_state, times = seq(0, 30000), func = model, parms = params)

sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()

pre_pert_eql = spec.abuns[spec.abuns$time == 19000,][-1]

eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])

plot(total_overlap, eql_abuns)
plot(min_overlap, eql_abuns)
print(p)