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

nspec = 4
nres = 4

I = matrix(c(1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0), nrow = 4, ncol = 4)

overlap = makeOverlap(I)
total_overlap = rowSums(overlap)

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(0, nspec),
  mu = rep(.05, nspec),
  rho = rep(.3, nres),
  delta = rep(.05, nres),
  # mu = rep(0, nspec),
  # rho = rep(0, nres),
  # delta = rep(0, nres),
  epsilon = 1,
  beta = 1,
  h = 1,
  I = I,
  p = array(rep(0, nspec * nres^2), dim=c(nspec, nres, nres))
)

init_state = c(rep(1, params$nspec), rep(1, params$nres))

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
    dNdt = params$alpha + gamma(params, R) * N - params$mu * N
    dRdt = params$rho - params$delta * R + colSums(omega(params, R) * N) - colSums(U(params, R) * N)
    
    return(list(c(dNdt, dRdt)))
  }

sim = ode(y = init_state, times = seq(0, 10000, by = 1), func = model, parms = params)

sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
print(p)

eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])