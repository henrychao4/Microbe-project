library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
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

    dNdt = with(parms, N * ((C %*% R) - m))
    dRdt = with(parms, R * (r * (1 - R / K) - t(C) %*% N))
    return(list(c(dNdt, dRdt)))
  }

nspec = 10
nres = 5
C = matrix(0, nrow = nspec, ncol = nres)
vec = rnorm(nres, mean = .7, sd = .1)
vec_2 = rnorm(nres, mean = 0, sd = .05)
#vec = c(rep(1, nres/2), rep(0, nres/2))
for (i in 1:nspec) {
  C[i,] = sample(vec, nres, replace = F) + sample(vec_2, nres, replace = F)
}


params = list(
  nspec = nspec,
  nres = nres,
  r = 50,
  K = 1,
  m = .2,
  C = C
  )

init_state = c(rep(1, params$nspec), rep(1, params$nres))

sim = ode(y = init_state, times = seq(0, 50000), func = MacArthur, parms = params)
sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])

p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')
print(p)