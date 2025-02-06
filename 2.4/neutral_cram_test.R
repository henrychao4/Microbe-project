library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(reshape2)
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

ar_sample = \(res_c, res_m) {
  cond = F
  while (cond == F) {
    sample_c1 = runif(1)
    sample_c2 = runif(1)
    sample_d = runif(1)
    req_d = sample_c1 * ((res_c[2,2] * res_m[1] - res_c[1,2] * res_m[2]) / (res_c[1,1]* res_c[2,2] - res_c[1,2] * res_c[2,1]))
    + sample_c2 * ((res_c[1,1] * res_m[2] - res_c[2,1] * res_m[1]) / (res_c[1,1]* res_c[2,2] - res_c[1,2] * res_c[2,1]))
    
    if (sample_d > req_d) {
      cond = T
      return(c(sample_c1, sample_c2, sample_d))
    }
  }
}

nspec = 100
nres = 2

res_c = matrix(c(.5, .2, .2, .5), nrow = 2, ncol = 2)
res_m = c(.2, .2)

invaders = matrix(0, nrow = nspec - 2, ncol = nres + 1)

for (i in 1:nspec - 2) {
  invaders[i,] = ar_sample(res_c, res_m)
}

C = matrix(0, nrow = nspec, ncol = nres)
C[1:2,] = res_c
C[3:nspec,] = invaders[,1:2]

m = rep(0, nspec)
m[1:2] = res_m
m[3:nspec] = invaders[,3]

params = list(
  nspec = nspec,
  nres = nres,
  alpha = 0.1,
  r = 3 * nspec,
  K = 1,
  m = m,
  C = C
)

init_abuns = rep(5, params$nspec)
init_res = rep(5, params$nres)
init_state = c(init_abuns, init_res)

sim = ode(y = init_state, times = seq(0, 5000), func = MacArthur, parms = params)
sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])

p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')
#print(p)

plot(C[,1], eql_abuns, type = 'h', xlab = 'Species trait 1', ylab = 'Equilibrium Abundance')

kmg_data = as.data.frame(C)
kmg_data$N = round(eql_abuns)
kmg_gap = KmeansGap(dat = kmg_data, multiD = T, mink = 1, maxk = 10)

plot(kmg_gap$data$k, kmg_gap$data$gap, type = 'b', xlab = 'k', ylab = 'Gap')

print(kmg_gap)