library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(reshape2)
library(lsa)
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

nspec = 30
nres = 30

ntraits = 2

res_traits = matrix(runif(nres * ntraits, min = 0, max = 1), nrow = nres, ncol = ntraits)
res_dists = as.matrix(dist(res_traits, upper = T, diag = T))
res_sim = exp(- res_dists^2 / .2)
C = res_sim

# res_traits = matrix(runif(nres * ntraits, min = 0, max = 1), nrow = ntraits, ncol = nres)
# css = cosine(res_traits)

# C = css
# C = C^10

# row_totals = rowSums(C)
# min_row_tot = min(row_totals)
# 
# for (j in 1:nrow(C)) {
#   C[j,] = C[j,] * (min_row_tot / sum(C[j,]))
# }
# 
# C = C / rowSums(C) * 5

params = list(
  nspec = nspec,
  nres = nres,
  alpha = 0.05,
  r = 50,
  K = 10,
  m = .2,
  C = C
)

init_abuns = rep(5, params$nspec)+ rnorm(params$nspec, mean = 0, sd = 0.01)
init_res = rep(5, params$nres)+ rnorm(params$nspec, mean = 0, sd = 0.01)
init_state = c(init_abuns, init_res)

sim = ode(y = init_state, times = seq(0, 1000, by = .01), func = MacArthur, parms = params, method = "rk4")
sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]

p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')
#print(p)

plot(eql_abuns, type = 'h')
#plot(spec_trait_1, eql_abuns, type = 'h')
# plot(spec_trait_2, eql_abuns, type = 'h')

kmg_data = as.data.frame(C)
kmg_data$N = round(eql_abuns)
kmg_gap = KmeansGap(dat = kmg_data, multiD = T, mink = 1, maxk = 10)

plot(kmg_gap$data$k, kmg_gap$data$gap, type = 'b', main = paste0('Number of resource dimensions = ', as.character(ntraits)))

print(kmg_gap)