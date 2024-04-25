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

nspec = 5
nres = 5

resources = seq(0, 1, length = nres)
traits = seq(0, 1, length = nspec)
dists = outer(resources, traits, FUN = \(x, y) abs(x - y))
dists[dists > .5] = 1 - dists[dists > .5]
C = exp(-(dists / .5) ^ 2)

params = list(
  nspec = nspec,
  nres = nres,
  immigration = rep(.1, nspec),
  death = rep(.1, nspec),
  inflow = rep(.1, nres),
  outflow = rep(.05, nres),
  assimilation = matrix(.8, nrow = nspec, ncol = nres),
  conversion = C,
  max_growth = matrix(.1, nrow = nspec, ncol = nres),
  half_saturation = matrix(.1, nrow = nspec, ncol = nres),
  byproduct = array(rep(.1, nspec * nres^2), dim=c(nspec, nres, nres))
)

init_state = c(rep(10, params$nspec), rep(5, params$nres))

mult = 
  \(params){
    ans = sapply(1:nrow(params$conversion), function(i) params$conversion[i,] %*% t(params$byproduct[i,,]))
    return(t(ans))
  }

index = 
  \(res){
    return(1)
  }

consumption = 
  \(params, R){
    ans = params$max_growth * (R / (params$half_saturation + R))
    return(ans)
  }

res_production = 
  \(params, R){
    ans = matrix(0,nrow(consumption(params, R)),ncol(consumption(params, R)))
    for (i in 1:nrow(consumption(params, R))) {
      ans[i,] = consumption(params, R)[i,] %*% params$byproduct[i,,]
    }
    return(ans)
  }

growth = 
  \(params, R){
    ans = rowSums((params$conversion * params$assimilation - mult(params)) * consumption(params, R))
    return(ans)
  }

model = 
  \(t, state, params){
    N = state[1:params$nspec]
    R = state[-(1:params$nspec)]
    dNdt = params$immigration + growth(params, R) * N - params$death * N
    dRdt = params$inflow - params$outflow * R + colSums(res_production(params, R)) - colSums(consumption(params, R) * N)
    return(list(c(dNdt,dRdt)))
  }

sim = ode(y = init_state, times = seq(0, 100, by = 1), func = model, parms = params)

sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 1))]
abuns.df = melt(spec.abuns, id.vars='time')
p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
print(p)
