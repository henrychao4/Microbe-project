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

simulation = 
  \(

  ){
    
    sim = ode(y = init_state, times = seq(0, 100, by = 1), func = model, parms = params)
    
    return(num_alpha)
  }

params = 
  expand_grid(
    seed = 0,
    S = 5,
    m = .2,
    g = 1,
    r = 10,
    K = 10,
    niche_width = .3,
    R = 1:1000,
    N = c(100),
    epsilon_R = 0.01,
    epsilon_N = 0.01
  )

params = expand_grid(
  nspec = nspec,
  nres = nres,
  immigration = rep(.1, nspec),
  death = rep(.1, nspec),
  inflow = rep(.1, nres),
  outflow = rep(.05, nres),
  assimilation = matrix(.8, nrow = nspec, ncol = nres),
  conversion = matrix(.9, nrow = nspec, ncol = nres),
  max_growth = matrix(.1, nrow = nspec, ncol = nres),
  half_saturation = matrix(.1, nrow = nspec, ncol = nres),
  byproduct = array(rep(.1, nspec * nres^2), dim=c(nspec, nres, nres))
)
