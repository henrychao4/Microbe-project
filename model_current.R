library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

nspec = 2
nres = 2
init_spec = rep(10, nspec)
init_res = rep(5, nres)

immigration = .1
death = .1
inflow = .1
outflow = .1
assimilation = matrix(.1, nrow = nspec, ncol = nres) 
conversion = matrix(.1, nrow = nspec, ncol = nres)
max_growth = matrix(.1, nrow = nspec, ncol = nres)
half_saturation = matrix(.1, nrow = nspec, ncol = nres)
byproduct = array(rep(1, nspec * nres^2), dim=c(nspec, nres, nres))

mult = 
  \(conversion, byproduct){
    ans = sapply(1:nrow(conversion), function(i) conversion[i,] %*% t(byproduct[i,,]))
    return(t(ans))
  }

index = 
  \(res){
    return(1)
  }

consumption = 
  \(max_growth, half_saturation, res, index){
    ans = max_growth * (res / (half_saturation + res))
    return(ans)
  }

res_production = 
  \(consumption, byproduct){
    ans = matrix(0,nrow(consumption),ncol(consumption))
    for (i in 1:nrow(consumption)) {
      ans[i,] = consumption[i,] %*% byproduct[i,,]
    }
    return(ans)
  }

growth = 
  \(conversion, assimilation, byproduct, consumption){
    ans = rowSums((conversion * assimilation - mult(conversion, byproduct)) * consumption)
    return(ans)
  }