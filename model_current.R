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
  \(conversion, assimilation, byproduct){
    ans = conversion * assimilation - conversion
    return(1)
  }

growth = 
  \(conversion, assimilation, byproduct){
    return(1)
  }

res_production = 
  \(conversion, assimilation, byproduct){
    return(1)
  }