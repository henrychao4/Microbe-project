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

set.seed(1)

nspec = 3
C = seq(0,1, length.out = nspec)

params = list(
  nspec = nspec,
  inflow = 1,
  outflow = .1,
  C = C,
  death = .05 * C
)

init_state = c(10, rep(10, params$nspec))

model = 
  \(t, state, params){
    R = state[1]
    N = state[-1]
    dRdt = params$inflow - params$outflow * R - sum(params$C * R * N)
    dNdt = N * (params$C * R - params$death)
    return(list(c(dRdt,dNdt)))
  }

sim = ode(y = init_state, times = seq(0, 1000, by = 1), func = model, parms = params)

sim.df = as.data.frame(sim)
spec.abuns = sim.df[-2]
abuns.df = melt(spec.abuns, id.vars='time')
p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
print(p)