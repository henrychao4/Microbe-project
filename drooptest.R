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

nspec = 6

params = list(
  nspec = nspec,
  inflow = 5,
  outflow = .1,
  death = .1,
  mu_max = rep(.1, nspec),
  q_min = rep(.1, nspec),
  v_max = rep(.1, nspec),
  K = rep(.1, nspec)
)

init_state = c(20, rep(10, params$nspec), rep(0, params$nspec))

growth = 
  \(Q, params){
    ans = params$mu_max * (1 - params$q_min / Q)
    return(ans)
  }

uptake =
  \(R, params){
    ans = params$v_max * (R / (params$K + R))
    return(ans)
  }

model = 
  \(t, state, params){
    R = state[1]
    Q = state[2:(params$nspec + 1)]
    N = state[-(1:(params$nspec + 1))]
    dRdt = params$inflow - params$outflow * R - sum(uptake(R, params) * N)
    dQdt = uptake(R, params) * N - growth(Q, params)
    dNdt = growth(Q, params) - params$death * N
    return(list(c(dRdt,dQdt,dNdt)))
  }

sim = ode(y = init_state, times = seq(0, 100, by = 1), func = model, parms = params)

sim.df = as.data.frame(sim)
abuns.df = melt(sim.df, id.vars='time')
p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
print(p)