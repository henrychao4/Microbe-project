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

nspec = 2

params = list(
  nspec = nspec,
  inflow = 1,
  outflow = .1,
  death = c(.17, .1),
  mu_max = c(.2, .1),
  q_min = c(1, 1),
  v_max = c(.4, .4),
  K = c(10, 10)
)

init_state = c(20, rep(1, params$nspec), rep(10, params$nspec))

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
    dQdt = uptake(R, params) - growth(Q, params)
    dNdt = N * (growth(Q, params) - params$death)
    # print(R)
    # print(Q)
    print(dNdt)
    return(list(c(dRdt,dQdt,dNdt)))
  }

sim = ode(y = init_state, times = seq(0, 1000, by = 1), func = model, parms = params)

sim.df = as.data.frame(sim)
spec.abuns = sim.df[-(2:4)]
abuns.df = melt(spec.abuns, id.vars='time')
p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
print(p)