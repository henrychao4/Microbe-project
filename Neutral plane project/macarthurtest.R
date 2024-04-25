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

nspec = 10
nclust = nspec / 2
variation = rnorm(nspec, 1, .5)

params = list(
  nspec = nspec,
  nclust = nclust,
  inflow_bad = 1,
  inflow_good = 10,
  outflow = .1,
  death = c(rep(.1, nclust), rep(.05, nclust)) * variation,
  C = c(rep(.23, nclust), rep(.1, nclust)) * variation,
  K = c(rep(1, nclust), rep(.7, nclust))
)

init_state = c(20, rep(10, params$nspec))

growth =
  \(R, params){
    ans = params$C * (R / (params$K + R))
    return(ans)
  }

model = 
  \(t, state, params){
    R = state[1]
    N = state[-1]
    if ((t %% 200) < 100) {
      dRdt = params$inflow_bad - params$outflow * R - sum(growth(R,params) * N)
    }
    else {
      dRdt = params$inflow_good - params$outflow * R - sum(growth(R,params) * N)
    }
    dNdt = N * (growth(R,params) - params$death)
    return(list(c(dRdt,dNdt)))
  }

sim = ode(y = init_state, times = seq(0, 100000, by = 1), func = model, parms = params)

sim.df = as.data.frame(sim)
spec.abuns = sim.df[-2]
abuns.df = melt(spec.abuns, id.vars='time')
p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
print(p)