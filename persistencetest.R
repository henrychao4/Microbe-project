library(ggplot2)
library(tensor)
library(deSolve)
library(tibble)
library(tidyr)
library(ggplot2)
library(dplyr)
library(reshape2)

m = .1
g = 1
I = 5
E = .2
nspec = 3
nres = 2
init_state = c(10,100,10,10,10)

resources = c(0,1)
traits = c(0,1,.5)
dists = outer(traits, resources, FUN = \(x, y) abs(x - y))
niche_width = .40
#niche_width = .42
C = exp(-(dists)^2 / niche_width)

MacArthur = 
  \(time, y, parms){
    R = with(parms, pmax(0, y[1:nres]))
    N = with(parms, pmax(0, y[(1 + nres):length(y)]))
    dRdt = with(parms, I - E * R - (t(C) %*% N) * R)
    dNdt = with(parms, N * (C %*% R - m))
    return(list(c(dRdt, dNdt)))
  }

sim = ode(y = init_state, times = seq(0, 1000), func = MacArthur, parms = list(m = m,  C = C, I = I, E = E, nspec = nspec, nres = nres))

sim.df = as.data.frame(sim)
abuns.df = sim.df[,c(1,4,5,6)]
abuns.plt.df = melt(abuns.df, id.vars='time')
p <- ggplot(abuns.plt.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
p + labs(title="Species Abundances Over Time") + theme(legend.position = "none")
print(p)