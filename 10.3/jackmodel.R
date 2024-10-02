library(deSolve)

set.seed(1)

model = 
  \(time, state, parms){
    N_AB = 0
    N_Ab = 0
    N_aB = 0
    N_ab = 0
    R = 0
    dNdt = with(parms, alpha + N * ((C %*% R) - m))
    dRdt = with(parms, R * (r * (1 - R / K) - t(C) %*% N))
    return(list(c(dNdt, dRdt)))
  }

props = c('AB' = .4, 'Ab' = .1, 'aB' = .25, 'ab' = .25)

new_props = \(props) {
  propA = props['AB'] + props['Ab']
  propa = 1 - propA
  propB = props['AB'] + props['aB']
  propb = 1 - propB
  
  newAB = as.numeric(propA * propB)
  newAb = as.numeric(propA * propb)
  newaB = as.numeric(propa * propB)
  newab = as.numeric(propa * propb)
  
  return(c('AB' = newAB, 'Ab' = newAb, 'aB' = newaB, 'ab' = newab))
}

init_abuns = c(1,2,3)
init_res = 5
init_state = c(init_abuns, init_res)


overall_growth = \(gen_abuns, res_abuns) {
  
}

a_1 = diag(3)

rho = c(90, 1, 9)
m = c(0, 0, 0)



init_state = c(rep(1, params$nspec), rep(1, params$nres))

sim = ode(y = init_state, times = seq(0, 15000), func = MacArthur, parms = params)
sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])

p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic() + ggtitle('MacArthur')
print(p)