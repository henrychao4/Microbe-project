library(ggplot2)
library(deSolve)
library(klaR)

set.seed(3)

nspec = 20
nres = 20

makeI = \(nclust, n, nres) {
  clust_traits = matrix(0, nrow = nclust, ncol = nres)
  nconsume = floor(nres/nclust)
  for (i in 1:nclust) {
    clust_traits[i, (nconsume * (i-1) + 1):(nconsume * i)] = rep(1, nconsume)
  }
  nspec = nclust * n
  I = matrix(0, nrow = nspec, nres)
  for (i in 1:nspec) {
    I[i,] = clust_traits[ceiling(i/n),]
    for (j in 1:nres) {
      flip = rbernoulli(1, .1)
      if (flip == 1) {
        if (I[i,j] == 0) {
          I[i,j] = 1
        }
        else {
          I[i,j] = 0
        }
      }
    }
  }
  return(I)
}

I = makeI(5, 4, 20)

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(.05, nspec),
  mu = rep(.05, nspec),
  rho = rep(.5, nres),
  delta = rep(.05, nres),
  epsilon = .9,
  beta = 1,
  h = 1,
  I = I,
  p = array(rep(0, nspec * nres^2), dim = c(nspec, nres, nres))
)

init_abuns = rep(1, params$nspec)
init_res = rep(1, params$nres)
init_state = c(init_abuns, init_res, 0)

U = 
  \(params, R){
    ans = t(params$beta * R / (params$h + R) * t(params$I))
    return(ans)
  }

omega = 
  \(params, R){
    ans = rowSums(params$p, dims = 2) * U(params, R)
    return(ans)
  }

gamma = 
  \(params, R){
    ans = params$epsilon * rowSums((1 - rowSums(params$p, dims = 2)) * U(params, R))
    return(ans)
  }

model = 
  \(t, state, params){
    N = state[1:params$nspec]
    R = state[(params$nspec + 1):(params$nspec + params$nres)]
    W = tail(state, 1)
    dNdt = params$alpha + gamma(params, R) * N - params$mu * N
    dRdt = params$rho - params$delta * R + colSums(omega(params, R) * N) - colSums(U(params, R) * N)
    dWdt = sum((1 - params$epsilon) * (1 - rowSums(params$p, dims = 2)) * U(params, R) * N)
    
    return(list(c(dNdt, dRdt, dWdt)))
  }

get_null_data = \(data) {
  n = nrow(data)
  m = ncol(data)
  null_data = matrix(0, nrow = n, ncol = m)
  
  for (j in 1:m) {
    null_data[,j] = sample(data[,j], n, replace = T)
  }
  
  return(null_data)
}

best_kmodes = \(data, modes, nruns) {
  df = list()
  tot.withindiff_vec = rep(0, nruns)
  for (i in 1:nruns) {
    k_modes = kmodes(data, modes = modes)
    df[[i]] = k_modes
    tot.withindiff_vec[i] = sum(k_modes$withindiff)
  }
  idx = which.min(tot.withindiff_vec)
  return(df[[idx]])
}

sim = ode(y = init_state, times = seq(0, 30000, by = 1), func = model, parms = params)

sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
#print(p)

eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])
df = c(eql_abuns, num_coexist) |> t() |> as_tibble()
colnames(df) = c(paste0('N', seq(nspec)), 'num_coexist')

init_data = 0
rounded_init_abuns = round(init_abuns)
for (i in 1:nspec) {
  rep_spec = matrix(rep(I[i,], rounded_init_abuns[i]), nrow = rounded_init_abuns[i], ncol = ncol(I), byrow = T)
  init_data = rbind(init_data, rep_spec)
}
init_data = init_data[-1,]

eql_data = 0
rounded_eql_abuns = round(eql_abuns)
for (i in 1:nspec) {
  rep_spec = matrix(rep(I[i,], rounded_eql_abuns[i]), nrow = rounded_eql_abuns[i], ncol = ncol(I), byrow = T)
  eql_data = rbind(eql_data, rep_spec)
}
eql_data = eql_data[-1,]

k_max = 15
init_clusgap = clusGap(init_data, FUNcluster = kmodes, K.max = k_max, B = 10, verbose = F)
eql_clusgap = clusGap(eql_data, FUNcluster = kmodes, K.max = k_max, B = 10)

plot(1:k_max, init_clusgap$Tab[,3], type = 'b', xlab = 'k', ylab = 'Gap', main = 'Initial Abundances of 1 each')
plot(1:k_max, eql_clusgap$Tab[,3], type = 'b', xlab = 'k', ylab = 'Gap', main = 'Equilibrium Abundances')