library(ggplot2)
library(tensor)
library(future)
library(parallel)
library(deSolve)
library(tibble)
library(tidyr)
library(dplyr)
library(reshape2)
library(cluster)
library(purrr)

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

set.seed(1)

nclust = 5 #number of clusters
n = 4 #number of species per cluster
nres = 10 #number of traits per species

clust_traits = matrix(rbinom(nclust * nres, 1, .5), nrow = nclust, ncol = nres)

nspec = nclust * n
I = matrix(0, nrow = nspec, nres)

makeI = \(n, nspec, nres) {
  for (i in 1:nspec) {
    I[i,] = clust_traits[ceiling(i/n),]
    for (j in 1:nres) {
      flip = rbernoulli(1, .1)
      if (flip == 1) {
        if (I[i,j] == 0) {
          I[i,j] = 1
        }
        if (I[i,j] == 1) {
          I[i,j] = 0
        }
      }
    }
  }
  return(I)
}

I = makeI(n, nspec, nres)

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(0.01, nspec),
  mu = rep(.05, nspec),
  rho = rep(.3, nres),
  delta = rep(.05, nres),
  epsilon = .9,
  beta = 1,
  h = 1,
  I = I,
  p = array(rep(0, nspec * nres^2), dim=c(nspec, nres, nres))
)

init_state = c(rep(1, params$nspec), rep(1, params$nres), 0)

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

sim = ode(y = init_state, times = seq(0, 30000, by = 1), func = model, parms = params)

sim.df = as.data.frame(sim)
spec.abuns = sim.df[-((nspec+2):(nspec + nres + 2))]
abuns.df = melt(spec.abuns, id.vars='time')
p <- ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
print(p)

eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])
df = c(eql_abuns, num_coexist) |> t() |> as_tibble()
colnames(df) = c(paste0('N', seq(nspec)), 'num_coexist')

weights = matrix(0, nrow = nspec, ncol = nspec)
for (i in 1:nspec) {
  for (j in 1:nspec) {
    weights[i,j] = eql_abuns[i] * eql_abuns[j]
  }
}

exps = seq(0, 1, .05)
optimal_clusters = rep(0, length(exps))

for (i in 1:length(exps)) {
  coex_weight_d = as.matrix(dist(params$I, method = "binary")) / weights^(exps[i])
  coex_weight_hc = hclust(as.dist(coex_weight_d))
  
  calculate_silhouette = function(d, hc, k) {
    clusters = cutree(hc, k = k)
    sil_width = silhouette(clusters, dist = d)
    avg_silhouette = mean(sil_width[, "sil_width"])
    return(avg_silhouette)
  }
  
  max_clusters = 5
  silhouette_scores = rep(-Inf, max_clusters)
  for (k in 2:max_clusters) {
    silhouette_scores[k] = calculate_silhouette(coex_weight_d, coex_weight_hc, k)
  }
  
  optimal_clusters[i] = which.max(silhouette_scores)
}