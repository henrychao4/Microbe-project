library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(bipartite)

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

nspec = 10
nres = 10

init_state = c(rep(1, nspec), rep(1, nres), 0)

makeI = \(nspec, nres) {
  I = matrix(0, nrow = nspec, ncol = nres)
  cons_vec = c(rep(0, nres/2), rep(1, nres/2))
  for (i in 1:nspec) {
    I[i,] = sample(cons_vec, nres, replace = F)
  }
  return(I)
}

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

hamming_dist = \(x,y) {
  ans = sum(x != y) / length(x)
  return(ans)
}

simulation = 
  \(
    numsim,
    nspec,
    nres
  ){
    set.seed(numsim)
    
    alpha = 0
    mu = .05
    rho = .3
    delta = .05
    epsilon = .9
    beta = .1
    h = 1
    I = makeI(nspec, nres)
    p = array(rep(0, nspec * nres^2), dim=c(nspec, nres, nres))
    
    init_state = c(rep(1, nspec), rep(1, nres), 0)
    sim = ode(y = init_state, times = seq(0, 30000, by = 1), func = model, parms = list(nspec = nspec, nres = nres, alpha = alpha, mu = mu, rho = rho, delta = delta, epsilon = epsilon, beta = beta, h = h, I = I, p = p))
    
    eql = tail(sim, 1)[-1]
    eql_abuns = eql[0:nspec]
    num_coexist = length(eql_abuns[eql_abuns > .1])
    df = c(eql_abuns, num_coexist) |> t()
    colnames(df) = c(paste0('N', seq(nspec)), 'num_coexist')
    df = as_tibble(df)
    df$I = list(I)
    
    overall_dists = matrix(0, nrow = nspec, ncol = nspec)
    
    for (i in 1:nspec) {
      for (j in 1:nspec) {
        overall_dists[i,j] = hamming_dist(I[i,], I[j,])
      }
    }
    
    coexist_traits = I[eql_abuns > .1,eql_abuns > .1]
    coex_dists = coex_dists = matrix(0, nrow = nrow(coexist_traits), ncol = nrow(coexist_traits))
    
    for (i in 1:nrow(coexist_traits)) {
      for (j in 1:nrow(coexist_traits)) {
        coex_dists[i,j] = hamming_dist(coexist_traits[i,], coexist_traits[j,])
      }
    }
    
    df$avg_ham = mean(overall_dists)
    df$avg_coex_ham = mean(coex_dists)
    df$modularity = DIRT_LPA_wb_plus(I)$modularity
    df$nestedness = nested(I)

    
    
    
    return(df)
  }

params = 
  expand_grid(
    numsim = 1:300,
    nspec = nspec,
    nres = nres
  )

result =
  params |>
  future_pmap_dfr(
    .f = simulation,
    .options = furrr_options(seed = NULL)
  )

hist(result$num_coexist)

ham_vec = rep(0, nspec)
mod_vec = rep(0, nspec)
nest_vec = rep(0, nspec)

for (i in 1:nspec) {
  ham_vec[i] = mean(result[result$num_coexist == i,]$avg_ham)
  mod_vec[i] = mean(result[result$num_coexist == i,]$modularity)
  nest_vec[i] = mean(result[result$num_coexist == i,]$nestedness)
}

plot(1:nspec, ham_vec, xlab = 'Number of Coexisting Species', ylab = 'Average Pairwise Hamming Distance')
plot(1:nspec, mod_vec, xlab = 'Number of Coexisting Species', ylab = 'Average Modularity')
plot(1:nspec, nest_vec, xlab = 'Number of Coexisting Species', ylab = 'Average Nestedness')
plot(result$avg_ham, result$avg_coex_ham, xlab = 'Average Initial Pairwise Hamming Distance', ylab = 'Average Coexisting Pairwise Hamming Distance')
abline(a = 0, b = 1)