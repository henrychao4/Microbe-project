library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

nspec = 5
nres = 5

init_state = c(rep(1, nspec), rep(1, nres), 0)

U = 
  \(params, R){
    ans = matrix(0, nrow = params$nspec, ncol = params$nres)
    for (j in 1:params$nres) {
      for (i in 1:params$nspec) {
        ans[i,j] = params$beta * R[j] / (params$h + R[j]) * params$I[i,j]
      }
    }
    return(ans)
  }

omega = 
  \(params, R){
    ans = matrix(0, nrow = params$nspec, ncol = params$nres)
    for (i in 1:params$nspec) {
      ans[i,] = U(params, R)[i,] %*% params$p[i,,]
    }
    return(ans)
  }

gamma = 
  \(params, R){
    ans = rep(0, params$nspec)
    for (i in 1:params$nspec) {
      ans[i] = params$epsilon * sum(U(params, R)[i,])
    }
    return(ans)
  }

hamming_dist = \(x,y) {
  ans = sum(x != y)
  return(ans)
}

model = 
  \(t, state, params){
    N = state[1:params$nspec]
    R = state[(params$nspec + 1):(params$nspec + params$nres)]
    W = tail(state, 1)
    dNdt = params$alpha + gamma(params, R) * N - params$mu * N
    dRdt = params$rho - params$delta * R + colSums(omega(params, R)) - colSums(U(params, R) * N)
    dWdt = (1 - params$epsilon) * sum(U(params, R) * N)
    return(list(c(dNdt, dRdt, dWdt)))
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
    I = matrix(rbinom(nspec * nres, 1, 0.5), nrow = nspec, ncol = nres)
    p = array(rep(0, nspec * nres^2), dim=c(nspec, nres, nres))
    
    init_state = c(rep(1, nspec), rep(1, nres), 0)
    sim = ode(y = init_state, times = seq(0, 3000, by = 1), func = model, parms = list(nspec = nspec, nres = nres, alpha = alpha, mu = mu, rho = rho, delta = delta, epsilon = epsilon, beta = beta, h = h, I = I, p = p))
    
    eql = tail(sim, 1)[-1]
    eql_abuns = eql[0:nspec]
    num_coexist = length(eql_abuns[eql_abuns > .1])
    df = c(eql_abuns, num_coexist) |> t()
    colnames(df) = c(paste0('N', seq(nspec)), 'num_coexist')
    df = as_tibble(df)
    df$I = list(I)
    
    dists = matrix(0, nrow = nspec, ncol = nspec)
    
    for (i in 1:nspec) {
      for (j in 1:nspec) {
        dists[i,j] = hamming_dist(I[i,], I[j,])
      }
    }
    df$avg_ham = mean(dists)
    
    
    return(df)
  }

params = 
  expand_grid(
    numsim = 1:1000,
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

plot(result$avg_ham, result$num_coexist)

ham_vec = rep(0, nspec)
for (i in 1:nspec) {
  ham_vec[i] = mean(result[result$num_coexist == i,]$avg_ham)
}
plot(1:5, ham_vec, xlab = 'Number of Coexisting Species', ylab = 'Average Pairwise Hamming Distance')
