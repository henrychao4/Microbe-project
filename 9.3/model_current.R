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

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

set.seed(1)

nspec = 20
nres = 20

makeI = \(nspec, nres) {
  I = matrix(0, nrow = nspec, ncol = nres)
  cons_vec = c(rep(0, nres/2), rep(1, nres/2))
  for (i in 1:nspec) {
    I[i,] = sample(cons_vec, nres, replace = F)
  }
  return(I)
}

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(0.01, nspec),
  mu = rep(.05, nspec),
  rho = rep(.3, nres),
  delta = rep(.05, nres),
  # mu = rep(0, nspec),
  # rho = rep(0, nres),
  # delta = rep(0, nres),
  epsilon = .9,
  beta = 1,
  h = 1,
  I = makeI(nspec, nres),
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

find_num_alpha =
  \(eql, eql_abuns, model, params){
    dt = .001
    dNj = .001
    fwd_result = ode(y = eql, times = seq(0, dt, by = dt/100), func = model, parms = params)
    fwd_state = as.numeric(tail(fwd_result, n = 1)[-1])
    fwd_abuns = fwd_state[0:nspec]
    num_alpha = matrix(0, params$nspec, params$nspec)
    for (j in 1:nspec) {
      pert_state = eql
      pert_state[j] = pert_state[j] + dNj
      pert_abuns = pert_state[0:nspec]
      fwd_pert_result = ode(y = pert_state, times = seq(0,dt, by = dt/100), func = model, parms = params)
      fwd_pert_state = as.numeric(tail(fwd_pert_result, n = 1)[-1])
      fwd_pert_abuns = fwd_pert_state[0:nspec]
      num_alpha[,j] = (fwd_pert_abuns - pert_abuns) / (dt * dNj * pert_abuns) - (fwd_abuns - eql_abuns) / (dt * dNj * eql_abuns)
    }
    return(num_alpha)
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

num_alpha = find_num_alpha(eql, eql_abuns, model, params)

weights = matrix(0, nrow = nspec, ncol = nspec)
for (i in 1:nspec) {
  for (j in 1:nspec) {
    weights[i,j] = eql_abuns[i] * eql_abuns[j]
  }
}

coex_weight_d = as.matrix(dist(params$I, method = "binary")) * weights
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

optimal_clusters = print(which.max(silhouette_scores))

membership = cutree(coex_weight_hc, k = optimal_clusters)
plot(coex_weight_hc)


nboot = 1000
boot_optimal_clusters = rep(0, nboot)
boot_silhouette_peaks = rep(0, nboot)
for (i in 1:nboot) {
  boot_I = matrix(0, nrow = nspec, ncol = nres)
  for (j in 1:nres) {
    col_j = params$I[,j]
    boot_I[,j] = sample(col_j, replace = F)
  }
  boot_d = as.dist(as.matrix(dist(boot_I, method = "binary")) * weights)
  boot_hc = hclust(boot_d)
  boot_silhouette_scores = rep(-Inf, max_clusters)
  for (k in 2:max_clusters) {
    boot_silhouette_scores[k] = calculate_silhouette(as.dist(boot_d), boot_hc, k)
  }
  boot_optimal_clusters[i] = which.max(boot_silhouette_scores)
  boot_silhouette_peaks[i] = max(boot_silhouette_scores)
}

hist(boot_optimal_clusters, xlab = "Optimal Amount of Clusters", ylab = "Frequency")
hist(boot_silhouette_peaks, xlab = "Peak Silhouette Score", ylab = "Frequency")

lb_nclust_ci = quantile(boot_optimal_clusters, probs = .025)
ub_nclust_ci = quantile(boot_optimal_clusters, probs = .975)
print(paste0('The 95% confidence interval for the number of clusters found by silhouette score is (', as.character(lb_nclust_ci), ', ', as.character(ub_nclust_ci), ')'))
print(paste0('The number of clusters found in the true data is ', as.character(optimal_clusters)))

lb_peak_ci = quantile(boot_silhouette_peaks, probs = .025)
ub_peak_ci = quantile(boot_silhouette_peaks, probs = .975)
print(paste0('The 95% confidence interval for the peak silhouette score is (', as.character(lb_peak_ci), ', ', as.character(ub_peak_ci), ')'))
print(paste0('The peak silhouette score found in the true data is ', as.character(max(silhouette_scores))))