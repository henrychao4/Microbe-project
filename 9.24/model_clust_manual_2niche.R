library(purrr)

set.seed(1)

weighted_hamming_dist = function(vec1, vec2, weight) {
  distance = sum(vec1 != vec2) * weight
  return(distance)
}

weighted_mode = \(x, weights) {
  unique_vals = unique(x)
  weighted_freq = sapply(unique_vals, function(val) {
    sum(weights[x == val])
  })
  unique_vals[which.max(weighted_freq)]
  #unique_vals[which.max(rank(weighted_freq, ties.method = "random"))]
}

init_centers = \(data, weights, k) {
  n = nrow(data)
  centers = matrix(0, nrow = k, ncol = ncol(data))
  first_center_idx = sample(1:n, 1)
  centers[1, ] = data[first_center_idx, ]
  
  for (i in 2:k) {
    min_dist = apply(data, 1, function(pt) {
      min(sapply(1:(i-1), function(c) {
        weighted_hamming_dist(pt, centers[c, ], weights)
      }))
    })
    
    prob = min_dist^2
    
    next_center_idx = sample(1:n, 1, prob = prob)
    centers[i, ] = data[next_center_idx, ]
  }
  
  return(centers)
}

weighted_k_modes = \(data, weights, k, max_iter) {
  centers = init_centers(data, weights, k)
  #centers = data[sample(1:nrow(data), k), ]
  clusters = rep(0, nrow(data))
  convergence = 0
  iter = 1
  
  while (convergence == 0) {
    for (i in 1:nrow(data)) {
      centers_dists = rep(0, k)
      for (j in 1:k) {
        centers_dists[j] = weighted_hamming_dist(data[i,], centers[j,], weights[i])
      }
      clusters[i] = which.min(centers_dists)
    }
    
    new_centers = centers
    for (i in 1:k) {
      cluster_pts = data[clusters == i, , drop = FALSE]
      clusters_wts = weights[clusters == i]
      new_centers[i,] = apply(cluster_pts, 2, weighted_mode, weights = weights)
    }
    
    if (all(centers == new_centers) || iter == max_iter) {
      convergence = 1
    }
    
    centers = new_centers
    iter = iter + 1
  }
  list(clusters = clusters, centers = centers)
}

calculate_within_cluster_error = \(data, clusters, centers, weights) {
  total_error = 0
  
  for (i in unique(clusters)) {
    cluster_data = data[clusters == i, , drop = FALSE]
    cluster_weights = weights[clusters == i]
    cluster_center = centers[i, ]
    
    for (j in 1:nrow(cluster_data)) {
      dissimilarity = sum((cluster_data[j, ] != cluster_center) * cluster_weights[j])
      total_error = total_error + dissimilarity
    }
  }
  
  return(total_error)
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

best_wKModes = \(data, weights, modes, max_iter, nruns) {
  df = list()
  tot.withindiff_vec = rep(0, nruns)
  for (i in 1:nruns) {
    wk_modes = weighted_k_modes(data = data, weights = weights, k = modes, max_iter = 100)
    err = calculate_within_cluster_error(data = data, clusters = wk_modes$clusters, centers = wk_modes$centers, weights)
    df[[i]] = wk_modes
    tot.withindiff_vec[i] = err
  }
  idx = which.min(tot.withindiff_vec)
  return(df[[idx]])
}

find_optimal_clusters = \(data, weights, max_iter, n_bootstrap, k_max, nruns) {
  k_vec = 2:k_max
  clust_dispersion = rep(0, length(k_vec))
  for (i in 1:length(k_vec)) {
    errs = rep(0, 10)
    for (j in 1:10) {
      kmode = best_wKModes(data, weights, k_vec[i], max_iter = max_iter, nruns = nruns)
      errs[j] = calculate_within_cluster_error(data, kmode$clusters, kmode$centers, weights)
    }
    clust_dispersion[i] = log(min(errs))
  }
  
  n_bootstrap = n_bootstrap
  null_clust_dispersion = matrix(0, nrow = n_bootstrap, ncol = length(k_vec))
  
  for (i in 1:n_bootstrap) {
    null_data = get_null_data(data)
    for (j in 1:length(k_vec)) {
      errs = rep(0, 10)
      for (k in 1:10) {
        kmode = best_wKModes(null_data, weights, k_vec[j], max_iter = 100, nruns = nruns)
        errs[k] = calculate_within_cluster_error(null_data, kmode$clusters, kmode$centers, weights)
      }
      null_clust_dispersion[i,j] = log(min(errs))
    }
  }
  
  null_clust_dispersion_avg = colSums(null_clust_dispersion) / n_bootstrap
  
  k1_center = apply(data, 2, weighted_mode, weights = weights)
  k1_dispersion = 0
  
  for (j in 1:nrow(data)) {
    dissimilarity = sum((data[j, ] != k1_center) * weights[j])
    k1_dispersion = k1_dispersion + dissimilarity
  }
  k1_dispersion = log(k1_dispersion)
  
  null_k1_dispersion = rep(0, 10)
  for (i in 1:10) {
    null_data = get_null_data(data)
    null_k1_center = apply(null_data, 2, weighted_mode, weights = weights)
    disp = 0
    for (j in 1:nrow(data)) {
      dissimilarity = sum((null_data[j, ] != null_k1_center) * weights[j])
      disp = disp + dissimilarity
      null_k1_dispersion[i] = log(disp)
    }
  }
  
  null_k1_dispersion_avg = mean(null_k1_dispersion)
  
  k_vec = c(1, k_vec)
  clust_dispersion = c(k1_dispersion, clust_dispersion)
  null_clust_dispersion_avg = c(null_k1_dispersion_avg, null_clust_dispersion_avg)
  
  gap = null_clust_dispersion_avg - clust_dispersion
  
  plot(k_vec, gap, type = 'b')
  
  return(which.max(gap))
}

set.seed(8)

nspec = 10
nres = 10

makeI = \(nspec, nres) {
  I = matrix(0, nrow = nspec, ncol = nres)
  cons_vec = c(rep(0, nres/2), rep(1, nres/2))
  for (i in 1:nspec) {
    I[i,] = sample(cons_vec, nres, replace = F)
  }
  return(I)
}

I = makeI(nspec, nres)

params = list(
  nspec = nspec,
  nres = nres,
  alpha = rep(0, nspec),
  mu = rep(.05, nspec),
  rho = c(rep(0, nres-2), rep(.6, 2)),
  delta = rep(.05, nres),
  epsilon = .9,
  beta = 1,
  h = 1,
  I = I,
  p = array(rep(0, nspec * nres^2), dim = c(nspec, nres, nres))
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
p = ggplot(abuns.df, aes(time, value, color = variable)) + geom_line() + theme_classic()
#print(p)

eql = tail(sim, 1)[-1]
eql_abuns = eql[0:nspec]
num_coexist = length(eql_abuns[eql_abuns > .1])
df = c(eql_abuns, num_coexist) |> t() |> as_tibble()
colnames(df) = c(paste0('N', seq(nspec)), 'num_coexist')

init_abun_gap = find_optimal_clusters(data = I, weights = rep(1, nspec), max_iter = 100, n_bootstrap = 10, k_max = 5, nruns = 10)

eql_abun_gap = find_optimal_clusters(data = I, weights = eql_abuns + .01, max_iter = 100, n_bootstrap = 10, k_max = 5, nruns = 10)

#this is actually demonstrating clustering