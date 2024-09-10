library(purrr)

set.seed(1)

nclust = 5 #number of clusters
n = 6 #number of species per cluster
nspec = nclust * n
nres = 20 #number of traits per species

makeI = \(nclust, nspec, nres) {
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
      flip = rbernoulli(1, 0.1)
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

data = makeI(nclust, nspec, nres)
#weights = runif(nspec)
weights = c(rep(.01, nspec/5), rep(.9, 4 * nspec / 5))
#weights = rep(1, nspec)

weighted_hamming_dist = function(vec1, vec2, weight) {
  distance <- sum(vec1 != vec2) * weight
  return(distance)
}

weighted_mode = function(x, weights) {
  unique_vals = unique(x)
  weighted_freq = sapply(unique_vals, function(val) {
    sum(weights[x == val])
  })
  unique_vals[which.max(weighted_freq)]
  #unique_vals[which.max(rank(weighted_freq, ties.method = "random"))]
}

init_centers = function(data, weights, k) {
  n = nrow(data)
  centers = matrix(nrow = k, ncol = ncol(data))
  first_center_idx = sample(1:n, 1)
  centers[1, ] = data[first_center_idx, ]
  
  for (i in 2:k) {
    min_diss = apply(data, 1, function(obs) {
      min(sapply(1:(i-1), function(c) {
        weighted_hamming_dist(obs, centers[c, ], weights)
      }))
    })
    
    prob <- min_diss^2 / sum(min_diss^2)
    
    next_center_idx = sample(1:n, 1, prob = prob)
    centers[i, ] = data[next_center_idx, ]
  }
  
  return(centers)
}

weighted_k_modes = function(data, weights, k, max_iter) {
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

calculate_within_cluster_error = function(data, clusters, centers, weights) {
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

get_null_data = function(data) {
  n = nrow(data)
  m = ncol(data)
  null_data = matrix(0, nrow = n, ncol = m)
  
  for (j in 1:m) {
    null_data[,j] = sample(data[,j], n, replace = T)
  }
  
  return(null_data)
}

k_vec = 2:8
clust_dispersion = rep(0, length(k_vec))
for (i in 1:length(k_vec)) {
  errs = rep(0, 10)
  for (j in 1:10) {
    kmode = weighted_k_modes(data, weights, k_vec[i], max_iter = 100)
    errs[j] = calculate_within_cluster_error(data, kmode$clusters, kmode$centers, weights)
  }
  clust_dispersion[i] = log(min(errs))
}

n_bootstrap = 10
null_clust_dispersion = matrix(0, nrow = n_bootstrap, ncol = length(k_vec))

for (i in 1:n_bootstrap) {
  null_data = get_null_data(data)
  for (j in 1:length(k_vec)) {
    errs = rep(0, 10)
    for (k in 1:10) {
      kmode = weighted_k_modes(null_data, weights, k_vec[j], max_iter = 100)
      errs[k] = calculate_within_cluster_error(null_data, kmode$clusters, kmode$centers, weights)
    }
    null_clust_dispersion[i,j] = log(min(errs))
  }
}

null_clust_dispersion_avg = colSums(null_clust_dispersion) / n_bootstrap

gap = null_clust_dispersion_avg - clust_dispersion

plot(k_vec, clust_dispersion, type = 'b', col = 'blue')
lines(k_vec, null_clust_dispersion_avg, type = 'b', col = 'red')
legend(5,5.5,c("Actual log W","Bootstrapped log W"), lwd=c(5,2), col=c("blue","red"), y.intersp=1.5)

plot(k_vec, gap, type = 'b')

#now have to implement to the actual model