set.seed(2)

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
weights = runif(nspec)
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
}

set.seed(1)

weighted_k_modes <- function(data, weights, k, max_iter) {
  centers = data[sample(1:nrow(data), k), ]
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

kmode = weighted_k_modes(data, weights, k = 5, max_iter = 100)
calculate_within_cluster_error(data, kmode$clusters, kmode$centers, weights)

k_vec = 2:8
error = rep(0, length(k_vec))
for (i in 1:length(k_vec)) {
  errs = rep(0, 10)
  for (j in 1:10) {
    kmode = weighted_k_modes(data, weights, k_vec[i], max_iter = 100)
    errs[j] = calculate_within_cluster_error(data, kmode$clusters, kmode$centers, weights)
  }
  error[i] = min(errs)
}

plot(k_vec, error, type = 'b')
