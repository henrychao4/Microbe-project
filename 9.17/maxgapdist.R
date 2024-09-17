library(cluster)

set.seed(1)

nclust = 3
similarity = 5
n = 10
m = 5
sd = .1

get_data = \(nclust, similarity, n, m, sd) {
  data = matrix(0, nrow = nclust * n, ncol = m)
  means = seq(0, similarity * nclust, similarity)
  for (i in 0:(nclust-1)) {
    for (j in 1:m) {
      data[((n*i)+1):(n*(i+1)), j] = rnorm(n, mean = means[i+1], sd = sd)
    }
  }
  
  return(data)
}

data = get_data(nclust, similarity, n, m, sd)

get_null_data = \(data) {
  n = nrow(data)
  m = ncol(data)
  null_data = matrix(0, nrow = n, ncol = m)
  
  for (j in 1:m) {
    col_j = data[,j]
    null_data[,j] = runif(n, min = min(col_j), max = max(col_j))
  }
  
  return(null_data)
}

n_bootstrap = 1000
boot_k = rep(0, n_bootstrap)
boot_max_gaps = rep(0, n_bootstrap)
k_max = 5

for (i in 1:n_bootstrap) {
  null_data = get_null_data(data)
  clusgap = clusGap(null_data, FUNcluster = kmeans, K.max = k_max, B = 10, verbose = F)
  boot_k[i] = which.max(clusgap$Tab[,3])
  boot_max_gaps[i] = max(clusgap$Tab[,3])
}

clusgap = clusGap(data, FUNcluster = kmeans, K.max = k_max, B = 100, verbose = F)
true_max_gap = max(clusgap$Tab[,3])

hist(boot_k)
hist(boot_max_gaps)

hist(boot_max_gaps[boot_k == 1])
hist(boot_max_gaps[boot_k == 2])
hist(boot_max_gaps[boot_k == 3])
hist(boot_max_gaps[boot_k == 4])
