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

tibshirani = \(kmeans) {
  withinss = kmeans$withinss
  size = kmeans$size
  W = 0
  k = length(withinss)
  for (i in 1:k) {
    W = W + withinss[i] / (2 * size[i])
  }
  return(W)
}

yanye = \(kmeans) {
  withinss = kmeans$withinss
  size = kmeans$size
  W = 0
  k = length(withinss)
  for (i in 1:k) {
    if (size[i] > 1) {
      W = W + withinss[i] / (2 * size[i] * (size[i] - 1))
    }
  }
  return(W)
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

best_k_means = \(data, k, nruns) {
  df = list()
  tot.withinss_vec = rep(0, nruns)
  for (i in 1:nruns) {
    k_means = kmeans(data, k)
    df[[i]] = k_means
    tot.withinss_vec[i] = k_means$tot.withinss
  }
  idx = which.min(tot.withinss_vec)
  return(df[[idx]])
}


nruns = 10
k_vec = 1:10
tibs_logW_vec = rep(0, length(k_vec))
yanye_logW_vec = rep(0, length(k_vec))

for (i in 1:length(k_vec)) {
  k_means = best_k_means(data, k_vec[i], nruns)
  tibs_logW_vec[i] = log(tibshirani(k_means))
  yanye_logW_vec[i] = log(yanye(k_means))
}

n_bootstrap = 30
boot_tibs_logW = matrix(0, nrow = n_bootstrap, ncol = length(k_vec))
boot_yanye_logW = matrix(0, nrow = n_bootstrap, ncol = length(k_vec))

for (j in 1:length(k_vec)) {
  for (i in 1:n_bootstrap){
    null_data = get_null_data(data)
    k_means_null = best_k_means(null_data, k_vec[j], nruns)
    boot_tibs_logW[i,j] = log(tibshirani(k_means_null))
    boot_yanye_logW[i,j] = log(yanye(k_means_null))
  }
}

avg_boot_tibs_logW_vec = colSums(boot_tibs_logW) / n_bootstrap
avg_boot_yanye_logW_vec = colSums(boot_yanye_logW) / n_bootstrap

plot(k_vec, avg_boot_tibs_logW_vec, type = 'b', col = 'red', main = 'Tibshirani W', ylim = c(-5, 10))
lines(k_vec, tibs_logW_vec, type = 'b', col = 'blue', main = 'Yan & Ye W')

tibs_gap = avg_boot_tibs_logW_vec - tibs_logW_vec
plot(k_vec, tibs_gap, type = 'b')

yanye_gap = avg_boot_yanye_logW_vec - yanye_logW_vec
plot(k_vec, yanye_gap, type = 'b')