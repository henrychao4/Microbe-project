library(data.table)
library(cluster)
library(factoextra)
library(purrr)

set.seed(1)

nclust = 5
n = 10
nres = 20
niche_width = .3

clusts = seq(0, 1, length.out = nclust)
resources = seq(0, 1, length.out = nres)

dists = outer(clusts, resources, FUN = \(x, y) abs(x - y))
clust_traits = exp(-(dists / niche_width) ^ 2)

N = nclust * n
I = matrix(0, nrow = N, ncol = nres)

makeI = \(n, N, nres) {
  for (i in 1:N) {
    I[i,] = jitter(clust_traits[ceiling(i/n),], amount = .1)
  }
  return(I)
}

cont_I = makeI(n, N, nres)
disc_I = (cont_I > .5) * 1

calculate_silhouette = function(d, hc, k) {
  clusters = cutree(hc, k = k)
  sil_width = silhouette(clusters, dist = d)
  avg_silhouette = mean(sil_width[, "sil_width"])
  return(avg_silhouette)
}

max_clusters = 10

nsample = 500
known_traits = seq(1, nres)
cont_optimal_clusters = rep(0, length(known_traits))
disc_optimal_clusters = rep(0, length(known_traits))
for (i in 1:nres) {
  sample_cont_optimal_clusters = rep(0, nsample)
  sample_disc_optimal_clusters = rep(0, nsample)
  for (j in 1:nsample) {
    sample_cont_I = cont_I[,sample(ncol(I), size = i)]
    sample_disc_I = disc_I[,sample(ncol(I), size = i)]
    sample_cont_d = dist(sample_cont_I, method = "euclidean")
    sample_disc_d = dist(sample_disc_I, method = "binary")
    sample_cont_hc = hclust(sample_cont_d)
    sample_disc_hc = hclust(sample_disc_d)
    sample_cont_silhouette_scores = rep(-Inf, max_clusters)
    sample_disc_silhouette_scores = rep(-Inf, max_clusters)
    for (k in 2:max_clusters) {
      sample_cont_silhouette_scores[k] = calculate_silhouette(as.dist(sample_cont_d), sample_cont_hc, k)
      sample_disc_silhouette_scores[k] = calculate_silhouette(as.dist(sample_disc_d), sample_disc_hc, k)
    }
    sample_cont_optimal_clusters[j] = which.max(sample_cont_silhouette_scores)
    sample_disc_optimal_clusters[j] = which.max(sample_disc_silhouette_scores)
  }
  cont_optimal_clusters[i] = mean(sample_cont_optimal_clusters)
  disc_optimal_clusters[i] = mean(sample_disc_optimal_clusters)
}

plot(known_traits, cont_optimal_clusters, type = "b", main = 'Continuous traits', xlab = "Number of Known Traits", ylab = "Average Number of Optimal Clusters Found")
plot(known_traits, disc_optimal_clusters, type = "b", main = 'Discrete traits', xlab = "Number of Known Traits", ylab = "Average Number of Optimal Clusters Found")