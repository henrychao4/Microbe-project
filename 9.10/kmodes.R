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
library(klaR)
library(MEDseq)

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

set.seed(1)

nclust = 5 #number of clusters
n = 6 #number of species per cluster
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
      flip = rbernoulli(1, .05)
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

I = makeI(nclust, nspec, nres)

I = as.data.frame(I)

rep

clust = wKModes(I, modes = 5)
print(clust)

k_vec = 1:8
W_vec = rep(0, length(k_vec))
for (k in 1:length(k_vec)) {
  k = k_vec[k]
  W = rep(0, 10)
  for (i in 1:length(W)) {
    W[i] = wKModes(I, modes = k)$tot.withindiff
  }
  W_vec[k] = min(W)
}

#elbow method shows that k = 5 is optimal
plot(k_vec, W_vec, type = 'b')

#trying to implement gap statistic
set.seed(4)
k_vec = 1:7
logW_mat = matrix(0, nrow = 10, ncol = length(k_vec))
for (i in 1:nrow(logW_mat)) {
  gaps = clusGap(I, FUN = wKModes, K.max = length(k_vec), B = 2)
  tbl = gaps$Tab
  logW_mat[i,] = tbl[,1]
}
W_vec = apply(logW_mat, 2, min)

gap_stat = clusGap(I, FUN = wKModes, K.max = 7, B = 50)



tbl = gap_stat$Tab

gap = tbl[,2] - W_vec

plot(W_vec, type = 'b', col = 'red')
lines(tbl[,2], type = 'b', col = 'blue')

plot(gap, type = 'b', col = 'red')

