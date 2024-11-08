library(pracma)

nspec = 500
nres = 500

ntraits = 4

circ_dist = \(vec1, vec2) {
  y = abs(outer(vec1, vec2, '-'))
  idx = which(y > 0.5)
  y[idx] = 1 - y[idx]
  return(y)
}

res_traits = matrix(runif(nres * ntraits, min = 0, max = 1), nrow = nres, ncol = ntraits)
res_dists = matrix(0, nrow = nres, ncol = nres)

dists_list = list()
for (k in 1:ntraits){
  dists_list[[k]] = circ_dist(res_traits[,k], res_traits[,k])
}
euclid_circ_dist = matrix(0, nrow = nres, ncol = nres)
inf_norm_circ_dist = matrix(0, nrow = nres, ncol = nres)
neg_inf_norm_circ_dist = matrix(0, nrow = nres, ncol = nres)

for (i in 1:nres) {
  for (j in 1:nres) {
    elements_vec = unlist(lapply(dists_list, function(mat) mat[i, j]))
    euclid_circ_dist[i,j] = Norm(elements_vec, p = 2)
    inf_norm_circ_dist[i,j] = Norm(elements_vec, p = Inf)
    neg_inf_norm_circ_dist[i,j] = Norm(elements_vec, p = -Inf)
    
  }
}

mod_euclid_circ_dist = euclid_circ_dist
diag(mod_euclid_circ_dist) = Inf

mod_inf_norm_circ_dist = inf_norm_circ_dist
diag(mod_inf_norm_circ_dist) = Inf

mod_neg_inf_norm_circ_dist = neg_inf_norm_circ_dist
diag(mod_neg_inf_norm_circ_dist) = Inf

nnd_2 = apply(mod_euclid_circ_dist, 1, min)
nnd_inf = apply(mod_inf_norm_circ_dist, 1, min)
nnd_neg_inf = apply(mod_neg_inf_norm_circ_dist, 1, min)


hist(nnd_2, xlab = "Nearest Neighbor Distance 2 Norm", main = paste0('Dimension = ', as.character(ntraits)))
hist(nnd_inf, xlab = "Nearest Neighbor Distance Inf Norm", main = paste0('Dimension = ', as.character(ntraits)))
hist(nnd_neg_inf, xlab = "Nearest Neighbor Distance Negative Inf Norm", main = paste0('Dimension = ', as.character(ntraits)))

connectance_2 = 