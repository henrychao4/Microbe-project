n_vec = 1:4
nreps = 10000

dists = matrix(0, nrow = length(n_vec), ncol = nreps)

euclidean = \(vect1, vect2) {
  ans = sqrt(sum((vect1 - vect2)^2))
  return(ans)
}

for (i in 1:length(n_vec)) {
  for (j in 1:nreps) {
    pt1 = runif(n_vec[i], min = 0, max = 1)
    pt2 = runif(n_vec[i], min = 0, max = 1)
    dists[i,j] = euclidean(pt1, pt2)
  }
}

hist(dists[1,])
hist(dists[2,])
hist(dists[3,])
hist(dists[4,])