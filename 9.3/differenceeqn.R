set.seed(1)

tmax = 50
nspec = 4
nres = 4

N = matrix(0, nrow = tmax, ncol = nspec)
R = matrix(0, nrow = tmax, ncol = nres)


makeI = \(nspec, nres) {
  I = matrix(0, nrow = nspec, ncol = nres)
  cons_vec = c(rep(1, nres/2), rep(0, nres/2))
  for (i in 1:nspec) {
    I[i,] = sample(cons_vec, nres, replace = F)
  }
  while (sum(duplicated(I)) > 0) {
    for (i in 1:nspec) {
      I[i,] = sample(cons_vec, nres, replace = F)
    }
  }
  return(I)
}

I = makeI(nspec, nres)

N[1,] = rep(10, nspec)
R[1,] = rep(10, nspec)

spec_immigration = 
  \(N){
    index = sample(1:length(N), 1)
    N[index] = N[index] + 1
    return(N)
  }

spec_death = 
  \(N){
    index = sample(1:length(N), 1, prob = N)
    N[index] = N[index] - 1
    return(N)
  }

res_supply = 
  \(R){
    index = sample(1:length(R), 1)
    R[index] = R[index] + 1
    return(R)
  }

res_decay = 
  \(R){
    index = sample(1:length(R), 1, prob = R)
    R[index] = R[index] - 1
    return(R)
  }

consumption = 
  \(N, R, I){
    consumer_index = sample(1:length(N), 1, prob = N)
    row_I = I[consumer_index,]
    selected_res_abuns = R[row_I == 1]
    res_index = sample(1:length(selected_res_abuns), 1, prob = selected_res_abuns)
    R[res_index] = R[res_index] - 1
    N[consumer_index] = N[consumer_index] + 1
    return(list(N,R))
  }

for (i in 1:(tmax - 1)) {
  event = sample(1:5, 1)
  if (event == 1) {
    N[(i+1), ] = spec_immigration(N[i, ])
    R[(i+1), ] = R[i, ]
  }
  
  if (event == 2) {
    N[(i+1), ] = spec_death(N[i, ])
    R[(i+1), ] = R[i, ]
  }
  
  if (event == 3) {
    N[(i+1), ] = N[i, ]
    R[(i+1), ] = res_supply(R[i, ])
  }
  
  if (event == 4) {
    N[(i+1), ] = N[i, ]
    R[(i+1), ] = res_decay(R[i, ])
  }
  
  if (event == 5) {
    N[(i+1), ] = consumption(N[i, ], R[i, ], I)[1]
    R[(i+1), ] = consumption(N[i, ], R[i, ], I)[2]
  }
  
}

