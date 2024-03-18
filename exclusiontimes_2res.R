library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)
library(plotly)

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

model = 
  \(t, state, params){
    R = state[1:2]
    N = state[-(1:2)]
    dRdt = with(params, R * (r * (1 - R / K) - t(C) %*% N))
    dNdt = with(params, N * (C %*% R - d))
    return(list(c(dRdt,dNdt)))
  }

simulation = 
  \(
    inv_c1,
    inv_c2
  ){
    
    r = c(1, 1)
    K = c(10, 10)
    res_C = matrix(c(.2, .1, .1, .2), nrow = 2, ncol = 2)
    res_d = c(.2, .2)
    
    inv_d = inv_c1 * (res_C[2,2] * res_d[1] - res_C[1,2] * res_d[2]) / (res_C[1,1] * res_C[2,2] - res_C[1,2] * res_C[2,1]) + inv_c2 * (res_C[1,1] * res_d[2] - res_C[2,1] * res_d[1]) / (res_C[1,1] * res_C[2,2] - res_C[1,2] * res_C[2,1])
    
    # pert_inv_c1 = runif(1, min = 0, max = .01)
    # pert_inv_c2 = runif(1, min = 0, max = .01)
    # pert_inv_d = runif(1, min = 0, max = .01)
    
    pert_inv_c1 = .000
    pert_inv_c2 = .000
    pert_inv_d = .000
    
    inv_c1 = inv_c1 - pert_inv_c1
    inv_c2 = inv_c2 - pert_inv_c2
    inv_d = inv_d + pert_inv_d
    
    inv_C = c(inv_c1, inv_c2)
    C = rbind(res_C, inv_C)
    d = c(res_d, inv_d)
    
    init_state = c(rep(10, length(r)), rep(10, length(res_d)), rep(10, length(inv_d)))
    
    sim = ode(y = init_state, times = seq(0, 3000, by = .1), func = model, parms = list(r = r, K = K, C = C, d = d)) |>
      as.data.frame()
    
    colnames(sim) = c('time', 'R1', 'R2', 'N1', 'N2', 'N3')
    
    t_extinct = min(filter(sim, N3 < .01)$time)
    
    A_plane = (C[2,2] * d[1] - C[1,2] * d[2]) / (C[1,1] * C[2,2] - C[1,2] * C[2,1])
    B_plane = (C[1,1] * d[2] - C[2,1] * d[1]) / (C[1,1] * C[2,2] - C[1,2] * C[2,1])
    C_plane = -1
    D_plane = 0
    dist_to_neutral = abs(A_plane * C[3,1] + B_plane * C[3,2] + C_plane * d[3] + D_plane) / sqrt(A_plane^2 + B_plane^2 + C_plane^2)
    
    
    output_vec = t(c(inv_c1, inv_c2, inv_d, dist_to_neutral, t_extinct)) |>
      as.data.frame()
    colnames(output_vec) = c('c_31', 'c_32', 'd_3', 'Distance_to_neutrality', 'Exclusion_Time')
    
    return(output_vec)
  }

params = expand_grid(
  inv_c1 = seq(.1, .3, by = .05),
  inv_c2 = seq(.1, .3, by = .05)
)

results =
  params |>
  future_pmap_dfr(
    .f = simulation,
    .options = furrr_options(seed = NULL)
  )

p1 = ggplot(data = results, aes(x= Distance_to_neutrality, y = Exclusion_Time, group=1)) +
  geom_point()

print(p1)

neutral_plane = function(plot_c31, plot_c32){
  C = matrix(c(.2, .1, .1, .2), nrow = 2, ncol = 2)
  d = c(.2, .2)
  d3 = plot_c31 * (C[2,2] * d[1] - C[1,2] * d[2]) / (C[1,1] * C[2,2] - C[1,2] * C[2,1]) + plot_c32 * (C[1,1] * d[2] - C[2,1] * d[1]) / (C[1,1] * C[2,2] - C[1,2] * C[2,1])
  return(d3)
}

plot_c31 = seq(.05, .25, by = .01)
plot_c32 = seq(.05, .25, by = .01)
plot_d = outer(plot_c31, plot_c32, neutral_plane)

p2 = plot_ly(data = results, x = ~c_31, y = ~c_32, z = ~d_3, type = "scatter3d", mode = "markers", color = ~Exclusion_Time) %>%
  add_trace(x = .2, y = .1, z = .2, inherit = FALSE) %>%
  add_trace(x = .1, y = .2, z = .2, color = 'resident', inherit = FALSE) %>%
  add_surface(x = plot_c31, y = plot_c32, z = plot_d, inherit = FALSE)

print(p2)