library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(matlib)

plan(multisession(workers = detectCores() - 2))
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

model = 
  \(t, state, params){
    R = state[1]
    N = state[-1]
    dRdt = params$inflow - params$outflow * R - sum(params$C * R * N)
    dNdt = N * (params$C * R - params$death)
    return(list(c(dRdt,dNdt)))
  }

simulation = 
  \(
    nspec,
    inflow,
    outflow,
    c1,
    c2,
    d1,
    d2
  ){
    C = c(c1, c2)
    death = c(d1, d2)
    
    init_state = c(10, rep(10, nspec))
    sim = ode(y = init_state, times = seq(0, 10000, by = 1), func = model, parms = list(inflow = inflow, outflow = outflow, C = C, death = death)) |>
      as.data.frame() 
    colnames(sim) = c('time', 'R', 'N1', 'N2')
    
    t_extinct = min(filter(sim, N2 < 1)$time)
    output_vec = c(t_extinct)
    names_vec = c(c2)
    output_df = data.frame(names_vec, output_vec)
    colnames(output_df) = c('C_2', 'Exclusion_Time')
    
    return(output_df)
  }

params = expand_grid(
  nspec = 2,
  inflow = 1,
  outflow = .1,
  c1 = 1,
  c2 = c(seq(0, .46, by = .01), seq(.47, .499, by = .001)),
  d1 = .2,
  d2 = .1
)

results =
  params |>
  future_pmap_dfr(
    .f = simulation,
    .options = furrr_options(seed = NULL)
  )

plot(results)

p = ggplot(data=results, aes(x= C_2, y = Exclusion_Time, group=1)) +
  geom_line()+
  geom_point()

print(p)