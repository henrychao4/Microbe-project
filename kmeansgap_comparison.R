#libraries
library(deSolve)
library(tidyverse)
library(cowplot)
library(SWKM)
library(furrr)
library(parallel)
library(tibble)


theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

plan(multisession(workers = detectCores() - 2))

source("Gap_Functions.r")

## Create pseudo data
set.seed(0)

num_clusters = 4 

# pseudo_data = 
#   sapply(seq(num_clusters), \(mu) rnorm(21e3/num_clusters, mu, mu/10)) |>
#   as.vector()


pseudo_data = 
  sapply(seq(num_clusters), \(mu) rnorm(21e3/num_clusters, mu, .1)) |>
  as.vector()

y = findInterval(pseudo_data, seq(-5, 10, .01))

df = 
  tibble(trait = y) |> 
  count(trait) |> 
  mutate(
    trait = scale(trait)[, 1],
    sp = seq(1:n())
  )

# Test
test = with(df, gap_function(trait, n, 2, 10, .05, 100))

#plots significance curve
plots = plot_gap(test)

plots$plot_clusters |> show()

plots$plot_significance |> show()