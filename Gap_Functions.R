## Functions

# WCSS: for a set of traits and abundances,  
#       calculates within-cluster sum of squares using 
#       weighted k-means clustering from package SKWM.

# gap_function: for a given data frame and a range of cluster counts,
#               calculates the gap statistic (and its variants from YanYe 2007),
#               finds the best number of clusters, determines significance.

# plot_gap: for a given object returned from `gap_function`, plots
#           the data colored by the estimated number of clusters 
#           according to each criterion, as well as the gap curve(s).

##new wcss function from the papers
WCSS = 
  \(
    trait, 
    abundance, 
    number_of_clusters, 
    seed
  ){
    if(seed > 0){
      set.seed(seed)
      abundance = sample(abundance)
    } 
    
    kW =
      kmeans.weight(
        x = matrix(c(trait, rep(0, length(trait))), ncol = 2),
        weight = 1e-6 + abundance,
        K = number_of_clusters
      )
    
    species = 
      tibble(
        trait = trait, 
        abundance = abundance,
        cluster = kW$cluster[,1]
      )
    
    centers = 
      species |> 
      group_by(cluster) |> 
      summarize(center = mean(trait))
    
    result = 
      species |>
      inner_join(centers, by = 'cluster') |>
      mutate(dic = abundance * (trait - center) ^ 2) |>
      group_by(cluster) |>
      summarize(
        Wk = sum(dic), 
        Wkbar = Wk / max(1, (sum(abundance) - 1)),
        .groups = 'drop'
      ) |>
      summarize(
        W = sum(Wk),
        Wbar = sum(Wkbar)
      )
    
    return(
      bind_cols(
        result,
        seed = seed,
        K = number_of_clusters
      )
    )
    
  }

gap_function = 
  \(
    trait, 
    abundance, 
    min_num_clusters, 
    max_num_clusters, 
    significance_threshold, 
    num_seeds, 
    ...
  ){
    parms = 
      expand_grid(
        seed = 0:(num_seeds - 1), 
        number_of_clusters = seq(min_num_clusters, max_num_clusters)
      )
    
    result = 
      parms |>
      future_pmap_dfr(
        .f = WCSS,
        trait = trait,
        abundance = abundance,
        .options = furrr_options(seed = NULL)
      ) |>
      mutate(
        lW = log(W),
        lWbar = log(Wbar)
      ) 
    
    expected = 
      result |>
      group_by(K) |>
      summarize(
        ElW = mean(lW),
        sdk = sd(lW),
        ElWbar = mean(lWbar),
        .groups = 'drop'
      )
    
    result_summary = 
      result |>
      full_join(expected, by = 'K') |>
      mutate(
        gapbar = ElWbar - lWbar,
        gap = ElW - lW,
        sk = sqrt(1 + 1 / num_seeds) * sdk
      ) |>
      group_by(seed) |>
      mutate(
        Dgap = c(NA, diff(gap)),        # Dgap[k] = gap[k] - gap[k-1]
        DDgap = -c(diff(Dgap), NA),     # DDgap[k] = Dgap[k] - Dgap[k+1]
        gapdiff = c(diff(gap), NA),     # gapdiff[k] = gap[k+1] - gap[k]
        test = gapdiff <= c(sk[-1], NA) # test[k] = is (gapdiff[k] <= sk[k+1])
        #         = is (gap[k+1] - gap[k] <= sk[k+1])
        #         = is (gap[k] >= gap[k+1] - sk[k+1])
      ) |>
      ungroup()
    
    tibshirani_results = 
      result_summary |>
      group_by(seed) |>
      summarize(
        khat = K[which(test == TRUE)[1]],
        statistic = gap[which(test == TRUE)[1]],
        .groups = 'drop'
      ) |>
      select(
        seed,
        khat,
        statistic
      ) |>
      mutate(criterion = 'Tibshirani2001')
    
    max_results =
      result_summary |>
      group_by(seed) |>
      slice_max(gap) |>
      mutate(khat = K) |>
      ungroup() |>
      select(
        seed,
        khat,
        statistic = gap
      ) |>
      mutate(criterion = 'MaxGap')
    
    yanye_gapbar_results = 
      result_summary |>
      group_by(seed) |>
      slice_max(gapbar) |>
      mutate(khat = K) |>
      ungroup() |>
      select(
        seed,
        khat,
        statistic = gapbar
      ) |>
      mutate(criterion = 'YanYe2007-Gapbar')
    
    yanye_DDgap_results = 
      result_summary |>
      group_by(seed) |>
      slice_max(DDgap) |>
      mutate(khat = K) |>
      ungroup() |>
      select(
        seed,
        khat,
        statistic = DDgap
      ) |>
      mutate(criterion = 'YanYe2007-DDGap')
    
    
    combined_results =
      tibshirani_results |>
      bind_rows(max_results) |>
      bind_rows(yanye_gapbar_results) |>
      bind_rows(yanye_DDgap_results)
    
    result_significance =
      combined_results |>
      group_by(criterion) |>
      summarize(
        significance_cutoff = quantile(statistic, 1 - significance_threshold), 
        .groups = 'drop'
      )
    
    final_report =
      combined_results |>
      filter(seed == 0) |>
      full_join(result_significance, by = 'criterion') |>
      mutate(significant = (statistic > significance_cutoff)) |>
      select(-seed)
    
    return(
      list(
        data = tibble(trait = trait, abundance = abundance),
        final_report = final_report, 
        result_summary = result_summary
      )
    )
  }

plot_gap =
  \(
    result, 
    criteria = 
      c(
        'Tibshirani2001', 
        'MaxGap', 
        'YanYe2007-Gapbar', 
        'YanYe2007-DDGap'
      )
  ){
    
    trait = result$data$trait
    abundance = result$data$abundance
    
    data = NULL
    for(.criterion in criteria){
      foo = 
        kmeans.weight(
          x = matrix(c(trait, rep(0, length(trait))), ncol = 2),
          weight = 1e-6 + abundance,
          K = result$final_report |> filter(criterion == .criterion) |> pull(khat)
        )
      
      data = 
        data |>
        bind_rows(
          tibble(
            method = .criterion,
            trait = trait,
            abundance = abundance,
            id = foo$cluster[,1],
            cluster = order(unique(id))[id]
          )
        )
    }
    
    plot_clusters =
      data |>
      mutate(cluster = factor(cluster)) |>
      ggplot() +
      geom_segment(
        aes(x = trait, xend = trait, y = 0, yend = abundance, color = cluster)
      )+
      facet_grid(~ method)
    
    plot_significance =
      result$result_summary |>
      filter(seed == 0) |>
      ggplot() +
      geom_line(aes(K, gap)) +
      geom_point(aes(K, gap)) +
      geom_line(aes(K, gapbar), color = 'red') +
      geom_point(aes(K, gapbar), color = 'red') +
      geom_line(aes(K, DDgap), color = 'blue') +
      geom_point(aes(K, DDgap), color = 'blue') +
      geom_hline(
        aes(
          yintercept = significance_cutoff, color = criterion
        ), 
        data = result$final_report
      )
    
    
    
    return(list(plot_clusters = plot_clusters, plot_significance = plot_significance))
    
  }