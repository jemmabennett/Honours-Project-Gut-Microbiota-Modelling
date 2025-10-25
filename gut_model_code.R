# Gut Microbiome Modelling R Script 
# Jemma Bennett (25136283)

# libraries needed

library(deSolve)
library(tidyverse)
library(igraph)
library(ggraph)
library(patchwork)
library(bipartite)
library(rootSolve)
library(numDeriv)

# BASIC MODEL SETUP

# The main ODE function
gut_model <- function(t, y, parms) {
  M <- y[1:parms$n_m]
  B <- y[(parms$n_m + 1):(parms$n_m + parms$n_b)]
  
  M[M < 0] <- 0
  B[B < 0] <- 0
  
  # Build the interaction matrix
  c_matrix <- rbind(parms$c_prod, parms$c_cons)
  c_nonmono_current <- matrix(0, nrow = parms$n_m_nonmono, ncol = parms$n_b)
  M_nonmono <- M[(parms$n_m_prod + parms$n_m_cons + 1):parms$n_m]
  
  # Calculate non-monotonic interactions
  for (i in 1:parms$n_m_nonmono) {
    for (j in 1:parms$n_b) {
      ratio <- ifelse(B[j] > 1e-9, M_nonmono[i] / B[j], 0)
      hill_term <- (parms$w[i, j]^parms$K) / (parms$w[i, j]^parms$K + ratio^parms$K)
      c_nonmono_current[i, j] <- parms$c_nonmono_L[i, j] + 
        (parms$c_nonmono_U[i, j] - parms$c_nonmono_L[i, j]) * hill_term
    }
  }
  
  c_matrix <- rbind(c_matrix, c_nonmono_current)
  
  # Calculate derivatives
  dM_dt <- parms$s_i - M * (c_matrix %*% B)
  growth_term <- parms$c0 * (t(c_matrix) %*% M)
  logistic_factor <- max(0, (parms$k - sum(B)))
  dB_dt <- B * (-parms$d + growth_term) * logistic_factor
  
  return(list(c(dM_dt, dB_dt)))
}

# SINGLE SIMULATION

run_single_simulation <- function() {
  set.seed(42)
  
  # System size
  n_b <- 3
  n_m_prod <- 3
  n_m_cons <- 7
  n_m_nonmono <- 2
  n_m <- n_m_prod + n_m_cons + n_m_nonmono
  
  # Parameters
  c0 <- 0.1
  k <- 14.0
  d <- rep(0.1, n_b)
  K <- 2
  
  # Initial conditions
  B_init <- runif(n_b, 0.05, 0.2)
  M_prod_init <- runif(n_m_prod, exp(-3), exp(-2))
  M_cons_init <- runif(n_m_cons, 0.5, 1.0)
  M_nonmono_init <- runif(n_m_nonmono, exp(-3), exp(-2))
  y0 <- c(M_prod_init, M_cons_init, M_nonmono_init, B_init)
  names(y0) <- c(paste0("M", 1:n_m), paste0("B", 1:n_b))
  
  # Interaction matrices
  c_prod <- matrix(runif(n_m_prod * n_b, -1, 0), nrow = n_m_prod, ncol = n_b)
  c_cons <- matrix(runif(n_m_cons * n_b, 0, 1), nrow = n_m_cons, ncol = n_b)
  c_nonmono_L <- matrix(runif(n_m_nonmono * n_b, -1, 0), nrow = n_m_nonmono, ncol = n_b)
  c_nonmono_U <- matrix(runif(n_m_nonmono * n_b, 0, 1), nrow = n_m_nonmono, ncol = n_b)
  w <- matrix(runif(n_m_nonmono * n_b, 0.1, 1), nrow = n_m_nonmono, ncol = n_b)
  s_i <- rep(0, n_m)
  
  parms <- list(
    n_m = n_m, n_b = n_b, n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
    c_prod = c_prod, c_cons = c_cons, c_nonmono_L = c_nonmono_L, c_nonmono_U = c_nonmono_U,
    w = w, K = K, s_i = s_i, d = d, c0 = c0, k = k
  )
  
  # Run simulation
  times <- seq(0, 10, by = 0.1)
  out <- ode(y = y0, times = times, func = gut_model, parms = parms)
  out_df <- as.data.frame(out)
  
  # Plot results
  out_long <- pivot_longer(out_df, -time, names_to = "variable", values_to = "value")
  
  # Add categories
  metabolite_subtypes <- c(rep("Produced", n_m_prod), rep("Consumed", n_m_cons), rep("Non-monotonic", n_m_nonmono))
  subtype_map <- setNames(metabolite_subtypes, paste0("M", 1:n_m))
  out_long$type <- ifelse(grepl("^B", out_long$variable), "Bacteria", "Metabolite")
  out_long$subtype <- subtype_map[out_long$variable]
  out_long$subtype[out_long$type == "Bacteria"] <- "Bacteria"
  
  # Bacteria plot
  p_bacteria <- ggplot(subset(out_long, type == "Bacteria"), aes(x = time, y = value, color = variable)) +
    geom_line(linewidth = 1) +
    labs(title = "Bacterial Dynamics", x = "Time", y = "Abundance") +
    theme_minimal()
  
  # Metabolite plot
  p_metabolites <- ggplot(subset(out_long, type == "Metabolite"), aes(x = time, y = value, color = variable)) +
    geom_line(linewidth = 1) +
    facet_wrap(~subtype, scales = "free_y") + 
    labs(title = "Metabolite Dynamics", x = "Time", y = "Concentration") +
    theme_minimal()
  
  print(p_bacteria)
  print(p_metabolites)
  
  return(list(data = out_df, parms = parms))
}


#==============================================================================
# MULTIPLE RUNS (10 SIMULATIONS)
#==============================================================================

run_multiple_sims <- function(n_runs = 10) {
  set.seed(42)
  
  n_b <- 3
  n_m_prod <- 3
  n_m_cons <- 7
  n_m_nonmono <- 2
  n_m <- n_m_prod + n_m_cons + n_m_nonmono
  
  c0 <- 0.1
  k <- 14
  d <- rep(0.1, n_b)
  K <- 2
  times <- seq(0, 10, by = 0.1)
  
  all_results <- list()
  
  for (r in 1:n_runs) {
    # New initial conditions each time
    B_init <- runif(n_b, 0.05, 0.2)
    M_prod_init <- runif(n_m_prod, exp(-3), exp(-2))
    M_cons_init <- runif(n_m_cons, 0.5, 1.0)
    M_nonmono_init <- runif(n_m_nonmono, exp(-3), exp(-2))
    y0 <- c(M_prod_init, M_cons_init, M_nonmono_init, B_init)
    names(y0) <- c(paste0("M", 1:n_m), paste0("B", 1:n_b))
    
    c_prod <- matrix(runif(n_m_prod * n_b, -1, 0), nrow = n_m_prod)
    c_cons <- matrix(runif(n_m_cons * n_b, 0, 1), nrow = n_m_cons)
    c_nonmono_L <- matrix(runif(n_m_nonmono * n_b, -1, 0), nrow = n_m_nonmono)
    c_nonmono_U <- matrix(runif(n_m_nonmono * n_b, 0, 1), nrow = n_m_nonmono)
    w <- matrix(runif(n_m_nonmono * n_b, 0.1, 1), nrow = n_m_nonmono)
    s_i <- rep(0, n_m)
    
    parms <- list(
      n_m = n_m, n_b = n_b, n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
      c_prod = c_prod, c_cons = c_cons, c_nonmono_L = c_nonmono_L, c_nonmono_U = c_nonmono_U,
      w = w, K = K, s_i = s_i, d = d, c0 = c0, k = k
    )
    
    out <- ode(y = y0, times = times, func = gut_model, parms = parms)
    out_df <- as.data.frame(out)
    out_df$run <- r
    all_results[[r]] <- out_df
  }
  
  all_results_df <- bind_rows(all_results)
  
  # Calculate mean and range
  summary_df <- all_results_df %>%
    pivot_longer(-c(time, run), names_to = "variable", values_to = "value") %>%
    group_by(time, variable) %>%
    summarise(
      mean_value = mean(value),
      min_value = min(value),
      max_value = max(value),
      .groups = "drop"
    )
  
  metabolite_subtypes <- c(rep("Produced", n_m_prod), rep("Consumed", n_m_cons), rep("Non-monotonic", n_m_nonmono))
  subtype_map <- setNames(metabolite_subtypes, paste0("M", 1:n_m))
  summary_df$subtype <- ifelse(grepl("^B", summary_df$variable), "Bacteria", subtype_map[summary_df$variable])
  
  # Plot with ribbons
  p_produced <- ggplot(subset(summary_df, subtype == "Produced"), 
                       aes(x = time, y = mean_value, color = variable, fill = variable)) +
    geom_ribbon(aes(ymin = min_value, ymax = max_value), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    labs(title = "Produced Metabolites (10 runs)", x = "Time", y = "Concentration") +
    theme_minimal()
  
  p_consumed <- ggplot(subset(summary_df, subtype == "Consumed"), 
                       aes(x = time, y = mean_value, color = variable, fill = variable)) +
    geom_ribbon(aes(ymin = min_value, ymax = max_value), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    labs(title = "Consumed Metabolites (10 runs)", x = "Time", y = "Concentration") +
    theme_minimal()
  
  p_nonmono <- ggplot(subset(summary_df, subtype == "Non-monotonic"), 
                      aes(x = time, y = mean_value, color = variable, fill = variable)) +
    geom_ribbon(aes(ymin = min_value, ymax = max_value), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    labs(title = "Non-monotonic Metabolites (10 runs)", x = "Time", y = "Concentration") +
    theme_minimal()
  
  print(p_produced)
  print(p_consumed)
  print(p_nonmono)
  
  return(all_results_df)
}

# DIFFERENT NUTRIENT SCENARIOS

compare_scenarios <- function() {
  set.seed(42)
  
  n_b <- 3
  n_m_prod <- 3
  n_m_cons <- 7
  n_m_nonmono <- 2
  n_m <- n_m_prod + n_m_cons + n_m_nonmono
  
  influx_levels <- c("Nutrient-Poor" = 0.01, "Balanced" = 0.05, "Nutrient-Rich" = 0.1)
  n_runs <- 10
  times <- seq(0, 20, by = 0.1)
  
  all_results <- list()
  
  for (scenario_name in names(influx_levels)) {
    influx <- influx_levels[scenario_name]
    
    for (run in 1:n_runs) {
      B_init <- runif(n_b, 0.05, 0.2)
      M_prod_init <- runif(n_m_prod, exp(-3), exp(-2))
      M_cons_init <- runif(n_m_cons, 0.5, 1.0)
      M_nonmono_init <- runif(n_m_nonmono, exp(-3), exp(-2))
      y0 <- c(M_prod_init, M_cons_init, M_nonmono_init, B_init)
      names(y0) <- c(paste0("M", 1:n_m), paste0("B", 1:n_b))
      
      c_prod <- matrix(runif(n_m_prod * n_b, -1, 0), nrow = n_m_prod)
      c_cons <- matrix(runif(n_m_cons * n_b, 0, 1), nrow = n_m_cons)
      c_nonmono_L <- matrix(runif(n_m_nonmono * n_b, -1, 0), nrow = n_m_nonmono)
      c_nonmono_U <- matrix(runif(n_m_nonmono * n_b, 0, 1), nrow = n_m_nonmono)
      w <- matrix(runif(n_m_nonmono * n_b, 0.1, 1), nrow = n_m_nonmono)
      
      s_i <- rep(0, n_m)
      s_i[(n_m_prod + 1):(n_m_prod + n_m_cons)] <- influx
      
      parms <- list(
        n_m = n_m, n_b = n_b, n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
        c_prod = c_prod, c_cons = c_cons, c_nonmono_L = c_nonmono_L, c_nonmono_U = c_nonmono_U,
        w = w, K = 2, s_i = s_i, d = rep(0.1, n_b), c0 = 0.1, k = 14
      )
      
      out <- ode(y = y0, times = times, func = gut_model, parms = parms)
      out_df <- as.data.frame(out)
      out_df$run <- run
      out_df$scenario <- scenario_name
      all_results[[length(all_results) + 1]] <- out_df
    }
  }
  
  results <- bind_rows(all_results)
  
  # Plot bacteria only
  plot_data <- results %>%
    select(time, run, scenario, B1, B2, B3) %>%
    pivot_longer(cols = starts_with("B"), names_to = "Bacteria", values_to = "Abundance") %>%
    mutate(scenario = factor(scenario, levels = c("Nutrient-Poor", "Balanced", "Nutrient-Rich")))
  
  summary_data <- plot_data %>%
    group_by(time, scenario, Bacteria) %>%
    summarise(
      mean_abundance = mean(Abundance),
      min_abundance = min(Abundance),
      max_abundance = max(Abundance),
      .groups = "drop"
    )
  
  p <- ggplot(summary_data, aes(x = time, y = mean_abundance, color = Bacteria, fill = Bacteria)) +
    geom_ribbon(aes(ymin = min_abundance, ymax = max_abundance), alpha = 0.3, color = NA) +
    geom_line(size = 1.2) +
    facet_wrap(~scenario, ncol = 3) +
    labs(title = "Bacterial Dynamics: Different Nutrient Levels",
         x = "Time", y = "Bacterial Abundance") +
    theme_minimal()
  
  print(p)
  return(results)
}


# BIFURCATION ANALYSIS

run_bifurcation <- function() {
  set.seed(42)
  
  n_b <- 3
  n_m_prod <- 3
  n_m_cons <- 7
  n_m_nonmono <- 2
  n_m <- n_m_prod + n_m_cons + n_m_nonmono
  
  # Base parameters
  c_prod <- matrix(runif(n_m_prod * n_b, -1, 0), nrow = n_m_prod)
  c_cons <- matrix(runif(n_m_cons * n_b, 0, 1), nrow = n_m_cons)
  c_nonmono_L <- matrix(runif(n_m_nonmono * n_b, -1, 0), nrow = n_m_nonmono)
  c_nonmono_U <- matrix(runif(n_m_nonmono * n_b, 0, 1), nrow = n_m_nonmono)
  w <- matrix(runif(n_m_nonmono * n_b, 0.1, 1), nrow = n_m_nonmono)
  s_i <- rep(0, n_m)
  
  base_parms <- list(
    n_m = n_m, n_b = n_b, n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
    c_prod = c_prod, c_cons = c_cons, c_nonmono_L = c_nonmono_L, c_nonmono_U = c_nonmono_U,
    w = w, K = 2, s_i = s_i, d = rep(0.1, n_b), c0 = 0.1, k = 14
  )
  
  # Initial state
  M0 <- numeric(n_m)
  M0[1:n_m_prod] <- runif(n_m_prod, exp(-5), exp(-3))
  M0[(n_m_prod + 1):(n_m_prod + n_m_cons)] <- runif(n_m_cons, 0.5, 1)
  M0[(n_m_prod + n_m_cons + 1):n_m] <- runif(n_m_nonmono, exp(-5), exp(-3))
  B0 <- runif(n_b, 0.05, 0.2)
  init_state <- c(M0, B0)
  
  # Vary c_upper for metabolite 11, bacterium 1
  met_idx <- 11
  bac_idx <- 1
  param_values <- seq(0, 1, length.out = 50)
  
  results <- data.frame()
  
  for (val in param_values) {
    parms <- base_parms
    parms$c_nonmono_U[met_idx - n_m_prod - n_m_cons, bac_idx] <- val
    
    # Find equilibrium
    times <- seq(0, 200, by = 0.1)
    out <- ode(y = init_state, times = times, func = gut_model, parms = parms)
    eq <- as.numeric(tail(out, 1)[, -1])
    eq[eq < 1e-10] <- 0
    
    result_row <- data.frame(
      param_value = val,
      B1 = eq[n_m + 1],
      B2 = eq[n_m + 2],
      B3 = eq[n_m + 3]
    )
    
    results <- rbind(results, result_row)
  }
  
  # Plot
  plot_data <- results %>%
    pivot_longer(cols = c(B1, B2, B3), names_to = "Bacterium", values_to = "Abundance")
  
  p <- ggplot(plot_data, aes(x = param_value, y = Abundance, fill = Bacterium)) +
    geom_area(alpha = 0.8) +
    scale_fill_manual(values = c("B1" = "#E57373", "B2" = "#81C784", "B3" = "#64B5F6")) +
    labs(title = "Bifurcation: varying c_upper for M11-B1",
         x = "c_upper", y = "Bacterial abundance at equilibrium") +
    theme_minimal()
  
  print(p)
  return(results)
}


# NETWORK VISUALIZATION

plot_network <- function() {
  set.seed(42)
  
  n_b <- 3
  n_m_prod <- 3
  n_m_cons <- 7
  n_m_nonmono <- 2
  n_m <- n_m_prod + n_m_cons + n_m_nonmono
  
  # Run a simulation first
  B_init <- runif(n_b, 0.05, 0.2)
  M_prod_init <- runif(n_m_prod, exp(-3), exp(-2))
  M_cons_init <- runif(n_m_cons, 0.5, 1.0)
  M_nonmono_init <- runif(n_m_nonmono, exp(-3), exp(-2))
  y0 <- c(M_prod_init, M_cons_init, M_nonmono_init, B_init)
  names(y0) <- c(paste0("M", 1:n_m), paste0("B", 1:n_b))
  
  influx_value <- 0.05
  s_i <- c(rep(influx_value, n_m_prod), rep(0, n_m_cons), rep(influx_value, n_m_nonmono))
  
  parms <- list(
    c_prod = matrix(runif(n_m_prod * n_b, -1, 0), nrow = n_m_prod),
    c_cons = matrix(runif(n_m_cons * n_b, 0, 1), nrow = n_m_cons),
    c_nonmono_L = matrix(runif(n_m_nonmono * n_b, -1, 0), nrow = n_m_nonmono),
    c_nonmono_U = matrix(runif(n_m_nonmono * n_b, 0, 1), nrow = n_m_nonmono),
    w = matrix(runif(n_m_nonmono * n_b, 0.1, 1), nrow = n_m_nonmono),
    K = 2, d = rep(0.1, n_b), c0 = 0.1, k = 14,
    n_m = n_m, n_b = n_b, n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
    s_i = s_i
  )
  
  times <- seq(0, 20, by = 0.5)
  out <- ode(y = y0, times = times, func = gut_model, parms = parms)
  sim_data <- as.data.frame(out)
  
  # Make network at t=2
  time_point <- 2.0
  state_at_t <- sim_data[sim_data$time == time_point, ]
  M <- as.numeric(state_at_t[1, 2:(n_m + 1)])
  B <- as.numeric(state_at_t[1, (n_m + 2):(n_m + n_b + 1)])
  
  # Build c_matrix at this time
  c_matrix <- rbind(parms$c_prod, parms$c_cons)
  c_nonmono_current <- matrix(0, nrow = n_m_nonmono, ncol = n_b)
  M_nonmono <- M[(n_m_prod + n_m_cons + 1):n_m]
  
  for (i in 1:n_m_nonmono) {
    for (j in 1:n_b) {
      ratio <- ifelse(B[j] > 1e-9, M_nonmono[i] / B[j], 0)
      hill_term <- (parms$w[i, j]^parms$K) / (parms$w[i, j]^parms$K + ratio^parms$K)
      c_nonmono_current[i, j] <- parms$c_nonmono_L[i, j] + 
        (parms$c_nonmono_U[i, j] - parms$c_nonmono_L[i, j]) * hill_term
    }
  }
  c_matrix <- rbind(c_matrix, c_nonmono_current)
  
  # Create network
  nodes <- data.frame(
    name = c(paste0("M", 1:n_m), paste0("B", 1:n_b)),
    type = c(rep("Metabolite", n_m), rep("Bacteria", n_b)),
    subtype = c(rep("Produced", n_m_prod), rep("Consumed", n_m_cons),
                rep("Non-monotonic", n_m_nonmono), rep("Bacteria", n_b))
  )
  
  edges <- expand.grid(from = paste0("M", 1:n_m), to = paste0("B", 1:n_b)) %>%
    mutate(
      m_idx = as.numeric(str_extract(from, "\\d+")),
      b_idx = as.numeric(str_extract(to, "\\d+")),
      cij = mapply(function(r, c) c_matrix[r, c], m_idx, b_idx),
      weight = abs(M[m_idx] * B[b_idx] * cij),
      interaction_type = ifelse(cij > 0, "Beneficial", "Inhibitory")
    ) %>%
    filter(weight > 1e-9 & M[m_idx] > 1e-4 & B[b_idx] > 1e-4)
  
  graph <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  
  p <- ggraph(graph, layout = 'bipartite') + 
    geom_edge_fan(aes(width = weight, alpha = weight, color = interaction_type), 
                  arrow = arrow(length = unit(2, 'mm'))) +
    geom_node_point(aes(color = subtype, shape = type), size = 10) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_width_continuous(range = c(0.2, 2)) +
    scale_edge_alpha_continuous(range = c(0.2, 1)) +
    scale_edge_color_manual(values = c("Beneficial" = "steelblue", "Inhibitory" = "firebrick3")) +
    scale_shape_manual(values = c("Metabolite" = 16, "Bacteria" = 15)) + 
    labs(title = paste("Network at time =", time_point)) +
    theme_graph()
  
  print(p)
  
  # Some network stats
  cat("\nNetwork Statistics:\n")
  cat("Number of nodes:", vcount(graph), "\n")
  cat("Number of edges:", ecount(graph), "\n")
  cat("Network density:", round(edge_density(graph), 4), "\n")
  
  return(graph)
}

# PLOT NETWORKS AT MULTIPLE TIME POINTS

plot_networks_at_times <- function(times_to_plot = c(0, 2.0, 20.0)) {
  set.seed(42)
  
  n_b <- 3
  n_m_prod <- 3
  n_m_cons <- 7
  n_m_nonmono <- 2
  n_m <- n_m_prod + n_m_cons + n_m_nonmono
  
  # Run simulation
  B_init <- runif(n_b, 0.05, 0.2)
  M_prod_init <- runif(n_m_prod, exp(-3), exp(-2))
  M_cons_init <- runif(n_m_cons, 0.5, 1.0)
  M_nonmono_init <- runif(n_m_nonmono, exp(-3), exp(-2))
  y0 <- c(M_prod_init, M_cons_init, M_nonmono_init, B_init)
  names(y0) <- c(paste0("M", 1:n_m), paste0("B", 1:n_b))
  
  influx_value <- 0.05
  s_i <- c(rep(influx_value, n_m_prod), rep(0, n_m_cons), rep(influx_value, n_m_nonmono))
  
  parms <- list(
    c_prod = matrix(runif(n_m_prod * n_b, -1, 0), nrow = n_m_prod),
    c_cons = matrix(runif(n_m_cons * n_b, 0, 1), nrow = n_m_cons),
    c_nonmono_L = matrix(runif(n_m_nonmono * n_b, -1, 0), nrow = n_m_nonmono),
    c_nonmono_U = matrix(runif(n_m_nonmono * n_b, 0, 1), nrow = n_m_nonmono),
    w = matrix(runif(n_m_nonmono * n_b, 0.1, 1), nrow = n_m_nonmono),
    K = 2, d = rep(0.1, n_b), c0 = 0.1, k = 14,
    n_m = n_m, n_b = n_b, n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
    s_i = s_i
  )
  
  times <- seq(0, 20, by = 0.5)
  out <- ode(y = y0, times = times, func = gut_model, parms = parms)
  sim_data <- as.data.frame(out)
  
  # Plot network at each time point
  for (t_point in times_to_plot) {
    state_at_t <- sim_data[sim_data$time == t_point, ]
    M <- as.numeric(state_at_t[1, 2:(n_m + 1)])
    B <- as.numeric(state_at_t[1, (n_m + 2):(n_m + n_b + 1)])
    
    # Build c_matrix
    c_matrix <- rbind(parms$c_prod, parms$c_cons)
    c_nonmono_current <- matrix(0, nrow = n_m_nonmono, ncol = n_b)
    M_nonmono <- M[(n_m_prod + n_m_cons + 1):n_m]
    
    for (i in 1:n_m_nonmono) {
      for (j in 1:n_b) {
        ratio <- ifelse(B[j] > 1e-9, M_nonmono[i] / B[j], 0)
        hill_term <- (parms$w[i, j]^parms$K) / (parms$w[i, j]^parms$K + ratio^parms$K)
        c_nonmono_current[i, j] <- parms$c_nonmono_L[i, j] + 
          (parms$c_nonmono_U[i, j] - parms$c_nonmono_L[i, j]) * hill_term
      }
    }
    c_matrix <- rbind(c_matrix, c_nonmono_current)
    
    # Create network
    nodes <- data.frame(
      name = c(paste0("M", 1:n_m), paste0("B", 1:n_b)),
      type = c(rep("Metabolite", n_m), rep("Bacteria", n_b)),
      subtype = c(rep("Produced", n_m_prod), rep("Consumed", n_m_cons),
                  rep("Non-monotonic", n_m_nonmono), rep("Bacteria", n_b))
    )
    
    edges <- expand.grid(from = paste0("M", 1:n_m), to = paste0("B", 1:n_b)) %>%
      mutate(
        m_idx = as.numeric(str_extract(from, "\\d+")),
        b_idx = as.numeric(str_extract(to, "\\d+")),
        cij = mapply(function(r, c) c_matrix[r, c], m_idx, b_idx),
        weight = abs(M[m_idx] * B[b_idx] * cij),
        interaction_type = ifelse(cij > 0, "Beneficial", "Inhibitory")
      ) %>%
      filter(weight > 1e-9 & M[m_idx] > 1e-4 & B[b_idx] > 1e-4)
    
    graph <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
    
    p <- ggraph(graph, layout = 'bipartite') + 
      geom_edge_fan(aes(width = weight, alpha = weight, color = interaction_type), 
                    arrow = arrow(length = unit(2, 'mm'))) +
      geom_node_point(aes(color = subtype, shape = type), size = 10) +
      geom_node_text(aes(label = name), repel = TRUE, size = 3) +
      scale_edge_width_continuous(range = c(0.2, 2)) +
      scale_edge_alpha_continuous(range = c(0.2, 1)) +
      scale_edge_color_manual(values = c("Beneficial" = "steelblue", "Inhibitory" = "firebrick3")) +
      scale_shape_manual(values = c("Metabolite" = 16, "Bacteria" = 15)) + 
      labs(title = paste("Interaction Network at Time =", t_point)) +
      theme_graph()
    
    print(p)
  }
  
  cat("\nNetworks plotted at times:", paste(times_to_plot, collapse = ", "), "\n")
}


# PLOT LARGE NETWORKS AT MULTIPLE TIME POINTS

plot_large_networks <- function(times_to_plot = c(0, 2.0, 20.0)) {
  set.seed(42)
  
  n_b <- 10
  n_m_prod <- 14
  n_m_cons <- 6
  n_m_nonmono <- 4
  n_m <- n_m_prod + n_m_cons + n_m_nonmono
  
  # Make it sparse
  sparsity <- 0.5
  make_sparse <- function(mat) {
    num_zeros <- round(length(mat) * sparsity)
    indices_to_zero <- sample(length(mat), num_zeros)
    mat[indices_to_zero] <- 0
    return(mat)
  }
  
  c_prod <- make_sparse(matrix(runif(n_m_prod * n_b, -1, 0), nrow = n_m_prod))
  c_cons <- make_sparse(matrix(runif(n_m_cons * n_b, 0, 1), nrow = n_m_cons))
  c_nonmono_L <- make_sparse(matrix(runif(n_m_nonmono * n_b, -1, 0), nrow = n_m_nonmono))
  c_nonmono_U <- make_sparse(matrix(runif(n_m_nonmono * n_b, 0, 1), nrow = n_m_nonmono))
  w <- matrix(runif(n_m_nonmono * n_b, 0.1, 1), nrow = n_m_nonmono)
  
  influx_value <- 0.05
  s_i <- c(rep(influx_value, n_m_prod), rep(0, n_m_cons), rep(influx_value, n_m_nonmono))
  
  parms <- list(
    n_m = n_m, n_b = n_b, n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
    c_prod = c_prod, c_cons = c_cons, c_nonmono_L = c_nonmono_L, c_nonmono_U = c_nonmono_U,
    w = w, K = 2, s_i = s_i, d = rep(0.1, n_b), c0 = 0.1, k = 14
  )
  
  y0 <- c(runif(n_m_prod, exp(-3), exp(-2)), 
          runif(n_m_cons, 0.5, 1.0),
          runif(n_m_nonmono, exp(-3), exp(-2)), 
          runif(n_b, 0.05, 0.2))
  names(y0) <- c(paste0("M", 1:n_m), paste0("B", 1:n_b))
  
  times <- seq(0, 20, by = 0.5)
  out <- ode(y = y0, times = times, func = gut_model, parms = parms)
  sim_data <- as.data.frame(out)
  
  # Plot network at each time point
  for (t_point in times_to_plot) {
    state_at_t <- sim_data[sim_data$time == t_point, ]
    M <- as.numeric(state_at_t[1, 2:(n_m + 1)])
    B <- as.numeric(state_at_t[1, (n_m + 2):(n_m + n_b + 1)])
    
    # Build c_matrix
    c_matrix <- rbind(parms$c_prod, parms$c_cons)
    c_nonmono_current <- matrix(0, nrow = n_m_nonmono, ncol = n_b)
    M_nonmono <- M[(n_m_prod + n_m_cons + 1):n_m]
    
    for (i in 1:n_m_nonmono) {
      for (j in 1:n_b) {
        ratio <- ifelse(B[j] > 1e-9, M_nonmono[i] / B[j], 0)
        hill_term <- (parms$w[i, j]^parms$K) / (parms$w[i, j]^parms$K + ratio^parms$K)
        c_nonmono_current[i, j] <- parms$c_nonmono_L[i, j] + 
          (parms$c_nonmono_U[i, j] - parms$c_nonmono_L[i, j]) * hill_term
      }
    }
    c_matrix <- rbind(c_matrix, c_nonmono_current)
    
    # Create network
    nodes <- data.frame(
      name = c(paste0("M", 1:n_m), paste0("B", 1:n_b)),
      type = c(rep("Metabolite", n_m), rep("Bacteria", n_b)),
      subtype = c(rep("Produced", n_m_prod), rep("Consumed", n_m_cons),
                  rep("Non-monotonic", n_m_nonmono), rep("Bacteria", n_b))
    )
    
    edges <- expand.grid(from = paste0("M", 1:n_m), to = paste0("B", 1:n_b)) %>%
      mutate(
        m_idx = as.numeric(str_extract(from, "\\d+")),
        b_idx = as.numeric(str_extract(to, "\\d+")),
        cij = mapply(function(r, c) c_matrix[r, c], m_idx, b_idx),
        weight = abs(M[m_idx] * B[b_idx] * cij),
        interaction_type = ifelse(cij > 0, "Beneficial", "Inhibitory")
      ) %>%
      filter(weight > 1e-6)
    
    graph <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
    
    p <- ggraph(graph, layout = 'nicely') + 
      geom_edge_fan(aes(width = weight, alpha = weight, color = interaction_type), 
                    arrow = arrow(length = unit(2, 'mm'))) +
      geom_node_point(aes(color = subtype, shape = type), size = 8) +
      geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
      scale_edge_width_continuous(range = c(0.2, 1.5)) +
      scale_edge_alpha_continuous(range = c(0.3, 1)) +
      scale_edge_color_manual(values = c("Beneficial" = "steelblue", "Inhibitory" = "firebrick3")) +
      scale_shape_manual(values = c("Metabolite" = 16, "Bacteria" = 15)) + 
      labs(title = paste("Large Network at Time =", t_point)) +
      theme_graph()
    
    print(p)
  }
  
  cat("\nLarge networks plotted at times:", paste(times_to_plot, collapse = ", "), "\n")
}


# NETWORK ANALYSIS WITH STATS

network_analysis <- function() {
  set.seed(42)
  
  n_b <- 3
  n_m_prod <- 3
  n_m_cons <- 7
  n_m_nonmono <- 2
  n_m <- n_m_prod + n_m_cons + n_m_nonmono
  
  # Run simulation
  B_init <- runif(n_b, 0.05, 0.2)
  M_prod_init <- runif(n_m_prod, exp(-3), exp(-2))
  M_cons_init <- runif(n_m_cons, 0.5, 1.0)
  M_nonmono_init <- runif(n_m_nonmono, exp(-3), exp(-2))
  y0 <- c(M_prod_init, M_cons_init, M_nonmono_init, B_init)
  names(y0) <- c(paste0("M", 1:n_m), paste0("B", 1:n_b))
  
  influx_value <- 0.05
  s_i <- c(rep(influx_value, n_m_prod), rep(0, n_m_cons), rep(influx_value, n_m_nonmono))
  
  parms <- list(
    c_prod = matrix(runif(n_m_prod * n_b, -1, 0), nrow = n_m_prod),
    c_cons = matrix(runif(n_m_cons * n_b, 0, 1), nrow = n_m_cons),
    c_nonmono_L = matrix(runif(n_m_nonmono * n_b, -1, 0), nrow = n_m_nonmono),
    c_nonmono_U = matrix(runif(n_m_nonmono * n_b, 0, 1), nrow = n_m_nonmono),
    w = matrix(runif(n_m_nonmono * n_b, 0.1, 1), nrow = n_m_nonmono),
    K = 2, d = rep(0.1, n_b), c0 = 0.1, k = 14,
    n_m = n_m, n_b = n_b, n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
    s_i = s_i
  )
  
  times <- seq(0, 20, by = 0.5)
  out <- ode(y = y0, times = times, func = gut_model, parms = parms)
  sim_data <- as.data.frame(out)
  
  # Build network at t=2
  time_point <- 2.0
  state_at_t <- sim_data[sim_data$time == time_point, ]
  M <- as.numeric(state_at_t[1, 2:(n_m + 1)])
  B <- as.numeric(state_at_t[1, (n_m + 2):(n_m + n_b + 1)])
  
  c_matrix <- rbind(parms$c_prod, parms$c_cons)
  c_nonmono_current <- matrix(0, nrow = n_m_nonmono, ncol = n_b)
  M_nonmono <- M[(n_m_prod + n_m_cons + 1):n_m]
  
  for (i in 1:n_m_nonmono) {
    for (j in 1:n_b) {
      ratio <- ifelse(B[j] > 1e-9, M_nonmono[i] / B[j], 0)
      hill_term <- (parms$w[i, j]^parms$K) / (parms$w[i, j]^parms$K + ratio^parms$K)
      c_nonmono_current[i, j] <- parms$c_nonmono_L[i, j] + 
        (parms$c_nonmono_U[i, j] - parms$c_nonmono_L[i, j]) * hill_term
    }
  }
  c_matrix <- rbind(c_matrix, c_nonmono_current)
  
  nodes <- data.frame(
    name = c(paste0("M", 1:n_m), paste0("B", 1:n_b)),
    type = c(rep(FALSE, n_m), rep(TRUE, n_b))
  )
  
  edges <- expand.grid(from = paste0("M", 1:n_m), to = paste0("B", 1:n_b)) %>%
    mutate(
      m_idx = as.numeric(str_extract(from, "\\d+")),
      b_idx = as.numeric(str_extract(to, "\\d+")),
      cij = mapply(function(r, c) c_matrix[r, c], m_idx, b_idx),
      weight = abs(M[m_idx] * B[b_idx] * cij)
    ) %>%
    filter(weight > 1e-9)
  
  graph <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  g_undirected <- as.undirected(graph)
  
  # Network metrics
  cat("\n=== Network Analysis ===\n")
  cat("Nodes:", vcount(graph), "\n")
  cat("Edges:", ecount(graph), "\n")
  cat("Density:", round(edge_density(graph), 4), "\n")
  cat("Average Degree:", round(mean(degree(graph)), 2), "\n")
  cat("Diameter:", diameter(g_undirected, weights = NA), "\n")
  cat("Average Path Length:", round(mean_distance(g_undirected), 3), "\n")
  
  # Centrality
  centrality_df <- data.frame(
    Node = V(graph)$name,
    Degree = degree(graph),
    Betweenness = betweenness(graph, weights = E(graph)$weight),
    Closeness = closeness(graph, weights = E(graph)$weight)
  )
  
  cat("\nTop 5 Nodes by Degree:\n")
  print(head(centrality_df[order(-centrality_df$Degree), c("Node", "Degree")], 5))
  
  # Clustering
  cat("\nClustering Coefficient:", round(transitivity(graph, type = "global"), 4), "\n")
  
  # Communities
  communities <- cluster_fast_greedy(g_undirected, weights = E(g_undirected)$weight)
  cat("Modularity:", round(modularity(communities), 4), "\n")
  cat("Number of communities:", length(communities), "\n")
  
  # Similarity networks
  W <- outer(M, B) * c_matrix
  bacteria_similarity <- t(W) %*% W
  metabolite_similarity <- W %*% t(W)
  
  g_bac_sim <- graph_from_adjacency_matrix(bacteria_similarity, mode = "undirected", weighted = TRUE)
  g_met_sim <- graph_from_adjacency_matrix(metabolite_similarity, mode = "undirected", weighted = TRUE)
  
  threshold <- 1e-6
  g_bac_sim <- delete_edges(g_bac_sim, E(g_bac_sim)[weight < threshold])
  g_met_sim <- delete_edges(g_met_sim, E(g_met_sim)[weight < threshold])
  
  cat("\nBacteria-Bacteria Similarity:\n")
  cat("  Nodes:", vcount(g_bac_sim), "\n")
  cat("  Edges:", ecount(g_bac_sim), "\n")
  cat("  Density:", round(edge_density(g_bac_sim), 4), "\n")
  
  cat("\nMetabolite-Metabolite Similarity:\n")
  cat("  Nodes:", vcount(g_met_sim), "\n")
  cat("  Edges:", ecount(g_met_sim), "\n")
  cat("  Density:", round(edge_density(g_met_sim), 4), "\n")
  
  # Plot similarity networks
  par(mfrow=c(1,2))
  
  plot(g_bac_sim,
       main = "Bacteria-Bacteria Similarity",
       vertex.color = "orange",
       vertex.label.color = "black",
       vertex.size = 30,
       edge.width = E(g_bac_sim)$weight * 5,
       edge.color = "grey60")
  
  plot(g_met_sim,
       main = "Metabolite-Metabolite Similarity",
       vertex.color = "skyblue",
       vertex.label.color = "black",
       vertex.size = 10,
       edge.width = E(g_met_sim)$weight * 2,
       edge.color = "grey60")
  
  par(mfrow=c(1,1))
  
  return(list(main = graph, bacteria_sim = g_bac_sim, metabolite_sim = g_met_sim))
}


# STABILITY ANALYSIS

stability_analysis <- function() {
  set.seed(42)
  
  n_b <- 3
  n_m_prod <- 3
  n_m_cons <- 7
  n_m_nonmono <- 2
  n_m <- n_m_prod + n_m_cons + n_m_nonmono
  
  # Setup parameters
  c_prod <- matrix(runif(n_m_prod * n_b, -1, 0), nrow = n_m_prod)
  c_cons <- matrix(runif(n_m_cons * n_b, 0, 1), nrow = n_m_cons)
  c_nonmono_L <- matrix(runif(n_m_nonmono * n_b, -1, 0), nrow = n_m_nonmono)
  c_nonmono_U <- matrix(runif(n_m_nonmono * n_b, 0, 1), nrow = n_m_nonmono)
  w <- matrix(runif(n_m_nonmono * n_b, 0.1, 1), nrow = n_m_nonmono)
  
  # Initial state
  M0 <- numeric(n_m)
  M0[1:n_m_prod] <- runif(n_m_prod, exp(-3), exp(-2))
  M0[(n_m_prod + 1):(n_m_prod + n_m_cons)] <- runif(n_m_cons, 0.5, 1.0)
  M0[(n_m_prod + n_m_cons + 1):n_m] <- runif(n_m_nonmono, exp(-3), exp(-2))
  B0 <- runif(n_b, 0.05, 0.2)
  y0 <- c(M0, B0)
  
  # Closed system (no influx)
  s_closed <- rep(0, n_m)
  parms_closed <- list(
    n_m = n_m, n_b = n_b, n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
    c_prod = c_prod, c_cons = c_cons, c_nonmono_L = c_nonmono_L, c_nonmono_U = c_nonmono_U,
    w = w, K = 2, s_i = s_closed, d = rep(0.1, n_b), c0 = 0.1, k = 14
  )
  
  # Open system (with influx)
  s_open <- rep(0, n_m)
  s_open[(n_m_prod + 1):(n_m_prod + n_m_cons)] <- 0.05
  parms_open <- list(
    n_m = n_m, n_b = n_b, n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
    c_prod = c_prod, c_cons = c_cons, c_nonmono_L = c_nonmono_L, c_nonmono_U = c_nonmono_U,
    w = w, K = 2, s_i = s_open, d = rep(0.1, n_b), c0 = 0.1, k = 14
  )
  
  # Find equilibria
  rhs <- function(y, p) {
    as.numeric(gut_model(0, y, p)[[1]])
  }
  
  eq_closed <- multiroot(f = function(y) rhs(y, parms_closed), start = y0)$root
  eq_open <- multiroot(f = function(y) rhs(y, parms_open), start = y0)$root
  
  # Compute Jacobians
  J_closed <- jacobian(func = function(y) rhs(y, parms_closed), x = eq_closed)
  J_open <- jacobian(func = function(y) rhs(y, parms_open), x = eq_open)
  
  # Get eigenvalues
  eig_closed <- eigen(J_closed)$values
  eig_open <- eigen(J_open)$values
  
  # Print results
  cat("\n=== CLOSED SYSTEM ===\n")
  cat("Dominant eigenvalue:", round(max(Re(eig_closed)), 6), "\n")
  if (max(Re(eig_closed)) < 0) {
    cat("Status: STABLE\n")
  } else {
    cat("Status: UNSTABLE\n")
  }
  
  cat("\n=== OPEN SYSTEM ===\n")
  cat("Dominant eigenvalue:", round(max(Re(eig_open)), 6), "\n")
  if (max(Re(eig_open)) < 0) {
    cat("Status: STABLE\n")
  } else {
    cat("Status: UNSTABLE\n")
  }
  
  # Plot eigenvalues
  par(mfrow = c(1, 2))
  plot(Re(eig_closed), Im(eig_closed), pch = 19, col = "red",
       main = "Closed System", xlab = "Real(λ)", ylab = "Imag(λ)")
  abline(v = 0, lty = 2)
  
  plot(Re(eig_open), Im(eig_open), pch = 19, col = "blue",
       main = "Open System", xlab = "Real(λ)", ylab = "Imag(λ)")
  abline(v = 0, lty = 2)
  par(mfrow = c(1, 1))
  
  return(list(closed = eig_closed, open = eig_open))
}


# STABILITY ANALYSIS - LARGE SYSTEM

stability_large <- function() {
  set.seed(123)
  
  # Large system dimensions
  n_b <- 10
  n_m_cons <- 14
  n_m_prod <- 6
  n_m_nonmono <- 4
  n_m <- n_m_cons + n_m_prod + n_m_nonmono
  
  # Initial conditions
  M0 <- numeric(n_m)
  M0[1:n_m_cons] <- runif(n_m_cons, 0.5, 1)
  M0[(n_m_cons + 1):(n_m_cons + n_m_prod)] <- runif(n_m_prod, exp(-3), exp(-2))
  M0[(n_m_cons + n_m_prod + 1):n_m] <- runif(n_m_nonmono, exp(-3), exp(-2))
  B0 <- runif(n_b, 0.05, 0.2)
  y0 <- c(M0, B0)
  
  # Interaction matrices - note ordering: cons, prod, nonmono
  c_cons <- matrix(runif(n_m_cons * n_b, 0, 1), nrow = n_m_cons)
  c_prod <- matrix(runif(n_m_prod * n_b, -1, 0), nrow = n_m_prod)
  c_nonmono_L <- matrix(runif(n_m_nonmono * n_b, -1, 0), nrow = n_m_nonmono)
  c_nonmono_U <- matrix(runif(n_m_nonmono * n_b, 0, 1), nrow = n_m_nonmono)
  w <- matrix(runif(n_m_nonmono * n_b, 0.1, 1), nrow = n_m_nonmono)
  
  # Parameters for closed system
  s_closed <- rep(0, n_m)
  parms_closed <- list(
    n_m = n_m, n_b = n_b, 
    n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
    c_prod = c_prod, c_cons = c_cons,
    c_nonmono_L = c_nonmono_L, c_nonmono_U = c_nonmono_U,
    w = w, K = 2, s_i = s_closed,
    d = runif(n_b, 0.05, 0.1), c0 = 0.5, k = 1.0
  )
  
  # Parameters for open system
  s_open <- rep(0, n_m)
  s_open[1:n_m_cons] <- runif(n_m_cons, 0.05, 0.1)
  parms_open <- list(
    n_m = n_m, n_b = n_b,
    n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
    c_prod = c_prod, c_cons = c_cons,
    c_nonmono_L = c_nonmono_L, c_nonmono_U = c_nonmono_U,
    w = w, K = 2, s_i = s_open,
    d = runif(n_b, 0.05, 0.1), c0 = 0.5, k = 1.0
  )
  
  # Find equilibria
  rhs <- function(y, p) {
    as.numeric(gut_model(0, y, p)[[1]])
  }
  
  eq_closed <- multiroot(f = function(y) rhs(y, parms_closed), start = y0)$root
  eq_open <- multiroot(f = function(y) rhs(y, parms_open), start = y0)$root
  
  # Compute Jacobians
  J_closed <- jacobian(func = function(y) rhs(y, parms_closed), x = eq_closed)
  J_open <- jacobian(func = function(y) rhs(y, parms_open), x = eq_open)
  
  # Get eigenvalues
  eig_closed <- eigen(J_closed)$values
  eig_open <- eigen(J_open)$values
  
  # Calculate stability indices
  S_closed <- mean(Re(eig_closed))
  S_open <- mean(Re(eig_open))
  
  # Print results
  cat("\n=== LARGE SYSTEM: CLOSED ===\n")
  cat("Dominant eigenvalue:", round(max(Re(eig_closed)), 6), "\n")
  cat("Stability index:", round(S_closed, 6), "\n")
  if (max(Re(eig_closed)) < 0) {
    cat("Status: STABLE\n")
  } else {
    cat("Status: UNSTABLE\n")
  }
  
  cat("\n=== LARGE SYSTEM: OPEN ===\n")
  cat("Dominant eigenvalue:", round(max(Re(eig_open)), 6), "\n")
  cat("Stability index:", round(S_open, 6), "\n")
  if (max(Re(eig_open)) < 0) {
    cat("Status: STABLE\n")
  } else {
    cat("Status: UNSTABLE\n")
  }
  
  # Plot eigenvalues
  par(mfrow = c(1, 2))
  plot(Re(eig_closed), Im(eig_closed), pch = 19, col = "red",
       main = "Large System: Closed", xlab = "Real(λ)", ylab = "Imag(λ)")
  abline(v = 0, lty = 2)
  
  plot(Re(eig_open), Im(eig_open), pch = 19, col = "blue",
       main = "Large System: Open", xlab = "Real(λ)", ylab = "Imag(λ)")
  abline(v = 0, lty = 2)
  par(mfrow = c(1, 1))
  
  return(list(closed = eig_closed, open = eig_open))
}


# LARGE NETWORK 

large_network <- function() {
  set.seed(42)
  
  n_b <- 10
  n_m_prod <- 14
  n_m_cons <- 6
  n_m_nonmono <- 4
  n_m <- n_m_prod + n_m_cons + n_m_nonmono
  
  # Make it sparse (50% of interactions are zero)
  sparsity <- 0.5
  
  make_sparse <- function(mat) {
    num_zeros <- round(length(mat) * sparsity)
    indices_to_zero <- sample(length(mat), num_zeros)
    mat[indices_to_zero] <- 0
    return(mat)
  }
  
  c_prod <- make_sparse(matrix(runif(n_m_prod * n_b, -1, 0), nrow = n_m_prod))
  c_cons <- make_sparse(matrix(runif(n_m_cons * n_b, 0, 1), nrow = n_m_cons))
  c_nonmono_L <- make_sparse(matrix(runif(n_m_nonmono * n_b, -1, 0), nrow = n_m_nonmono))
  c_nonmono_U <- make_sparse(matrix(runif(n_m_nonmono * n_b, 0, 1), nrow = n_m_nonmono))
  w <- matrix(runif(n_m_nonmono * n_b, 0.1, 1), nrow = n_m_nonmono)
  
  influx_value <- 0.05
  s_i <- c(rep(influx_value, n_m_prod), rep(0, n_m_cons), rep(influx_value, n_m_nonmono))
  
  parms <- list(
    n_m = n_m, n_b = n_b, n_m_prod = n_m_prod, n_m_cons = n_m_cons, n_m_nonmono = n_m_nonmono,
    c_prod = c_prod, c_cons = c_cons, c_nonmono_L = c_nonmono_L, c_nonmono_U = c_nonmono_U,
    w = w, K = 2, s_i = s_i, d = rep(0.1, n_b), c0 = 0.1, k = 14
  )
  
  y0 <- c(runif(n_m_prod, exp(-3), exp(-2)), 
          runif(n_m_cons, 0.5, 1.0),
          runif(n_m_nonmono, exp(-3), exp(-2)), 
          runif(n_b, 0.05, 0.2))
  names(y0) <- c(paste0("M", 1:n_m), paste0("B", 1:n_b))
  
  times <- seq(0, 20, by = 0.5)
  out <- ode(y = y0, times = times, func = gut_model, parms = parms)
  
  cat("Large network simulation complete!\n")
  cat("Number of metabolites:", n_m, "\n")
  cat("Number of bacteria:", n_b, "\n")
  cat("Sparsity:", sparsity, "\n")
  
  return(as.data.frame(out))
}