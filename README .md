# Gut Microbiome Model

Combined R script for simulating gut microbiome dynamics with bacterial species and metabolites.

## Requirements

```r
library(deSolve)
library(tidyverse)
library(igraph)
library(ggraph)
library(patchwork)
library(bipartite)
library(rootSolve)
library(numDeriv)
```

## Model Description

The model simulates interactions between bacteria and metabolites in the gut, including:
- Produced metabolites (negative interactions with bacteria)
- Consumed metabolites (positive interactions with bacteria)
- Non-monotonic metabolites (Hill function switching between production and consumption)

## Available Functions

### Basic Simulations
- `run_single_simulation()` - Run one simulation and plot results
- `run_multiple_sims(n_runs = 10)` - Run multiple simulations and show mean ± range
- `compare_scenarios()` - Compare different nutrient influx levels

### Analysis
- `run_bifurcation()` - Bifurcation analysis varying interaction parameters
- `plot_network()` - Visualize single metabolite-bacteria interaction network
- `plot_networks_at_times()` - Plot small networks at multiple time points (default: t=0, 2, 20)
- `plot_large_networks()` - Plot large networks at multiple time points
- `network_analysis()` - Full network analysis with centrality, clustering, and similarity networks
- `stability_analysis()` - Linear stability analysis with eigenvalues (small system)
- `stability_large()` - Linear stability analysis for large system
- `large_network()` - Simulate larger system (10 bacteria, 24 metabolites)

## Usage

Load the script:
```r
source("combined_gut_model.R")
```

Run a basic simulation:
```r
result <- run_single_simulation()
```

Run multiple simulations:
```r
all_runs <- run_multiple_sims(10)
```

Compare nutrient scenarios:
```r
scenarios <- compare_scenarios()
```

Perform bifurcation analysis:
```r
bifurcation_data <- run_bifurcation()
```

Visualize network:
```r
network <- plot_network()
```

Plot networks at multiple times (small system):
```r
plot_networks_at_times(c(0, 2, 20))
```

Plot large networks:
```r
plot_large_networks(c(0, 2, 20))
```

Full network analysis:
```r
networks <- network_analysis()
```

Check stability (small system):
```r
eigenvalues <- stability_analysis()
```

Check stability (large system):
```r
eigenvalues_large <- stability_large()
```

## Model Parameters

Default system size:
- 3 bacterial species
- 3 produced metabolites
- 7 consumed metabolites  
- 2 non-monotonic metabolites

Key parameters:
- `c0 = 0.1` - Conversion coefficient
- `k = 14` - Carrying capacity
- `d = 0.1` - Death rate
- `K = 2` - Hill function exponent

## Output

Functions produce ggplot2 visualizations showing:
- Time series of bacterial abundances
- Time series of metabolite concentrations
- Bifurcation diagrams
- Interaction networks
- Eigenvalue spectra
