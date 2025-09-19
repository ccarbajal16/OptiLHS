# OptiLHS Technical Guide

<div align="center">
  <img src="inst/shiny-app/www/hex-Opti-LHS.png" alt="OptiLHS Logo" width="200" height="200">
</div>

**OptiLHS:** Optimized Latin Hypercube Sampling for Soil Sampling Designs

---
**Version:** 0.1.1  
**Author:** Carlos Carbajal  
**Date:** September 2025

## Table of Contents

1. [Overview](#overview)
2. [Module 1: KL Divergence Sample Size Optimization](#module-1-kl-divergence-sample-size-optimization)
3. [Module 2: Random Forest Optimization](#module-2-random-forest-optimization)
4. [Module 3: Alternative Site Selection](#module-3-alternative-site-selection)
5. [Utility Functions](#utility-functions)
6. [Mathematical Background](#mathematical-background)
7. [Implementation Examples](#implementation-examples)
8. [Performance Considerations](#performance-considerations)

## Overview

The `OptiLHS` package implements a comprehensive three-module workflow for optimizing soil sampling designs using conditioned Latin Hypercube Sampling (cLHS) methods. This technical guide provides detailed documentation of all functions, algorithms, and mathematical foundations.

### Core Methodology

The package is built around three main optimization approaches:

1. **KL Divergence Optimization** (Malone et al. 2019) - Determines optimal sample sizes
2. **Random Forest Enhancement** - Improves sampling designs using ML-based optimization
3. **Alternative Site Selection** - Handles inaccessible locations through buffer-based replacement

## Module 1: KL Divergence Sample Size Optimization

### Primary Function: `optimize_sample_size()`

**Purpose:** Determines the optimal number of samples using Kullback-Leibler divergence analysis.

**Mathematical Foundation:**
- Implements the Malone et al. (2019) methodology
- Fits exponential decay function: `KL = b1 * exp(-k * n) + b0`
- Uses cumulative density function to identify optimal sample size

#### Function Signature

```r
optimize_sample_size(
  population_data,
  min_samples = 10,
  max_samples = 500,
  step_size = 10,
  n_replicates = 10,
  n_bins = 25,
  probability_threshold = 0.95
)
```

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `population_data` | data.frame | - | Environmental covariate data for the entire study area |
| `min_samples` | integer | 10 | Minimum number of samples to test |
| `max_samples` | integer | 500 | Maximum number of samples to test |
| `step_size` | integer | 10 | Increment size for sample size sequence |
| `n_replicates` | integer | 10 | Number of replicate samplings per sample size |
| `n_bins` | integer | 25 | Number of histogram bins for KL divergence calculation |
| `probability_threshold` | numeric | 0.95 | CDF threshold for optimal size determination |

#### Return Value

Returns a list with the following components:

```r
list(
  raw_results = data.frame,      # All KL divergence results
  summary_results = data.frame,  # Mean KL by sample size
  fitted_model = nls,           # Exponential decay model
  fitted_curve = data.frame,    # Model predictions
  optimal_sample_size = integer, # Recommended sample size
  plot_kl = ggplot,            # KL divergence plot
  plot_cdf = ggplot            # Cumulative density plot
)
```

#### Algorithm Details

1. **Sample Size Testing Loop:**
   ```r
   for (n_samples in sample_sizes) {
     for (rep in 1:n_replicates) {
       clhs_result <- clhs::clhs(population_data, size = n_samples)
       kl_div <- calculate_kl_divergence(population_data, sample_data)
     }
   }
   ```

2. **Exponential Model Fitting:**
   ```r
   exp_model <- nlsLM(
     mean_kl ~ b1 * exp(-k * sample_size) + b0,
     start = list(b0 = min_kl, b1 = range_kl, k = 0.01)
   )
   ```

3. **Optimal Size Determination:**
   ```r
   cdf_values <- (max_kl - fitted_kl) / max_improvement
   optimal_idx <- which(cdf_values >= probability_threshold)[1]
   ```

### Supporting Function: `calculate_kl_divergence()`

**Purpose:** Computes Kullback-Leibler divergence between sample and population distributions.

#### Function Signature

```r
calculate_kl_divergence(population_data, sample_data, n_bins = 25)
```

#### Mathematical Implementation

```r
# For each continuous variable:
# 1. Create histogram bins
breaks <- seq(pop_range[1], pop_range[2], length.out = n_bins + 1)

# 2. Calculate population distribution (Ei)
pop_density <- pop_hist$counts / sum(pop_hist$counts)

# 3. Calculate sample distribution (Oi)  
sample_density <- sample_hist$counts / sum(sample_hist$counts)

# 4. Compute KL divergence: KL = Σ Oi * log(Oi/Ei)
kl <- sum(sample_density * (log(sample_density) - log(pop_density)))
```

**Key Features:**
- Handles multiple environmental variables
- Prevents log(0) errors with small constant addition
- Returns mean KL divergence across all variables

## Module 2: Random Forest Optimization

### Primary Function: `rf_optimize_sampling()`

**Purpose:** Enhances cLHS designs using Random Forest-based simulated annealing optimization.

#### Function Signature

```r
rf_optimize_sampling(
  covariates,
  n_samples,
  n_iterations = 500,
  seed = NULL
)
```

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `covariates` | SpatRaster | - | Environmental raster layers |
| `n_samples` | integer | - | Number of samples to optimize |
| `n_iterations` | integer | 500 | Number of optimization iterations |
| `seed` | integer | NULL | Random seed for reproducibility |

#### Return Value

```r
list(
  optimized_samples = data.frame,  # Final optimized sample locations
  initial_samples = data.frame,    # Initial cLHS sample locations
  initial_mse = numeric,          # Initial Random Forest MSE
  final_mse = numeric,            # Final optimized MSE
  improvement = numeric           # Percentage improvement
)
```

#### Algorithm: Simulated Annealing

1. **Initialization:**
   ```r
   initial_samples <- simple_clhs_sampling(covariates, n_samples, seed)
   current_mse <- calculate_rf_mse(initial_samples, full_data)
   temperature <- 1000
   cooling_rate <- 0.95
   ```

2. **Optimization Loop:**
   ```r
   for (i in 1:n_iterations) {
     # Generate candidate by swapping random point
     candidate_sample <- current_sample
     replace_idx <- sample(nrow(candidate_sample), 1)
     candidate_sample[replace_idx, ] <- random_replacement_point
     
     # Calculate MSE and acceptance probability
     candidate_mse <- calculate_rf_mse(candidate_sample, full_data)
     delta <- candidate_mse - current_mse
     
     # Accept/reject based on simulated annealing criterion
     if (delta < 0 || runif(1) < exp(-delta / temperature)) {
       current_sample <- candidate_sample
       current_mse <- candidate_mse
     }
     
     # Cool down temperature
     temperature <- temperature * cooling_rate
   }
   ```

### Supporting Function: `calculate_rf_mse()`

**Purpose:** Evaluates sampling design quality using Random Forest cross-validation MSE.

#### Implementation Details

```r
calculate_rf_mse(sample_points, full_data, n_folds = 5)
```

**Algorithm:**
1. **Target Variable Selection:** Uses third environmental variable as target
2. **Cross-Validation:** 5-fold CV to estimate prediction error
3. **Random Forest Training:** 100 trees per fold
4. **MSE Calculation:** Mean squared error across all folds

### Supporting Function: `simple_clhs_sampling()`

**Purpose:** Generates initial cLHS design for optimization.

#### Key Features

```r
simple_clhs_sampling(covariates, n_samples, seed = NULL)
```

- Converts SpatRaster to data.frame with coordinates
- Removes constant variables (variance < 1e-10)
- Handles missing data through `na.omit()`
- Returns sample with x,y coordinates + environmental values

## Module 3: Alternative Site Selection

### Primary Function: `process_inaccessible_sites()`

**Purpose:** Finds alternative sampling locations for inaccessible sites using buffer-based environmental similarity.

#### Function Signature

```r
process_inaccessible_sites(
  inaccessible_sites,
  raster_data,
  buffer_distance = 250,
  similarity_threshold = 0.90,
  categorical_vars = NULL,
  method = "mahalanobis",
  select_method = "random"
)
```

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `inaccessible_sites` | data.frame | - | Coordinates (x,y) of inaccessible locations |
| `raster_data` | SpatRaster | - | Environmental covariate layers |
| `buffer_distance` | numeric | 250 | Buffer radius in map units |
| `similarity_threshold` | numeric | 0.90 | Minimum environmental similarity |
| `categorical_vars` | character | NULL | Names of categorical variables |
| `method` | character | "mahalanobis" | Distance method for similarity |
| `select_method` | character | "random" | Alternative selection method |

#### Return Value

```r
list(
  results = list,              # Detailed results per site
  alternative_sites = data.frame, # All alternative locations
  success_rate = numeric       # Proportion of successful replacements
)
```

#### Algorithm Workflow

1. **Buffer Creation:**
   ```r
   for (each_inaccessible_site) {
     buffer_zone <- create_buffer_zone(site_coords, buffer_distance, crs)
     buffer_data <- extract_buffer_data(raster_data, buffer_zone)
   }
   ```

2. **Candidate Selection:**
   ```r
   n_candidates <- min(10, nrow(buffer_data))
   candidates <- buffer_data[sample(nrow(buffer_data), n_candidates), ]
   ```

3. **Alternative Selection:**
   - **"best"**: Closest to buffer center
   - **"nearest"**: Nearest to original site
   - **"random"**: Random selection from candidates

### Supporting Functions

#### `create_buffer_zone()`

```r
create_buffer_zone(point_coords, buffer_distance, crs = NULL)
```

**Implementation:**
- Creates `sf` point geometry
- Applies `st_buffer()` with specified distance
- Maintains coordinate reference system

#### `extract_buffer_data()`

```r
extract_buffer_data(raster_data, buffer_zone)
```

**Implementation:**
- Uses `terra::extract()` for raster value extraction
- Includes xy coordinates in output
- Removes NA values and ID columns
- Returns clean data.frame with coordinates + environmental values

## Utility Functions

### Application Launcher: `launch_optilhs_app()`

**Purpose:** Launches the Shiny web application interface.

```r
launch_clhs_app(..., port = NULL, host = "127.0.0.1", launch.browser = TRUE)
```

**Features:**
- Flexible app directory detection
- Development and production mode support
- Configurable port and host settings
- Error handling for missing dependencies

## Mathematical Background

### KL Divergence Theory

The Kullback-Leibler divergence measures how one probability distribution differs from another:

**Formula:**
```
KL(P||Q) = Σ P(i) * log(P(i)/Q(i))
```

Where:
- P(i) = sample distribution
- Q(i) = population distribution

**In our context:**
- Measures how well sample represents population
- Lower values indicate better representation
- Used to evaluate different sample sizes

### Exponential Decay Model

**Model:** `KL = b1 * exp(-k * n) + b0`

Where:
- `n` = sample size
- `b0` = asymptotic KL value
- `b1` = initial KL range
- `k` = decay rate

**Optimization Goal:** Find minimum `n` where improvement plateaus

### Simulated Annealing

**Acceptance Probability:**
```
P(accept) = exp(-ΔE / T)
```

Where:
- `ΔE` = change in energy (MSE difference)
- `T` = temperature (decreases over time)
- Higher `T` = more exploration
- Lower `T` = more exploitation

## Implementation Examples

### Basic Usage

```r
# Load package
library(OptiLHS)

# Launch web application
launch_optilhs_app()

# Or use functions directly:

# 1. Load environmental data
raster_data <- terra::rast(c("elevation.tif", "slope.tif", "ndvi.tif"))
pop_data <- as.data.frame(raster_data, xy = TRUE, na.rm = TRUE)

# 2. Optimize sample size
kl_result <- optimize_sample_size(
  population_data = pop_data,
  min_samples = 10,
  max_samples = 200,
  step_size = 10,
  n_replicates = 5
)

optimal_size <- kl_result$optimal_sample_size
cat("Optimal sample size:", optimal_size)

# 3. Optimize sampling design
rf_result <- rf_optimize_sampling(
  covariates = raster_data,
  n_samples = optimal_size,
  n_iterations = 300
)

cat("MSE improvement:", rf_result$improvement, "%")

# 4. Handle inaccessible sites
inaccessible <- data.frame(x = c(100, 200), y = c(150, 250))

alt_result <- process_inaccessible_sites(
  inaccessible_sites = inaccessible,
  raster_data = raster_data,
  buffer_distance = 500
)

cat("Success rate:", alt_result$success_rate * 100, "%")
```

### Advanced Configuration

```r
# Custom KL divergence parameters
kl_custom <- optimize_sample_size(
  population_data = pop_data,
  min_samples = 20,
  max_samples = 1000,
  step_size = 20,
  n_replicates = 15,
  n_bins = 30,
  probability_threshold = 0.98
)

# High-intensity RF optimization
rf_intensive <- rf_optimize_sampling(
  covariates = raster_data,
  n_samples = 150,
  n_iterations = 2000,
  seed = 12345
)

# Strict alternative site selection
alt_strict <- process_inaccessible_sites(
  inaccessible_sites = problem_sites,
  raster_data = raster_data,
  buffer_distance = 100,
  similarity_threshold = 0.95,
  method = "mahalanobis",
  select_method = "best"
)
```

## Performance Considerations

### Computational Complexity

| Function | Time Complexity | Memory Usage | Bottlenecks |
|----------|----------------|--------------|-------------|
| `optimize_sample_size()` | O(n × r × v) | Moderate | cLHS iterations |
| `rf_optimize_sampling()` | O(i × f × t) | High | Random Forest training |
| `process_inaccessible_sites()` | O(s × b) | Low | Buffer extractions |

Where:
- n = number of sample sizes tested
- r = number of replicates
- v = number of variables
- i = optimization iterations
- f = CV folds
- t = RF trees
- s = number of inaccessible sites
- b = buffer size

### Optimization Tips

1. **Sample Size Optimization:**
   - Use smaller `step_size` for fine-tuning
   - Reduce `n_replicates` for faster testing
   - Limit `max_samples` based on study requirements

2. **Random Forest Optimization:**
   - Start with 500 iterations, increase if needed
   - Use seeds for reproducible results
   - Monitor MSE improvement convergence

3. **Alternative Sites:**
   - Adjust `buffer_distance` based on accessibility constraints
   - Use "best" selection for optimal locations
   - Preprocess inaccessible sites to remove duplicates

### Memory Management

```r
# For large datasets:
# 1. Process in chunks
chunk_size <- 10000
n_chunks <- ceiling(nrow(pop_data) / chunk_size)

# 2. Clear intermediate results
rm(intermediate_results)
gc()

# 3. Use appropriate data types
pop_data <- pop_data[, sapply(pop_data, is.numeric)]
```

## References

- Malone, B.P., Minansy, B., Brungard, C. (2019). Some methods to improve the utility of conditioned Latin hypercube sampling. *PeerJ*, 7, e6451. https://doi.org/10.7717/peerj.6451
- Minasny, B., McBratney, A.B. (2006). A conditioned Latin hypercube method for sampling in the presence of ancillary information. *Computers & Geosciences*, 32(9), 1378–1388. https://doi.org/10.1016/j.cageo.2005.12.009
- Wadoux, A. M. J. C., Brus, D. J., & Heuvelink, G. B. M. (2019). Sampling design optimization for soil mapping with random forest. *Geoderma*, 355, 113913. https://doi.org/10.1016/j.geoderma.2019.113913.

## Technical Support

For technical issues or questions:

1. Check function documentation: `help("function_name")`
2. Review this technical guide
3. Consult the `CLAUDE.md` file for development guidance
4. Ensure all dependencies are properly installed

---

*This technical guide documents OptiLHS v0.1.1. For the latest updates and features, refer to the package documentation and DESCRIPTION file.*