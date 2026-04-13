# parallelPSO

An R package for Particle Swarm Optimization (PSO) with support for parallel processing using doParallel. This package is designed for computationally intensive fitness functions and long simulations.

## Features

- **Particle Swarm Optimization**: Implementation of the PSO algorithm with customizable parameters
- **Parallel Processing**: Run optimization in parallel using multiple CPU cores via doParallel
- **Flexible Fitness Functions**: Easily define your own objective function depending on the problem at hand
- **Neighborhood Topology**: Support for random informant-based neighborhood structure
- **Adaptive Parameters**: Time-varying inertia weight and acceleration coefficients
- **Comprehensive Tracing**: Optional detailed statistics collection during optimization

## Installation

### From GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")
library(devtools)

# Install parallelPSO from GitHub
remotes::install_github("ChatalovErick/parallelPSO")
```

### Manual Installation:

```r
# Clone or download the repository
# Then install from local directory
install.packages(".", repos = NULL, type = "source")
```

## Dependencies

The package requires the following packages:
- `foreach` - Framework for looping constructs
- `doParallel` - Parallel execution backend for foreach
- `parallel` - Base R package for parallel computing

Optional:
- `doRNG` - For reproducible parallel computations

## Usage

### Basic Example (Sequential):

```r
library(parallelPSO)

# Define a simple fitness function (Sphere function)
sphere <- function(x) sum(x^2)

# Run PSO sequentially
result <- pso(
  fitness_function = sphere,
  number_parameters = 5,
  number_of_partiples = 20,
  max_number_iterations = 100,
  parallel = FALSE
)

# View results
print(result$par)    # Best parameters found
print(result$value)  # Best fitness value
```

### Parallel Example:

```r
# Run PSO in parallel (automatically uses available cores)
result_parallel <- pso(
  fitness_function = sphere,
  number_parameters = 5,
  number_of_partiples = 20,
  max_number_iterations = 100,
  parallel = TRUE
)
```

### Using the Ackley Function:

```r
# Optimize the Ackley function (challenging test function)
result_ackley <- pso(
  fitness_function = ackley,
  number_parameters = 10,
  number_of_partiples = 30,
  max_number_iterations = 200,
  parameters_bounds = cbind(rep(-32, 10), rep(32, 10)),
  parallel = TRUE
)
```

### Custom Fitness Function with Additional Parameters:

```r
# Define a fitness function with additional parameters
custom_fitness <- function(x, alpha = 1, beta = 2) {
  sum((x - alpha)^2) + beta * sum(x^4)
}

# Run PSO with custom parameters
result_custom <- pso(
  fitness_function = custom_fitness,
  number_parameters = 5,
  number_of_partiples = 25,
  max_number_iterations = 150,
  alpha = 0.5,
  beta = 0.1,
  parallel = FALSE
)
```

## Algorithm Parameters

The PSO algorithm can be customized with the following parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `W.1` | Initial inertia weight | 0.9 |
| `W.2` | Final inertia weight | 0.4 |
| `C.1i`, `C.1f` | Cognitive coefficient (initial, final) | 0.5 + log(2) |
| `C.2i`, `C.2f` | Social coefficient (initial, final) | 0.5 + log(2) |
| `K` | Neighborhood size exponent | 3 |
| `parameters_bounds` | Lower and upper bounds for parameters | [-1, 1] |
| `trace` | Return optimization statistics | TRUE |
| `parallel` | Enable parallel processing | FALSE |

## Output

The `pso()` function returns a list containing:

- `par`: Best parameter values found
- `value`: Fitness value at best parameters
- `counts`: Number of iterations performed
- `convergence`: Convergence code
- `message`: Termination message
- `stats`: (if trace=TRUE) Iteration history including:
  - `it`: Iteration numbers
  - `error`: Best fitness at each iteration
  - `f`: Fitness values of all particles
  - `x`: Particle positions

## References

1. Alam, M. S., & Tokhi, M. O. (2007). Dynamic modelling of a single-link flexible manipulator system: a particle swarm optimisation approach. https://doi.org/10.1260/026309207781487466

2. Kennedy, J., & Eberhart, R. (1995). Particle swarm optimization. Proceedings of ICNN'95 - International Conference on Neural Networks.

3. Surjanovic, S., & Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. Simon Fraser University. http://www.sfu.ca/~ssurjano/ackley.html

## License

This package is licensed under GPL-2.

## Author

Erick Chatalov <chataloverick@gmail.com>

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.