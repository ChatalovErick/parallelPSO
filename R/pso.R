#' Particle Swarm Optimization (PSO) Algorithm
#'
#' Minimization of a fitness function using Particle Swarm Optimization.
#' The algorithm can be run sequentially or in parallel using doParallel,
#' making it suitable for computationally intensive fitness functions and long simulations.
#'
#' @param fitness_function A function to be minimized (or maximized), with first
#'    argument the vector of parameters over which minimization is to take place.
#'    It should return a scalar result. Additional arguments can be passed via ...
#' @param number_parameters Number of parameters used by the fitness function.
#' @param number_of_partiples The swarm size (number of particles).
#' @param max_number_iterations The maximum number of iterations.
#' @param W.1 Initial value of the inertia weight. Defaults to 0.9.
#' @param W.2 Final value of the inertia weight. Defaults to 0.4.
#' @param C.1i Initial value of the cognitive (local exploration) constant. Defaults to 0.5 + log(2).
#' @param C.1f Final value of the cognitive constant. Defaults to 0.5 + log(2).
#' @param C.2i Initial value of the social (global exploration) constant. Defaults to 0.5 + log(2).
#' @param C.2f Final value of the social constant. Defaults to 0.5 + log(2).
#' @param K The exponent for calculating the number of informants (neighborhood size). Defaults to 3.
#' @param parameters_bounds Matrix with two columns specifying lower and upper bounds for parameters.
#'    Defaults to [-1, 1] for all parameters.
#' @param Vmax Maximum velocity (currently not implemented).
#' @param trace Logical; if TRUE, return statistics collected during optimization.
#' @param parallel Logical; if TRUE, run PSO in parallel using doParallel. If FALSE, run sequentially.
#'    Alternatively, a cluster object can be provided.
#' @param ... Additional arguments passed to the fitness function.
#'
#' @return A list containing:
#' \describe{
#'   \item{par}{The best parameter values found.}
#'   \item{value}{The fitness value at the best parameters.}
#'   \item{counts}{A named vector with the number of iterations performed.}
#'   \item{convergence}{Convergence code (2 = maximal number of iterations reached).}
#'   \item{message}{A message describing the termination reason.}
#'   \item{stats}{If trace=TRUE, a list containing iteration history:
#'     \describe{
#'       \item{it}{Iteration numbers.}
#'       \item{error}{Best fitness value at each iteration.}
#'       \item{f}{Fitness values of all particles at each iteration.}
#'       \item{x}{Particle positions at each iteration.}
#'     }}
#' }
#'
#' @details
#' The PSO algorithm uses a neighborhood topology based on random informants.
#' Each particle is influenced by its own best position and the best position
#' among its informants. The inertia weight and acceleration coefficients can
#' be varied linearly over time to balance exploration and exploitation.
#'
#' When running in parallel mode, fitness function evaluations are distributed
#' across multiple cores using foreach and doParallel, which can significantly
#' speed up optimization for computationally expensive fitness functions.
#'
#' @references
#' M. S. Alam and M. O. Tokhi. Dynamic modelling of a single-link flexible
#' manipulator system: a particle swarm optimisation approach.
#' https://doi.org/10.1260/026309207781487466
#'
#' Kennedy, J., & Eberhart, R. (1995). Particle swarm optimization.
#' Proceedings of ICNN'95 - International Conference on Neural Networks.
#'
#' @examples
#' \dontrun{
#' # Define a simple fitness function (Sphere function)
#' sphere <- function(x) sum(x^2)
#'
#' # Run PSO sequentially
#' result_seq <- pso(fitness_function = sphere,
#'                   number_parameters = 5,
#'                   number_of_partiples = 20,
#'                   max_number_iterations = 100,
#'                   parallel = FALSE)
#'
#' # Run PSO in parallel (requires multiple cores)
#' result_par <- pso(fitness_function = sphere,
#'                   number_parameters = 5,
#'                   number_of_partiples = 20,
#'                   max_number_iterations = 100,
#'                   parallel = TRUE)
#'
#' # Print results
#' print(result_seq$par)
#' print(result_seq$value)
#' }
#'
#' @importFrom foreach %dopar% foreach
#' @importFrom doParallel registerDoParallel getDoParWorkers
#' @importFrom parallel stopCluster makeCluster detectCores clusterExport
#' @export
pso <- function(fitness_function, number_parameters, number_of_partiples, max_number_iterations,
                W.1 = 0.9, W.2 = 0.4, C.1i = .5+log(2), C.1f = .5+log(2),
                C.2i = .5+log(2), C.2f = .5+log(2), K = 3,
                parameters_bounds = cbind(rep(-1, number_parameters), rep(1, number_parameters)),
                Vmax = NULL, trace = TRUE, parallel = FALSE, ...) {
  
  callargs <- list(...)

  mrunif <- function(n, m, lower, upper) {
    return(matrix(runif(n * m, 0, 1), nrow = n, ncol = m) * (upper - lower) + lower)
  }
  
  # parameters intervals. 
  lower <- as.double(parameters_bounds[, 1])
  upper <- as.double(parameters_bounds[, 2])
  
  # % of informants. based in https://github.com/cran/pso.git. neighborhood method !?
  p.p <- 1 - (1 - 1 / number_parameters)^K

  # Start parallel computing (if needed)
  cluster_created <- FALSE
  if (isTRUE(parallel)) {
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' is required for parallel execution. Please install it.")
    }
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Package 'foreach' is required for parallel execution. Please install it.")
    }
    
    # Create cluster with number of available cores minus 1
    n_cores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    cluster_created <- TRUE
    
    # Export necessary variables and functions to cluster
    parallel::clusterExport(cl, c("fitness_function", "callargs"), envir = environment())
  }
  
  # Define the evaluation operator based on parallel mode
  if (isTRUE(parallel) && requireNamespace("doRNG", quietly = TRUE)) {
    `%DO%` <- doRNG::`%dorng%`
  } else if (isTRUE(parallel)) {
    `%DO%` <- foreach::`%dopar%`
  } else {
    `%DO%` <- foreach::`%do%`
  }
  
  # Initialization
  ### initial particles
  X <- mrunif(number_parameters, number_of_partiples, lower, upper)

  ### initial velocities
  V <- (mrunif(number_parameters, number_of_partiples, lower, upper) - X) / 2
  
  if (isTRUE(trace)) {
    trace.it <- c()
    trace.error <- c()
    trace.f <- NULL
    trace.x <- NULL
  }

  # first evaluations and initial population
  if (!isTRUE(parallel)) {
    f.x <- apply(X, 2, function(x) {
      do.call(fitness_function, c(list(x), callargs))
    }) # first evaluations  
  } else {
    message("Running in parallel mode with ", parallel::getDoParWorkers(), " workers")
    f.x <- foreach::foreach(i = seq_len(number_of_partiples), .combine = c) %DO% {
      do.call(fitness_function, c(list(X[, i]), callargs))
    }
  }

  P <- X
  f.p <- f.x
  error <- f.p[which.min(f.p)] # lower error found 
  P.best <- P[, which.min(f.p)] # parameterization with the lower fitness
  
  if (isTRUE(trace)) {
    message("It 1: fitness=", signif(error, 4), " | Mean fitness=", signif(mean(f.p), 4))
    trace.it <- c(trace.it, 1)
    trace.error <- c(trace.error, error)
    trace.f <- c(trace.f, list(f.x))
    trace.x <- c(trace.x, list(X))
  }
  
  # iterations 
  stats.iter <- 1
  iter <- 1
  while (stats.iter < max_number_iterations) {
    
    # Calculate the W, inertia. 
    if (W.1 == W.2) {
      # fixed inertia 
      W <- W.1
    } else {
      # acceleration over time for the inertia 
      W <- W.2 + (W.1 - W.2) * (((max_number_iterations - 1) - stats.iter) / (max_number_iterations - 1))
    }

    # Calculate the C.1
    if (C.1i == C.1f) {
      # fixed C.1
      C.1 <- C.1i
    } else {
      # acceleration over time for the C.1 
      C.1 <- C.1i + (C.1f - C.1i) * (stats.iter / (max_number_iterations - 1))
    }

    # Calculate the C.2
    if (C.2i == C.2f) {
      # fixed C.2
      C.2 <- C.2i
    } else {
      # acceleration over time for the C.2 
      C.2 <- C.2i + (C.2f - C.2i) * (stats.iter / (max_number_iterations - 1))
    }
    # iterations for the while loop.
    stats.iter <- stats.iter + 1 
    
    # Calculate the velocity
    # make a neighborhood for the pso based in https://github.com/cran/pso.git 

    if (K != 0) {
      neighborhood <- matrix(runif(number_of_partiples * number_of_partiples, 0, 1) <= p.p, 
                            number_of_partiples, number_of_partiples)
      diag(neighborhood) <- TRUE

      for (i in 1:number_of_partiples) {
        
        j <- which(neighborhood[, i])[which.min(f.p[neighborhood[, i]])] # best informant

        # velocity
        V[, i] <- W * V[, i] + runif(number_parameters, 0, C.2) * (P[, i] - X[, i])
        if (i != j) {
          V[, i] <- V[, i] + runif(number_parameters, 0, C.1) * (P[, j] - X[, i])
        }
      }

    } else {
      # velocity
      V <- W * V + runif(number_parameters, 0, C.2) * (-sweep(X, 1, P.best, "-")) + 
           runif(number_parameters, 0, C.1) * (P - X)
    }

    # new positions 
    X <- X + V
    
    ## check bounds ##
    ##################
    for (i in 1:number_of_partiples) {
      temp <- X[, i] < lower
      if (any(temp)) {
        X[temp, i] <- lower[temp]
        V[temp, i] <- 0
      }
      temp <- X[, i] > upper
      if (any(temp)) {
        X[temp, i] <- upper[temp]
        V[temp, i] <- 0
      }
    }
    ##################
    ##################
    
    #### functions evaluations ####
    ###############################
    
    if (!isTRUE(parallel)) {
      # loop for the particles 
      for (i in 1:number_of_partiples) {
        # evaluate function
        f.x[i] <- do.call(fitness_function, c(list(X[, i]), callargs))
      }
    } else {
      
      f.x <- foreach::foreach(i = seq_len(number_of_partiples), .combine = c) %DO% {
        do.call(fitness_function, c(list(X[, i]), callargs))
      }
      
    }

    ###############################
    ###############################
    
    ## improvements ##
    ##################
    
    # best parameterizations for iteration i 
    for (i in 1:number_of_partiples) {
      if (f.x[i] < f.p[i]) {
        P[, i] <- X[, i]
        f.p[i] <- f.x[i]
      }
    } 
    
    
    if (f.p[which.min(f.p)] < error) {
      P.best <- P[, which.min(f.p)] # best parameterization from all optimization
      error <- f.p[which.min(f.p)] # best fitness value found
    
    }
    
    ##################
    if (isTRUE(trace)) {
      message("It ", stats.iter, ": fitness=", signif(error, 4), 
              " | Mean fitness=", signif(mean(f.p), 4))
      trace.it <- c(trace.it, stats.iter)
      trace.error <- c(trace.error, error)
      trace.f <- c(trace.f, list(f.p))
      trace.x <- c(trace.x, list(P))
    }

  }
  
  if (stats.iter >= max_number_iterations) {
    msg <- "Maximal number of iterations reached"
    msgcode <- 2
  }
  
  if (isTRUE(trace)) message(msg)
  
  res <- list(par = P.best, value = error,
              counts = c("iteration" = stats.iter),
              convergence = msgcode, message = msg)
  
  if (isTRUE(trace)) {
    res <- c(res, list(stats = list(it = trace.it,
                                   error = trace.error,
                                   f = trace.f,
                                   x = trace.x)))
  }
  
  # Stop parallel cluster if created
  if (cluster_created) {
    parallel::stopCluster(cl)
    parallel::stopImplicitCluster()
  }
  
  invisible(gc())
  return(res)
}
