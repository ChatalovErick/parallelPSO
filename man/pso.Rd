\name{pso}
\alias{pso}

\title{pso}
\description{
  Minimization of a fitness function using PSO. The algorithm can be run sequentially or in parallel.}

\usage{
pso(fitness_function, number_parameters, number_of_partiples, max_number_iterations, W.1 = 0.9, W.2 = 0.4, C.1i = 0.5 + log(2), C.1f = 0.5 + log(2), C.2i = 0.5 + log(2), C.2f = 0.5 + log(2), K = 3, parameters_bounds = cbind(rep(-1, number_parameters), rep(1, number_parameters)), Vmax, trace = TRUE, parallel = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitness_function}{
    A function to be minimized (or maximized), with first
    argument the vector of parameters over which minimization is
    to take place.  It should return a scalar result.
}
  \item{number_parameters}{
    Number of parameters that are used to use the fitness fucntion.
}
  \item{number_of_partiples}{
    The swarm size.
}
  \item{max_number_iterations}{
    The maximum number of iterations.
}

\item{W}{
  The exploitation constant. Inertia

  \item{W.1}{
    Initial value of the inertia weight Defaults to 0.9.
  }
  \item{W.2}{
    Final value of the inertia weight Defaults to 0.4.
  }

}
\item{C.1}{
  The local exploration constant. Defaults to \code{.5+log(2)}.
  \item{C.1i}{
    Initial value of the C.1.
  }
    \item{C.1f}{
    Final value of the C.1.
  }
}
\item{C.2}{
  The global exploration constant. Defaults to \code{.5+log(2)}.
  \item{C.2i}{
    Initial value of the C.2.
  }
  \item{C.2f}{
    Final value of the C.2.
  }
}

  \item{K}{
  The exponent for calculating number of informants. Defaults to \code{3}
}
  \item{parameters_bounds}{
%%     ~~Describe \code{parameters_bounds} here~~
}

  \item{trace}{
  Return the statistics collecte during the optimization.
}
  \item{parallel}{
  An optional argument which allows to specify if PSO should be run sequentially or in parallel
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
M. S. Alam and M. O. Tokhi. Dynamic modelling of a single-link flexible
manipulator system: a particle swarm optimisation approach.
https://doi.org/10.1260/026309207781487466
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (fitness_function, number_parameters, number_of_partiples,
    max_number_iterations, W.1 = 0.9, W.2 = 0.4, C.1i = 0.5 +
        log(2), C.1f = 0.5 + log(2), C.2i = 0.5 + log(2), C.2f = 0.5 +
        log(2), K = 3, parameters_bounds = cbind(rep(-1, number_parameters),
        rep(1, number_parameters)), Vmax, trace = TRUE, parallel = FALSE,
    ...)
{
    callargs <- list(...)
    mrunif <- function(n, m, lower, upper) {
        return(matrix(runif(n * m, 0, 1), nrow = n, ncol = m) *
            (upper - lower) + lower)
    }
    lower = as.double(parameters_bounds[, 1], number_parameters)
    upper = as.double(parameters_bounds[, 2], number_parameters)
    p.p <- 1 - (1 - 1/number_parameters)^K
    if (is.logical(parallel)) {
        if (parallel) {
            parallel <- startParallel(parallel)
            stopCluster <- TRUE
        }
        else {
            parallel <- stopCluster <- FALSE
        }
    }
    else {
        stopCluster <- if (inherits(parallel, "cluster"))
            FALSE
        else TRUE
        parallel <- startParallel(parallel)
    }
    on.exit(if (parallel & stopCluster) stopParallel(attr(parallel,
        "cluster")))
    `\%DO\%` <- if (parallel && requireNamespace("doRNG", quietly = TRUE))
        doRNG::`\%dorng\%`
    else if (parallel)
        `\%dopar\%`
    else `\%do\%`
    X <- mrunif(number_parameters, number_of_partiples, lower,
        upper)
    V <- (mrunif(number_parameters, number_of_partiples, lower,
        upper) - X)/2
    if (isTRUE(trace)) {
        trace.it <- c()
        trace.error <- c()
        trace.f <- NULL
        trace.x <- NULL
    }
    if (!(parallel)) {
        f.x <- apply(X, 2, function(x) {
            do.call(fitness_function, c(list(x), callargs))
        })
    }
    else {
        print("going parallel")
        f.x <- rep(NA, number_of_partiples)
        f.x <- foreach(i = seq_len(number_of_partiples), .combine = "c") \%DO\%
            {
                do.call(fitness_function, c(list(X[, i]), callargs))
            }
    }
    P <- X
    f.p <- f.x
    error <- f.p[which.min(f.p)]
    P.best <- P[, which.min(f.p)]
    if (isTRUE(trace)) {
        message("It 1: fitness=", signif(error, 4))
        trace.it <- c(trace.it, 1)
        trace.error <- c(trace.error, error)
        trace.f <- c(trace.f, list(f.x))
        trace.x <- c(trace.x, list(X))
    }
    stats.iter <- 1
    iter <- 1
    while (stats.iter < max_number_iterations) {
        if (W.1 == W.2) {
            W = W.1
        }
        else {
            W = W.2 + (W.1 - W.2) * (((max_number_iterations -
                1) - stats.iter)/(max_number_iterations - 1))
        }
        if (C.1i == C.1f) {
            C.1 = C.1i
        }
        else {
            C.1 = C.1i + (C.1f - C.1i) * (stats.iter/(max_number_iterations -
                1))
        }
        if (C.2i == C.2f) {
            C.2 = C.2i
        }
        else {
            C.2 = C.2i + (C.2f - C.2i) * (stats.iter/(max_number_iterations -
                1))
        }
        stats.iter <- stats.iter + 1
        if (K != 0) {
            neighborhood <- matrix(runif(number_of_partiples *
                number_of_partiples, 0, 1) <= p.p, number_of_partiples,
                number_of_partiples)
            diag(neighborhood) <- TRUE
            for (i in 1:number_of_partiples) {
                j <- which(neighborhood[, i])[which.min(f.p[neighborhood[,
                  i]])]
                V[, i] <- W * V[, i] + runif(number_parameters,
                  0, C.2) * (P[, i] - X[, i])
                if (i != j) {
                  V[, i] <- V[, i] + runif(number_parameters,
                    0, C.1) * (P[, j] - X[, i])
                }
            }
        }
        else {
            V <- W * V + runif(number_parameters, 0, C.2) * (-sweep(X,
                1, P.best, "-")) + runif(number_parameters, 0,
                C.1) * (P - X)
        }
        X <- X + V
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
        if (!(parallel)) {
            for (i in 1:number_of_partiples) {
                f.x[i] <- do.call(fitness_function, c(list(X[,
                  i]), callargs))
            }
        }
        else {
            f.x <- rep(NA, number_of_partiples)
            f.x <- foreach(i = seq_len(number_of_partiples),
                .combine = "c") \%DO\% {
                do.call(fitness_function, c(list(X[, i]), callargs))
            }
        }
        for (i in 1:number_of_partiples) {
            if (f.x[i] < f.p[i]) {
                P[, i] <- X[, i]
                f.p[i] <- f.x[i]
            }
        }
        if (f.p[which.min(f.p)] < error) {
            P.best <- P[, which.min(f.p)]
            error <- f.p[which.min(f.p)]
        }
        if (isTRUE(trace)) {
            message("It ", stats.iter, ": fitness=", signif(error,
                4))
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
    if (isTRUE(trace))
        message(msg)
    res <- list(par = P.best, value = error, counts = c(iteration = stats.iter),
        convergence = msgcode, message = msg)
    if (isTRUE(trace)) {
        res <- c(res, list(stats = list(it = trace.it, error = trace.error,
            f = trace.f, x = trace.x)))
    }
    if (isTRUE(parallel)) {
        stopParallel(parallel)
    }
    invisible(gc())
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
