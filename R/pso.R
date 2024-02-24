pso <- function(fitness_function,number_parameters,number_of_partiples,max_number_iterations,
                W=1/(2*log(2)),C.1 = .5+log(2),C.2 = .5+log(2),parameters_bounds = cbind(rep(-1,number_parameters),rep(1,number_parameters)),
                Vmax,trace=TRUE,parallel=FALSE, ...) {
  
  
  callargs <- list(...)

  mrunif <- function(n,m,lower,upper) {
    return(matrix(runif(n*m,0,1),nrow=n,ncol=m)*(upper-lower)+lower)
  }
  
  # parameters intervals.
  lower = as.double(parameters_bounds[,1],number_parameters)
  upper = as.double(parameters_bounds[,2],number_parameters)
  
  # Start parallel computing (if needed)
  if(is.logical(parallel))
    { if(parallel) 
      { parallel <- startParallel(parallel)
        stopCluster <- TRUE }
      else
      { parallel <- stopCluster <- FALSE } 
  }
  else
    { stopCluster <- if(inherits(parallel, "cluster")) FALSE else TRUE
      parallel <- startParallel(parallel) 
    }
  on.exit(if(parallel & stopCluster)
    stopParallel(attr(parallel, "cluster")) )
  # define operator to use depending on parallel being TRUE or FALSE
  `%DO%` <- if(parallel && requireNamespace("doRNG", quietly = TRUE)) 
                            doRNG::`%dorng%`
            else if(parallel) `%dopar%` else `%do%`
  
  
  # Initialization
  ### initial particles
  X <- mrunif(number_parameters,number_of_partiples,lower,upper)

  ### initial velocities
  V <- (mrunif(number_parameters,number_of_partiples,lower,upper)-X)/2
  
  if (isTRUE(trace)) {
      trace.it <- c()
      trace.error <- c()
      trace.f <- NULL
      trace.x <- NULL
  }

  # first evaluations and initial population
  if(!(parallel)){
    f.x <- apply(X,2,function(x){ do.call(fitness_function,c(list(x),callargs)) }) # first evaluations  
  } else {
    print("going parallel")
    f.x <- rep(NA, number_of_partiples)
    f.x <- foreach(i = seq_len(number_of_partiples), .combine= "c") %DO%
      {
        do.call(fitness_function,c(list(X[,i]),callargs)) 
      }
  }

  P <- X
  f.p <- f.x
  error <- f.p[which.min(f.p)] # lower error found 
  P.best  <- P[,which.min(f.p)] # parameterization with the lower fitness
  
  if (isTRUE(trace)) {
    message("It 1: fitness=",signif(error,4))
    trace.it <- c(trace.it,1)
    trace.error <- c(trace.error,error)
    trace.f <- c(trace.f,list(f.x))
    trace.x <- c(trace.x,list(X))
  }
  
  # iterations 
  stats.iter <- 1
  while(stats.iter<max_number_iterations) {
    
    stats.iter <- stats.iter+1 

    # velocity
    V <- W*V + runif(number_parameters,0,C.2)*(-sweep(X,1,P.best,"-")) + runif(number_parameters,0,C.1)*(P-X)
      
    # new positions 
    X <- X+V
    
    ## check bounds ##
    ##################
    for(i in 1:number_of_partiples){
      
      temp <- X[,i]<lower
      if (any(temp)) {
        X[temp,i] <- lower[temp]
        V[temp,i] <- 0
      }
      temp <- X[,i]>upper
      if (any(temp)) {
        X[temp,i] <- upper[temp]
        V[temp,i] <- 0
      }
    }
    ##################
    ##################
    
    #### functions evaluations ####
    ###############################
    
    if (!(parallel)){
      # loop for the particles 
      for(i in 1:number_of_partiples){
        # evaluate function
        f.x[i] <- do.call(fitness_function,c(list(X[,i.]),callargs)) 
      }
    } else {
      
      f.x <- rep(NA, number_of_partiples)
      f.x <- foreach(i= seq_len(number_of_partiples),.combine= "c") %DO% {
        do.call(fitness_function,c(list(X[,i]),callargs)) 
      }
      
    }

    ###############################
    ###############################
    
    ## improvements ##
    ##################
    
    # best parameterizations for iteration i 
    for (i in 1:number_of_partiples){
      if(f.x[i]< f.p[i]){
        P[,i] <- X[,i]
        f.p[i] <- f.x[i]
      }
    } 
    
    
    if (f.p[which.min(f.p)] < error){
      P.best <- P[,which.min(f.p)] # best parameterization from all optimization
      error <- f.p[which.min(f.p)] # best fitness value found 
    
    }
    
    ##################
    if (isTRUE(trace)) {
        message("It ",stats.iter,": fitness=",signif(error,4))
        trace.it <- c(trace.it,stats.iter)
        trace.error <- c(trace.error,error)
        trace.f <- c(trace.f,list(f.p))
        trace.x <- c(trace.x,list(P))
    }

  }
  
  if (stats.iter>=max_number_iterations) {
    msg <- "Maximal number of iterations reached"
    msgcode <- 2
  }
  
  if (isTRUE(trace)) message(msg)
  
  res <- list(par=P.best,value=error,
              counts=c("iteration"=stats.iter),
              convergence=msgcode,message=msg)
  
  if (isTRUE(trace)) {
    res <- c(res,list(stats=list(it=trace.it,
                               error=trace.error,
                               f=trace.f,
                               x=trace.x)))
  }
  
  if (isTRUE(parallel)){
    stopParallel(parallel)
  }
  
  invisible(gc())
  return(res)
}
