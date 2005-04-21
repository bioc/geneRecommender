.packageName <- "geneRecommender"

"gr.main" <-
function(normalized.dataset, query, fun = median, ngenes = NULL, extra = FALSE){

  #check for errors in input command
  proceed <- FALSE
  if (is(normalized.dataset, "matrix") == TRUE){
    data.format <- "matrix"
    proceed <- TRUE
    if (is.null(dimnames(normalized.dataset)[[1]])){
      warning("Input dataset did not have row names.  ")
      proceed <- FALSE
    }
    if (mode(normalized.dataset) != "numeric"){
      warning("Input dataset was not numeric.  ")
      proceed <- FALSE
    } 
  }
  if (is(normalized.dataset, "exprSet") == TRUE){
    data.format <- "exprSet"
    proceed <- TRUE
    if (is.null(dimnames(exprs(normalized.dataset))[[1]])){
      warning("Input dataset did not have row names.  ")
      proceed <- FALSE
    }
    if (mode(exprs(normalized.dataset)) != "numeric"){
      warning("Input dataset was not numeric.  ")
      proceed <- FALSE
    } 
  }
  if (proceed == FALSE){
    warning("Input dataset was not of class matrix or class exprSet.  ")
  }

  if (is.vector(query) != TRUE){
    warning("Input query was not a vector.  ")
    proceed <- FALSE
  }
  if (length(query) < 2){
    warning("Input query contained less than 2 elements.  ")
    proceed <- FALSE
  }
  if (is.function(fun) != TRUE){
    warning("Argument fun was not a function.  ")
    proceed <- FALSE
  } 
  if (is.null(ngenes) != TRUE){
    if (is.finite(ngenes) != TRUE){
      warning("Argument ngenes was not a real number.  ")
      proceed <- FALSE
    }
  } 
  if (is.logical(extra) != TRUE){
    warning("Argument extra was not logical.  ")
    proceed <- FALSE
  } 
  if (proceed != TRUE){
    warning("Function aborted due to invalid arguments.  ")
  } 

  #proceed if no errors were found
  if (proceed == TRUE){

    #define y.i.j
    if (data.format == "matrix"){
      y.i.j <- normalized.dataset
      rm(normalized.dataset)
    }
    if (data.format == "exprSet"){
      y.i.j <- exprs(normalized.dataset)
      rm(normalized.dataset)
    }


    #initialize variables
    n <- dim(y.i.j)[1]
    p <- dim(y.i.j)[2] 
    r.i.j.index.i <- array(NA, dim = c(p))
    k.j <- array(NA, dim = c(p))
    y.b.q.j <- array(NA, dim = c(p))
    v.h.q.j <- array(NA, dim = c(p))
    z.e.j <- array(NA, dim = c(p))
    s.g.i <- array(NA, dim = c(n))
    pre.s.g.i <- array(0, dim = c(n))
    total.observations <- array(0, dim = c(n))

    #define q and y.q.i.j
    q <- match(query, dimnames(y.i.j)[[1]])
    y.q.i.j <- y.i.j[q,]

    #define z.e.j
    for (index.j in 1:p){
      k.j[index.j] <- length(na.omit(y.q.i.j[,index.j]))
      y.b.q.j[index.j] <- median(y.q.i.j[,index.j], na.rm = TRUE)
      #There are two different definitions of the term "sample variance."
      #This program uses the same definition used in Art Owen's code.  
      v.h.q.j[index.j] <-  ((k.j[index.j] - 1)/k.j[index.j]) * var(y.q.i.j[,index.j], na.rm = TRUE)
    }

    min.k.j <- min(5,length(q))
    k.j <- ifelse(k.j >= min.k.j, k.j, NA)
    z.e.j <- (k.j^(1/2))*y.b.q.j/(v.h.q.j + 1/(3*p^2))^(1/2)
 
    #rank z.e.j
    z.e.j <- abs(z.e.j)
    z.e.j <- ifelse(is.finite(z.e.j), z.e.j, -Inf)
    ranked.j <- sort(z.e.j, index.return = TRUE, decreasing = TRUE, na.last = NA)[[2]]

    #find the optimal number of experiments.  
    #The number of experiments which minimizes score is optimal
    score.best <- Inf
    for (index in 1:(length(ranked.j))){
      index.j <- ranked.j[index]
      pre.s.g.i <- pre.s.g.i + ifelse(is.finite(y.i.j[,index.j] * y.b.q.j[index.j]), y.i.j[,index.j] * y.b.q.j[index.j], 0)
      total.observations <- total.observations + ifelse(is.finite(y.i.j[,index.j] * y.b.q.j[index.j]), 1, 0)
      s.g.i <- pre.s.g.i / total.observations
      s.g.i <- ifelse(is.finite(s.g.i), s.g.i, -Inf)
      ranked.s.g.i <- sort(s.g.i, index.return = TRUE, decreasing = TRUE, na.last = NA)[[2]]
      score <- fun(match(q, ranked.s.g.i))
      if (score <= score.best){
        score.best <- score
        s.g.i.best <- s.g.i
        index.best <- index
      }
    }

    #calculate number of genes found at 50% recall
    sorted.s.g.i <- sort(s.g.i.best, index.return = TRUE, decreasing = TRUE, na.last = NA)
    counter <- 0
    index <- 0
    while (counter < as.integer((1 + length(q))/2)){
      index <- index + 1
      if (is.finite(match(sorted.s.g.i[[2]][index], q))){
        counter <- counter + 1
      }
    }
    fifty.percent.recall <- index
    if (is.null(ngenes) == TRUE){
      ngenes <- fifty.percent.recall
    }

    #only the ngenes highest scoring genes will be listed
    listed.genes <- sorted.s.g.i[[2]][1:ngenes]
    
    #generate result
    main.result <- cbind( dimnames(y.i.j)[[1]][listed.genes],  ifelse(is.na(match(listed.genes, q)), "Not within query", "Within query")  )

    #compute the extra outputs
    if (extra == TRUE){

      #determine included/excluded experiments and compute s.g.i and z.g.i and contribution
      experiments.included <- c()
      experiments.excluded <- c()
      sorted.s.g.i.result <- array(sorted.s.g.i[[1]][1:ngenes], dim = c(ngenes))
      z.g.i.numerator <- array(0, dim = c(ngenes))
      z.g.i.denominator <- array(0, dim = c(ngenes))
      contribution <- array(NA, dim = c(ngenes,index.best))
      index <- 0
      for (index.j in 1:p){

        #check if experiment index.j was included.  If so, then use it in computing z.g.i and contribution      
        if (is.finite(match(index.j, ranked.j[1:index.best]))){
          index <- index + 1
          experiments.included <- c(experiments.included, index.j)
          present <- ifelse(is.finite(y.b.q.j[index.j]*y.i.j[listed.genes,index.j]), 1, 0)
          z.g.i.numerator <- z.g.i.numerator + ifelse(present == 1, y.b.q.j[index.j]*y.i.j[listed.genes, index.j], 0)
          z.g.i.denominator <- z.g.i.denominator + ifelse(present == 1, y.b.q.j[index.j]^2, 0)
          contribution[,index] <- ifelse(present == 1, y.b.q.j[index.j] * y.i.j[listed.genes,index.j], 0)
        }

        #otherwise, just add it to the list of excluded experiments
        else{
          experiments.excluded <- c(experiments.excluded, index.j)
        }
      }

      #if the experiments have names, then write experiment.included and experiment.excluded in terms of these names
      if (is.null(dimnames(y.i.j)[[2]]) != TRUE){
        experiments.included <- dimnames(y.i.j)[[2]][experiments.included]
        if (is.null(experiments.excluded) != TRUE){
          experiments.excluded <- dimnames(y.i.j)[[2]][experiments.excluded]
        }
      }

      #finish computing z.g.i
      z.g.i.denominator <- ((1/3)*z.g.i.denominator)^(1/2)
      sorted.z.g.i.result <- array(z.g.i.numerator/z.g.i.denominator, dim = c(ngenes))

      #label results
      dimnames(sorted.s.g.i.result) <- list(main.result[,1])
      dimnames(sorted.z.g.i.result) <- list(main.result[,1])
      dimnames(contribution) <- list(main.result[,1],experiments.included)

      #put all the results into one big list
      result <- list(main.result = main.result, fifty.percent.recall = fifty.percent.recall, experiments.included = experiments.included, experiments.excluded = experiments.excluded, s.g.i = sorted.s.g.i.result, z.g.i = sorted.z.g.i.result, contribution = contribution)
    }

    #otherwise, just put the result into a one item list
    else{
      result <- list(main.result = main.result)
    }

    #return result
    result
  }

  #in the case of bad input, return NA
  else{
    NA
  }
}

"gr.normalize" <-
function(unnormalized.dataset){

  #check for errors in input command
  proceed <- FALSE
  if (is(unnormalized.dataset, "matrix") == TRUE){
    data.format <- "matrix"
    proceed <- TRUE
    if (is.null(dimnames(unnormalized.dataset)[[1]])){
      warning("Input dataset did not have row names.  ")
      proceed <- FALSE
    }
    if (mode(unnormalized.dataset) != "numeric"){
      warning("Input dataset was not numeric.  ")
      proceed <- FALSE
    } 
  }
  if (is(unnormalized.dataset, "exprSet") == TRUE){
    data.format <- "exprSet"
    proceed <- TRUE
    if (is.null(dimnames(exprs(unnormalized.dataset))[[1]])){
      warning("Input dataset did not have row names.  ")
      proceed <- FALSE
    }
    if (mode(exprs(unnormalized.dataset)) != "numeric"){
      warning("Input dataset was not numeric.  ")
      proceed <- FALSE
    } 
  }
  if (proceed == FALSE){
    warning("Input dataset was not of class matrix or class exprSet.  ")
  }

  #proceed if no errors were found
  if (proceed == TRUE){

    if (data.format == "matrix"){
      y.i.j <- unnormalized.dataset
      rm(unnormalized.dataset)
    }
    if (data.format == "exprSet"){
      y.i.j <- exprs(unnormalized.dataset)
      exprs(unnormalized.dataset) <- array(0, dim = c(1,1))
    }

    n <- dim(y.i.j)[1]
    p <- dim(y.i.j)[2]
    for (index.j in 1:p){
      y.i.j[,index.j] <- ifelse(is.finite(y.i.j[,index.j]), y.i.j[,index.j], NA)
    }
    for (index.i in 1:n){
      r.i.j.index.i <- rank(y.i.j[index.i,], na.last = "keep")
      p.i.index.i <- length(na.omit(r.i.j.index.i))
      y.i.j[index.i,] <- (r.i.j.index.i - (p.i.index.i + 1)/2)/(p.i.index.i/2)
    }

    #prepare output
    if (data.format == "matrix"){
      normalized.dataset <- y.i.j
      rm(y.i.j)
    }
    if (data.format == "exprSet"){
      normalized.dataset <- unnormalized.dataset
      rm(unnormalized.dataset)
      exprs(normalized.dataset) <- y.i.j
      rm(y.i.j)
    }

    #return result
    normalized.dataset
  }
 
  #in the case of bad input, return NA
  else{
    NA
  }
}

"gr.cv" <-
function(normalized.dataset, query, fun = median){
  result <- c()
  for (index in 1:length(query)){
    left.out <- query[index]
    kept.in <- query
    kept.in[index] <- NA
    kept.in <- as.vector(na.omit(kept.in))
    gr.result <- gr.main(normalized.dataset, kept.in, fun = fun, ngenes = dim(normalized.dataset)[1])
    result <- cbind(result, match(left.out,gr.result[[1]][,1]))
  }
  as.vector(result)
}

require("Biobase")

