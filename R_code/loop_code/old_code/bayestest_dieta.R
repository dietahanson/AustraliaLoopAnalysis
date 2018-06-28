## Code to produce loop models and use the bayes.test function and workflow
## written by Melbourne-Thomas et al. (2012 Ecol Mon) to choose the best model.
## The only two required inputs are x: a table of all interactions, with 3 
## columns; From = first species in interaction, To = second species in 
## interaction, and Type = category of interaction (predatorprey, habitat, etc) 
## and y: a table of models, with each row being a model and each column being 
## an interaction type to include in the model. Interaction types are included 
## with a 1 and excluded with a 0.



##------------------------------------------------------------------------------
## All required functions and packages
##------------------------------------------------------------------------------


library(data.table)

# Function to make a matrix symmetrical (used for trophic matrices)
makeSymm = function(m) {
  m[upper.tri(m)] = t(m)[upper.tri(m)]
  return(m)
}



# Function to populate a trophic matrix with interaction strengths. Requires 
# inputs n = number of nodes, node_names, A = sign matrix, x = master table of 
# interactions
trophic = function(n, node_names, A, x) {    
  
  # first fill the matrix with random beta values
  Bt = makeSymm(matrix(rbeta(n^2,  # make symmetric so ij is a function of ji
                             shape1 = 1, shape2 = 4),  # parameters of beta dist
                       nrow = n))*.99 + 0.01 # make all values between 0.1 and 1  
  colnames(Bt) = node_names
  rownames(Bt) = node_names
  
  # then add signs to strengths
  Ct=Bt*A  
  
  # last, populate diagonal elements 
  for (m in (1:length(node_names))) {
    
    sp = node_names[m]  
    
    # only predators will be in the From column (for predatorprey interactions)
    if (sp %in% x[x$Type=="Predatorprey",]$From) {     
      Ct[sp, sp] = -.1  # if a predator  
    } else { Ct[sp, sp] = -1}  # if basal
  }
  return(Ct)
}



# Function to populate a non-trophic matrix with interaction strengths. Requires 
# inputs n = number of nodes, node_names, A = sign matrix, x = master table of 
# interactions
nontrophic = function(n, node_names, A, x) {  
  
  # first fill the matrix with random uniform values
  Bnt = matrix(runif(n^2), nrow=n)*0.99+0.01  # values between 0.1 and 1
  colnames(Bnt) = node_names
  rownames(Bnt) = node_names
  
  # then add signs to strengths
  Cnt=Bnt*A  
  
  # last, populate diagonal elements
  diag(Cnt) = 0  # diagonal will be populated from the trophic matrix
  return(Cnt)
}



# Function to make a unique identifier for each species pair in edges table
pairnamer = function (df) {return(paste(unlist(sort(df)), collapse = ""))
}



# Function to generate the total community matrix, which is the sum of the 
# various interaction type community matrices required by the model. Requires 
# need = a list of interaction types (predatorprey, habitat, etc.) required by
# the model, n = number of nodes, node_names, am = a list of sign matrices for 
# each interaction type, named by interaction type, and x = master table of 
# interactions
sampler = function(need, n, node_names, am, x) {
  
  sm = list()  # list to contain the interaction type matrices    
  
  # for each interaction type required by the model, populate a matrix with 
  # random values (will use the trophic function if a predatorprey interaction,
  # otherwise will use the nontrophic function) and put into the sm list
  for (g in 1:length(need)) {  
    
    if (grepl("predatorprey", need[g], ignore.case = T)) {
      tm = trophic(n, node_names, am[[need[g]]], x)  # appropriate sign matrix
      sm[[g]] = tm
    } else {
      ntm = nontrophic(n, node_names, am[[need[g]]])
      sm[[g]] = ntm
    }
  } 
  
  # then add all type matrices
  popmod = Reduce("+", sm)  
  
  # now create a table "edges" that has columns; From = first species in 
  # interaction, To = second species in interaction, Type = "N" for negative or
  # "P" for positive, and Pair = unique pair name where ij and ji are both 
  # called ij. Note this is not exactly the same format as the master 
  # interactions table x as it will only contain the interactions called by the 
  # particular model, the Type is only negative or positive, and Pairs is added.
  # Also, each two-way interaction is expanded to two rows, so one predatorprey 
  # interaction in x will here be one positive interaction from prey to predator
  # and one negative interaction for predator to prey
  edges = data.frame(as.table(popmod), stringsAsFactors = F)
  edges$Pair = apply(edges[,c("Var1", "Var2")],1,pairnamer) 
  edges = edges[abs(edges$Freq)>0,]  # we only need non-zero interactions
  setnames(edges,
           old = c("Var1", "Var2", "Freq", "Pair"),
           new = c("From", "To", "Type", "Pair"))
  edges$Type = sapply(edges$Type,
                      function(x) if (x < 0) { x = "N"} else {x = "P"})
  
  # create a list containing the community matrix and the edges table
  res = list(cmat = popmod, edges = edges)
  
}



# Function to estimate posterior probabilities of model correctness, 
# assuming uniform priors. Requires inputs; model = row number of y (model 
# table) for the focal model, perturb = node being perturbed, monitor = node 
# being monitored, n.samples = # of times to create a stable matrix from each 
# model to compare to the model of interest, y = model table
Bayes.test <- function(model,perturb,monitor,n.samples, y) {

    # get a list of the types of interactions the focal model needs
    need = as.numeric(y[model,])
    names(need) = colnames(y)
    need = names(need[need > 0])
    
    
    # simulate a stable weight matrix W for the focal model
    H <- sampler(need, n, node_names, am, x)
    W <- H$cmat
    while(!stable.community(W)) {  
      H <- sampler(need, n, node_names, am, x)
      W <- H$cmat  
    }

    # determine the outcome of the press perturbation. First generate a function
    # to determine the impact (returned function depends on if there are 
    # monitored nodes)
    impact <- press.impact(H$edges,perturb=perturb,monitor=monitor)
    
    # call the function, and get the sign of the impact on the monitored nodes.
    # Call this result (the sign of the impact) "observed"
    observed <- signum(impact(W))
    names(observed) <- names(monitor)

    
    # empty vector for results of model comparison loop below
    prop <- double(nrow(y))    
    
    # for each model, compare the outcome of the press perturbation with the 
    # "observed" outcome for the focal model above 
    for(k in 1:nrow(y)) { 
        n.stable <- 0
        n.valid <- 0
        
        # get a list of the types of interactions this model needs
        need = as.numeric(y[k,])
        names(need) = colnames(y)
        need = names(need[need > 0])

        while(n.stable < n.samples) {

          # simulate a community matrix W
          H = sampler(need, n, node_names, am, x)
          W = H$cmat
          
          # generate a function to see if the impact of a press perturbation on
          # monitored nodes is the same as "observed" above for the focal model
          press = press.validate(H$edges,perturb = perturb, monitor = observed)

            # check stability of W
            if(!stable.community(W)) next  # skip if not stable
            n.stable = n.stable+1  # count if stable

            # check press condition 
            if(press(W)) n.valid = n.valid+1  # count if valid
        }
        # store the ratio of valid matrices to stable and valid matrices for the
        # n.samples of each model
        prop[k] = n.valid/n.stable  
    }
    
    # output the proportion for each model over the sum of proportions for all 
    # models (should be a vector with length = number of models) as p   
    list(p = prop/sum(prop), 
         perturb = perturb,
         monitor = monitor)
}



# Function to check the stability of a simulated community matrix W. If all 
# eigenvalues have a negative real part, matrix is stable
stable.community = function(W) {
  all(Re(eigen(W,symmetric = F,only.values = T)$values)<0)
}



# Function to return sign of s, with values of magnitude less than epsilon
# rounded down to zero
signum = function(s,epsilon = 1.0E-5) {
  (s > epsilon) - (s < -epsilon)
}



# Function to generate a function to check a perturbation outcome against some 
# expectation. Must supply a vector of named elements that specify the relative 
# magnitude of the press perturbation, and a vector of named elements that 
# specify the signs of the change in the monitored nodes
press.validate = function(edges,perturb,monitor,epsilon = 1.0E-5) {
  
  # first make a function that calls the name of a node and retrieves its index
  # in the edges "From" levels
  index = function(name) {   
    k = match(name,levels(edges$From))
    if(any(is.na(k)))
      warning("Unknown nodes:",paste(name[is.na(k)],collapse = " "))
    k
  }
  
  # get the indices of perturbed and monitored nodes
  k.perturb = index(names(perturb))
  k.monitor = index(names(monitor))
  
  # make a vector the length of the edges table and make the the pertubed node 
  # element the negative of the perturbation
  S.press = double(nlevels(edges$From))
  S.press[k.perturb] = -perturb
  monitor = sign(monitor) 
  
  # return a function to check the perturbation condition, and a logical 
  # response asking if the outcome matches the expectation of monitor
  function(W) {
    s = tryCatch(solve(W,S.press),error=function(e) NULL)
    !is.null(s) && all(signum(s[k.monitor],epsilon)==monitor)
  }
}


## Generate a function to determine the impact of a press perturbation
press.impact <- function(edges,perturb,monitor=NULL) {
  
  index <- function(name) {
    k <- match(name,levels(edges$From))
    if(any(is.na(k)))
      warning("Unknown nodes:",paste(name[is.na(k)],collapse=" "))
    k
  }
  
  ## Indices of perturb
  k.perturb <- index(names(perturb))
  S.press <- double(nlevels(edges$From))
  S.press[k.perturb] <- -perturb
  
  if(length(monitor)==0) {
    impact <- function(W) solve(W,S.press)
  } else {
    k.monitor <- index(names(monitor))
    impact <- function(W) solve(W,S.press)[k.monitor]
  }
  
  ## Return function to compute impact
  impact
}



cat(sprintf('Main loop started at: %s\n',date()))

n.sims <- 5
P <- array(0,c(nrow(y),n.sims,nrow(y)))

perturb = setNames(1, "human")
monitor = NULL
for(model in 1:nrow(y)) {
  for(k in 1:n.sims) {
    ## Compute posterior probabilities
    fit <- Bayes.test(model,perturb,monitor,n.samples=10, y)
    P[,k,model] <- fit$p
  }
}

## Generate a confusion matrix
CM <- matrix(0,dim(P)[3],dim(P)[1])
for(model in 1:(dim(P)[3]))
  CM[model,] <- tabulate(unlist(apply(P[,,model],2,function(x) which(x==max(x)))),dim(P)[1])
## Columns: correct model, Rows: highest ranking model
t(CM)

## Determine how often the correct model is most likely, second most
## likely etc.
Q <- matrix(0,dim(P)[3],dim(P)[1])
for(model in 1:(dim(P)[3]))
  Q[model,] <- tabulate(apply(-P[,,model],2,
	function(x) which(order(x,c(model,rep(0,dim(P)[1]-1)))==model)),dim(P)[1])

## Plot how often the correct model is most likely, second most likely
## etc, separately for each correct model and overall.
opar <- par(mfrow=c(3,2))
for(model in 1:nrow(Q)) {
  p <- Q[model,]
  names(p) <- 1:ncol(Q)
  barplot(p/sum(p),main=model,ylim=c(0,1))
}
p <- colSums(Q)
names(p) <- 1:ncol(Q)
barplot(p/sum(p),main="all",ylim=c(0,1))

## Use clustering to show model similarity
opar <- par(mfrow=c(3,2))
method <- "complete"
for(model in 1:(dim(P)[3])) {
  plot(hclust(dist(P[,,model]),method=method),ann=F)
  title(model)
}
plot(hclust(dist(matrix(P,dim(P)[1],prod(dim(P)[-1]))),method=method),ann=F)
title("All")
par(opar)

