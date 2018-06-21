## This example of our Bayes test for model correctness uses the
## set of 5 five-variable models (fivevariable1.dia - fivevariable5.dia),
## originally developed by Dambacher et al. (2003, Ecologial
## Modelling).


source("community.r")
source("preamble.R")



## Function to estimate posterior probabilities of model correctness, 
## assuming uniform priors
Bayes.test <- function(model,perturb,monitor,n.samples=1000, y) {

  
    # get a list of the types of interactions this model needs
    need = as.numeric(y[model,])
    names(need) = colnames(y)
    need = names(need[need > 0])
    
    
    ## Simulate a stable weight matrix from the selected model
    H <- sampler(need, n, node_names, am, x)
    W <- H$cmat
    while(!stable.community(W)) {  
      H <- sampler(need, n, node_names, am, x)
      W <- H$cmat  
    }

    ## Determine the outcome of the press perturbation
    impact <- press.impact(H$edges,perturb=perturb,monitor=monitor)
    observed <- signum(impact(W))
    names(observed) <- names(monitor)

    prop <- double(nrow(y)) #edges here is a list, so length is # of models
    
    ## Loop over models
    for(k in 1:nrow(y)) { #edges here is a list, so length is # of models
        n.stable <- 0
        n.valid <- 0
        
        # get a list of the types of interactions this model needs
        need = as.numeric(y[k,])
        names(need) = colnames(y)
        need = names(need[need > 0])



        while(n.stable < n.samples) {

            ## Sample community matrix
          H <- sampler(need, n, node_names, am, x)
          W <- H$cmat
          
          ## Function to validate press condition
          press <- press.validate(H$edges,perturb=perturb, monitor=observed)

            ## Check stability
            if(!stable.community(W)) next
            n.stable <- n.stable+1

            ## Check press condition
            if(press(W)) n.valid <- n.valid+1
        }
        prop[k] <- n.valid/n.stable
    }
    list(p=prop/sum(prop),
         perturb=perturb,
         monitor=monitor)
}


n.sims <- 5
P <- array(0,c(nrow(y),n.sims,nrow(y)))

perturb = setNames(1, "human")
monitor = setNames(1, "human")
for(model in 1:5) {
  for(k in 1:n.sims) {
    ## Compute posterior probabilities
    fit <- Bayes.test(model,perturb,monitor,n.samples=10, y)
    P[,k,model] <- fit$p
  }
}

# ## Generate a confusion matrix
# CM <- matrix(0,dim(P)[3],dim(P)[1])
# for(model in 1:(dim(P)[3]))
#   CM[model,] <- tabulate(unlist(apply(P[,,model],2,function(x) which(x==max(x)))),dim(P)[1])
# ## Columns: correct model, Rows: highest ranking model
# t(CM)
# 
# ## Determine how often the correct model is most likely, second most
# ## likely etc.
# Q <- matrix(0,dim(P)[3],dim(P)[1])
# for(model in 1:(dim(P)[3]))
#   Q[model,] <- tabulate(apply(-P[,,model],2,
# 	function(x) which(order(x,c(model,rep(0,dim(P)[1]-1)))==model)),dim(P)[1])
# 
# ## Plot how often the correct model is most likely, second most likely
# ## etc, separately for each correct model and overall.
# opar <- par(mfrow=c(3,2))
# for(model in 1:nrow(Q)) {
#   p <- Q[model,]
#   names(p) <- 1:ncol(Q)
#   barplot(p/sum(p),main=model,ylim=c(0,1))
# }
# p <- colSums(Q)
# names(p) <- 1:ncol(Q)
# barplot(p/sum(p),main="all",ylim=c(0,1))
# 
# ## Use clustering to show model similarity
# opar <- par(mfrow=c(3,2))
# method <- "complete"
# for(model in 1:(dim(P)[3])) {
#   plot(hclust(dist(P[,,model]),method=method),ann=F)
#   title(model)
# }
# plot(hclust(dist(matrix(P,dim(P)[1],prod(dim(P)[-1]))),method=method),ann=F)
# title("All")
# par(opar)
# 
