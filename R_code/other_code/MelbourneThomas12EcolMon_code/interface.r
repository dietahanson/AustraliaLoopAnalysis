## This example of our 'node selector' interface for Hulot et al.'s
## (2000, Nature) eight-variable lake mesocosm model (mesocosm1.dia)
## allows the user to explore prediction space and derive estimates of
## marginal likelihoods for various sets of observations.

source("dia.r")
source("community.r")
source("tk.r")

## Read model specification
edges <- model.dia("mesocosm1.dia")

## Examine unweighted adjacency matrix
A <- adjacency.matrix(edges)
A

## Function to generate the community matrix
s <- community.sampler(edges)

## Use 10000 simulations
n.sims <- 10000
W <- s$community()
nodes <- levels(edges$From)
Ws <- vector("list",n.sims)
n.stable <- 0
while(n.stable < n.sims) {

  ## Randomly choose edges to retain
  z <- s$select(runif(1))
  ## Sample community matrix
  W <- s$community()

  ## Check stability
  if(!stable.community(W)) next
  n.stable <- n.stable+1

  Ws[[n.stable]] <- W
}

node.selector(nodes,Ws)













