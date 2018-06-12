## This example demonstrates the calculation of mutual information
## for two alternative versions of Hulot et al.'s (2000, Nature) eight-variable
## lake mesocosm model ("mesocosm1.dia" and "mesocosm2.dia") under a 
## perturbation scenario of increased phosphorus.

source("dia.r")
source("community.r")

files <- c("mesocosm1.dia","mesocosm2.dia")

## Read model specifications
edges <- lapply(files,model.dia)

## Create adjacency matrices
A <- lapply(edges,adjacency.matrix)

## Function to generate the community matrices
sampler <- lapply(edges,community.sampler)

## Function to define the perturbation scenario
impact <- lapply(edges,press.impact,perturb=c("Phos"=1))

labels <- node.labels(edges[[1]])
n.models <- length(files)

## Use 100000 simulations
n.samples <- 100000

## Model indicator
model <- integer(n.samples)

## Outcomes
impacts <- matrix(0,n.samples,length(labels))
accepted <- 0
while(accepted < n.samples) {

  ## Sample community matrix
  m <- sample(n.models,1)
  W <- sampler[[m]]$community()

  ## Check press condition and stability
  if(!stable.community(W)) next
  accepted <- accepted+1

  ## Monitor impact post press
  imp <- impact[[m]](W)
  model[accepted] <- m
  impacts[accepted,] <- imp
}
colnames(impacts) <- labels

## Mutual information of each node with the model indicator
MI <- apply(signum(impacts),2,function(x) mutual.info(model,x))
round(MI[order(MI)],digits=3)

## Posterior probabilities of model correctness
table(model)/n.samples