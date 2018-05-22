## This example simulation for Raymond et al.'s (2011, Journal of Applied
## Ecology ) Subantarctic Macquarie Island model ("macquarie.dia") 
## uses an increase in rabbits and a decrease in tall tussock vegetation
## under a positive press perturbation to rabbits as the validation
## criterion, and examines outcomes under a perturbation scenario 
## where rabbits, rats and mice decrease.

source("dia.r")
source("community.r")

## Read model specification
edges <- model.dia("macquarie.dia")

## Examine unweighted adjacency matrix
A <- adjacency.matrix(edges)
A

## Function to generate the community matrix
s <- community.sampler(edges)

## Function to check the validation condition
press <- press.validate(edges,
                         perturb=c("Rabbits"=1),
                         monitor=c("Rabbits"=1,"Tall tussock vegetation"=-1))

## Function to define the perturbation scenario
impact <- press.impact(edges,perturb=c("Rabbits"=-1,"Rats"=-1,"Mice"=-1))

## Use 10000 simulations
n.sims <- 10000
results <- 0
i <- 0
while(i < n.sims) {

  ## Randomly choose edges to retain
  z <- s$select(runif(1))
  ## Sample community matrix
  W <- s$community()

  ## Check press condition and stability
  if(!(press(W) && stable.community(W))) next

  ## Monitor impact post press
  imp <- impact(W)
  results <- results + outer(sign(imp),-1:1,'==')
  i <- i+1
}

## Print results
rownames(results) <- levels(edges$From)
colnames(results) <- c('-','0','+')
results

## Plot outcomes
library(RColorBrewer)
pal <- brewer.pal(n=5,"RdBu")[4:2]
opar <- par(mar=c(5,10,1,1)+0.1)
prop <- results/rowSums(results)
r <- colSums(t(prop)*(-1:1))
barplot(t(prop[order(r),]),
		horiz=T,cex.names=0.8,cex.axis=0.8,las=2,border=F,col=pal,xlab="Proportion")
par(opar)
