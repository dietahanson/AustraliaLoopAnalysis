## This example uses boosted regression trees to examine the relative
## importance of model edges (linkages) in determining the response
## of Alg1 to and increase in Phos for simulations of Hulot et al.'s
## (2000, Nature) eight-variable lake mesocosm model (mesocosm1.dia).

source('dia.r')
source('community.r')

## Read model specification
edges <- model.dia('mesocosm1.dia')

## Examine unweighted adjacency matrix
A <- adjacency.matrix(edges)
A

## Function to generate the community matrix
s <- community.sampler(edges)

## Function to define the perturbation scenario
impact <- press.impact(edges,perturb=c("Phos"=1))

## Keep track of edge weights in data frame
Wlog <- matrix(nrow=1,ncol=length(edges$From))
Wlog <- as.data.frame(Wlog)
whichW <- which(A!=0)
tempedgenames <- matrix("",nrow=dim(A)[1],ncol=dim(A)[2])
for (i in 1:dim(edges)[1]) {
    tempedgenames[edges$To[i],edges$From[i]] <- paste(edges$From[i],edges$To[i],sep='_')
}
colnames(Wlog) <- tempedgenames[whichW]

## Use 1000 simulations
n.sims <- 1000
n.stable <- 0
allimp <- matrix(0,n.sims,length(levels(edges$From)))
while(n.stable < n.sims) {

  ## Randomly choose edges to retain
  z <- s$select(runif(1))
  ## Sample community matrix
  W <- s$community()

  ## Check stability
  if(!stable.community(W)) next
  n.stable <- n.stable+1

  ## Monitor impact post press
  imp <- impact(W)
  allimp[n.stable,] <- imp
  Wlog[n.stable,] <- W[whichW]
}

## Fit classification tree model (with Alg1 as the response)
library(dismo)
temp <- sign(allimp)
colnames(temp) <- levels(edges$From)
temp <- cbind(temp,Wlog)
d <- cbind(temp["Alg1"],temp[names(Wlog)])
fit3 <- gbm.step(data=d,gbm.x=2:24,gbm.y=1,family="gaussian",tree.complexity=5,
				learning.rate = 0.01,bag.fraction = 0.5)
dev.new()
opar <- par(mar=c(5,10,1,1)+0.1)
summary(fit3,las=2)
dev.new()
gbm.plot(fit3, n.plots=12, write.title = FALSE)

