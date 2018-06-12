source("dia.r")
source("community.r")

## Read specification
edges <- model.dia("wapextended.dia")
edges <- enforce.limitation(edges)

## Examine unweighted adjacency matrix
A <- adjacency.matrix(edges)
A

## Functions to generate the community matrix
s <- community.sampler(edges)

## Function to define the perturbation scenario
impact <- press.impact(edges,perturb=c("Warming"=1,"Krill fishery"=1))

n.sims <- 10000
results <- 0
i <- 0
while(i < n.sims) {

  ## Sample community matrix
  z <- s$select(runif(1))
  W <- s$community()

  ## Check press condition and stability
  if(!stable.community(W)) next

  ## Monitor impact post press
  imp <- impact(W)
  results <- results + outer(sign(imp),-1:1,'==')
  i <- i+1
}

rownames(results) <- levels(edges$From)
colnames(results) <- c('-','0','+')
results

## Proportions
prop <- results/rowSums(results)
## Extract the name of the entity
entity <- gsub(" *\\([MSN]\\)$","",rownames(prop))
## Enforce ordering
entity <- match(entity,c("Warming","Sea ice","Krill fishery","Small phytoplankton","Large phytoplankton",
                         "Salps","Larval krill","Adult krill","Fish",
                         "Chinstrap penguins","Adelie penguins","Trophic competitors"))

library(RColorBrewer)
pal <- brewer.pal(n=5,"RdBu")[4:2]

opar <- par(mar=c(5,10,1,1)+0.1)
barplot(t(prop[order(entity),]),
        horiz=T,cex.names=1.1,cex.axis=1.1,cex.lab=1.1,las=1,border=F,col=pal,xlab="Proportion")
par(opar)
