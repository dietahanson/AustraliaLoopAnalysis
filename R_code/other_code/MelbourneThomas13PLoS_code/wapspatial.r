source("dia.r")
source("community.r")

## Read specification
edges <- model.dia("wapspatial.dia")
edges <- enforce.limitation(edges)

## Examine unweighted adjacency matrix
A <- adjacency.matrix(edges)
A

## Functions to generate the community matrix
s <- community.sampler(edges)

## Function to define the perturbation scenario
impact <- press.impact(edges,perturb=c("Regional warming"=1,"Krill fishery"=1))

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
## Extract the name of the entity and the regional specifier
entity <- gsub(" *\\([MSN]\\)$","",rownames(prop))
region <- gsub(".*\\(([MSN])\\)$|.*","\\1",rownames(prop))
## Enforce ordering
entity <- match(entity,rev(c("Trophic competitors (pelagic)","Trophic competitors","Adelie penguins","Chinstrap penguins","Fish",
                         "Adult krill","Larval krill","Salps",
                         "Large phytoplankton","Small phytoplankton","Krill fishery","Sea ice","Regional warming")))
region <- match(region,c("","S","M","N"))

library(RColorBrewer)
pal <- brewer.pal(n=5,"RdBu")[4:2]

opar <- par(mar=c(5,10,1,1)+0.1)
barplot(t(prop[order(entity,region),]),
        space=c(0,diff(sort(entity)))+0.4,
        horiz=T,cex.names=0.8,cex.axis=1,las=1,border=F,col=pal,xlab="Proportion")
par(opar)


