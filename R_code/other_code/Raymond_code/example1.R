## This is a simplified example of the modelling technique used in:
## Raymond, B. et al. (2010) Qualitative modelling of invasive species eradication on subantarctic Macquarie Island. Journal of Applied Ecology
##
## In this example, we ignore "uncertain" links in our model, and evaluate only one model (using just the links that we are certain about).
## This demonstrates the method used to evaluate a single model. See example2.R to see how this is embedded within the evaluation of
##  structural uncertainty (i.e. taking "uncertain" links into account).

##------------------------------------------------------------------------------
## read in the interactions data
## adjust the path here to match where you put the data file

x <- read.table("macquarie_interactions.csv", sep=",", header=TRUE)


##------------------------------------------------------------------------------
## initialise variables and settings

nwrand <- 1000 ## number of randomisations with different random weights on interaction strengths
max_nwrand <- 20*nwrand ##try a maximum of this many realisations to get the nwrand results


## validation data: what can we use to ground-truth our models?
## these are responses to suppression of cats and myxoma (i.e. rabbits increased, tall tussock decreased, and cats decreased)
model_validation <- matrix(c("Rabbits", 1, "Tall tussock vegetation", -1, "Cats", -1), nrow=3, ncol=2, byrow=TRUE)
colnames(model_validation) <- c("response_node", "response_value")

##> model_validation
##     response_node             response_value
##[1, ] "Rabbits"                 "1"
##[2, ] "Tall tussock vegetation" "-1"
##[3, ] "Cats"                    "-1"


## get the list of unique names within these interactions
node_names <- unique(union(x$To, x$From))

## initialise some variables that will keep track of various things as we go
##  suffixes of _cm indicate that this variable applies to the cat/myxo suppression simulations
##  suffixes of _pcm indicate that this variable applied to the eradication project simulations (i.e. post cats and myxo)

stable_count_cm <- 0 ## number of realisations that were stable

## summary of predictions for each element in the model
rsummary_cm <- matrix(0, nrow=3, ncol=length(node_names))
rownames(rsummary_cm) <- c("Negative", "Zero", "Positive")
colnames(rsummary_cm) <- node_names

cm_idx <- grep("(cats|myxoma)", node_names, ignore.case=TRUE)
pcm_idx <- grep("[^(cats|myxoma)]", node_names, ignore.case=TRUE) ##which elements in the model are relevant for the eradication project (not cats or myxo)

pcm_node_names <- node_names[pcm_idx]

stable_count_pcm <- 0 ## number of realisations that passed the additional stability test (after the removal of cats and myxo from the model)

pcm_results_cols <- grep("(rabbits|rats|mice)", node_names, ignore.case=TRUE) ## the columns in the A matrix that give the responses to perturbations of our eradication target species

## summary of predictions for each element in the model
rsummary_pcm <- matrix(0, nrow=3, ncol=length(pcm_node_names))
rownames(rsummary_pcm) <- c("Negative", "Zero", "Positive")
colnames(rsummary_pcm) <- pcm_node_names

## summary of predictions for each element in the model, but only for the subset of simulations in which the target species were actually suppressed
rsummary_pcm_subset <- matrix(0, nrow=3, ncol=length(pcm_node_names)) ##summary for successful runs
rownames(rsummary_pcm_subset) <- c("Negative", "Zero", "Positive")
colnames(rsummary_pcm_subset) <- pcm_node_names

diagval_max <- -0.25 ## maximum allowable value for self-limitation links

n <- length(node_names) ##number of nodes
n_pcm <- length(pcm_node_names) ##number of nodes once cats and myxo removed


##------------------------------------------------------------------------------
## do the actual computations

cat(sprintf("Main loop started at: %s\n", date()))

## initialise the interactions matrix A
A <- matrix(0, nrow=n, ncol=n)
colnames(A) <- node_names
rownames(A) <- node_names

## computations expect convention A[i, j] is the effect on i from j

for (k in (1:dim(x)[1])) {  ##from 1 to the total number of links
  if (grepl("(unknown|minor)", x[k, ]$Importance, ignore.case=TRUE)) {
    ## ignore "minor" or "unknown" links
    next
  }
  this_from <- as.character(x[k, ]$From) ## the "from" element of the interaction
  this_to <- as.character(x[k, ]$To) ## the "to" element of the interaction
  if (grepl("competition", x[k, ]$Type, ignore.case=TRUE)) {
    A[this_from, this_to] <- -1
    A[this_to, this_from] <- -1
  } else if (grepl("habitat", x[k, ]$Type, ignore.case=TRUE) || grepl("positive", x[k, ]$Type, ignore.case=TRUE)) {
    A[this_to, this_from] <- 1
  } else if (grepl("limiting", x[k, ]$Type, ignore.case=TRUE)) {
    A[this_to, this_from] <- -1
  } else if (grepl("negative", x[k, ]$Type, ignore.case=TRUE)) {
    A[this_to, this_from] <- -1
  } else if (grepl("predator-prey", x[k, ]$Type, ignore.case=TRUE)) {
    A[this_from, this_to] <- 1
    A[this_to, this_from] <- -1
  } else if (grepl("scavenging", x[k, ]$Type, ignore.case=TRUE)) {
    A[this_from, this_to] <- 1
  } else {
    stop(sprintf("unrecognised link type %s (%s to %s)", x[k, ]$Type, this_from, this_to))
  }
}

this_valid_count <- 0

## now loop many times, each time generating new (random) weights in A
for (twi in 1:max_nwrand) {
  ## randomisation of weights in A, in range 0.01-1
  wA <- matrix(runif(n^2), nrow=n)*0.99+0.01
  wA <- wA*sign(A)
  ## but make sure diagonals are between -1 and diagval_max
  temp <- diag(rep(1, dim(A)[1]))
  wA <- wA*(1-temp)+(-1-diagval_max)*diag(runif(n))+diagval_max*temp

  ## calculate the response, and check the validation criteria
  adjwA <- -solve(wA) ##negative inverse of wA
  adjwA[abs(adjwA)<1e-07] <- 0 ## set very small responses to zero

  ## check that response fits validation data
  this_valid <- TRUE
  ## check all validation criteria: if any fail, then this model is not valid
  for (k in 1:dim(model_validation)[1]) {
    temp_response <- -adjwA[, cm_idx]
    if (dim(temp_response)[2]>1) {
      ## sum response over the individual responses
      temp_response <- rowSums(temp_response)
    }
    temp_response <- sign(temp_response)
    if (temp_response[model_validation[k, "response_node"]] != model_validation[k, "response_value"]) {
      this_valid <- FALSE ##not valid
      break
    }
  }

  if (!this_valid) next ## not valid, so discard this realisation

  ## ok, this realisation is valid: is it also stable?
  if (!all(Re(eigen(wA, only.values=TRUE)$values)<0)) {
    ## not stable
    next
  }

  ## so if we got this far, this realisation is both valid and stable
  stable_count_cm <- stable_count_cm+1
  this_valid_count <- this_valid_count+1

  ## find the predicted response to cat/myxo suppression
  temp <- -adjwA[, cm_idx]
  if (dim(temp)[2]>1)   temp <- rowSums(temp)    ## sum response over the individual responses (i.e. total response is the sum of responses to cat and myxo suppression)

  temp <- sign(temp) ## only interested in signs of responses

  ## keep tabs on summary results
  rsummary_cm[1, ] <- rsummary_cm[1, ]+as.double(temp==-1)
  rsummary_cm[2, ] <- rsummary_cm[2, ]+as.double(temp==0)
  rsummary_cm[3, ] <- rsummary_cm[3, ]+as.double(temp==1)


  ## now simulate the eradication project

  wA <- wA[pcm_idx, pcm_idx] ##drop the columns and rows associated with cats and myxo

  ##check stability of model with cats and myxo and their associated links removed
  this_pcm_stable <- 1
  if (!all(Re(eigen(wA, only.values=TRUE)$values)<0)) {
    this_pcm_stable <- 0
  }

  if (this_pcm_stable) {
    stable_count_pcm <- stable_count_pcm+1

    adjwA <- -solve(wA) ##negative inverse of wA
    adjwA[abs(adjwA)<1e-07] <- 0 ## set very small responses to zero

    ## predicted response to eradication project
    temp <- -adjwA[, pcm_results_cols]
    if (dim(temp)[2]>1) temp <- rowSums(temp)
    temp <- sign(temp)

    rsummary_pcm[1, ] <- rsummary_pcm[1, ]+as.double(temp==-1)
    rsummary_pcm[2, ] <- rsummary_pcm[2, ]+as.double(temp==0)
    rsummary_pcm[3, ] <- rsummary_pcm[3, ]+as.double(temp==1)

    ## also keep separate count of results for models in which target species were actually suppressed
    if (all(temp[pcm_results_cols]<0)) {
      rsummary_pcm_subset[1, ] <- rsummary_pcm_subset[1, ]+as.double(temp==-1)
      rsummary_pcm_subset[2, ] <- rsummary_pcm_subset[2, ]+as.double(temp==0)
      rsummary_pcm_subset[3, ] <- rsummary_pcm_subset[3, ]+as.double(temp==1)
    }
  } ## if this_pcm_stable
  if (this_valid_count==nwrand) {
    ## we are trying up to max_nwrand times, but if we"ve achieved our target of nwrand valid configurations, bail out
    break
  }
} ##end twi loop

cat(sprintf("\n"))


##------------------------------------------------------------------------------
## show some outputs

## outcomes of cat/myxo suppression
par(mfrow=c(2, 1), fig=c(0.02, 0.98, 0.2, 0.8))
barplot(rsummary_cm/stable_count_cm, cex.names=0.7, las=2, legend=TRUE, xpd=TRUE, args.legend=list(horiz=TRUE, y=1.25), ylab="Proportion")
title("Cat/myxoma suppression", line=3)

## outcomes of eradication project
## all simulation runs
par(mfrow=c(2, 1), fig=c(0.02, 0.98, 0.2, 0.8))
barplot(rsummary_pcm/stable_count_pcm, cex.names=0.7, las=2, legend=TRUE, xpd=TRUE, args.legend=list(horiz=TRUE, y=1.25), ylab="Proportion")
title("Eradication project predictions", line=3)

## also show only those in which the target species were actually suppressed
par(mfrow=c(2, 1), fig=c(0.02, 0.98, 0.2, 0.8))
barplot(rsummary_pcm_subset/unique(colSums(rsummary_pcm_subset)), cex.names=0.7, las=2, legend=TRUE, xpd=TRUE, args.legend=list(horiz=TRUE, y=1.25), ylab="Proportion")
title("Eradication project predictions (subset)", line=3)

