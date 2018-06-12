###############################################################################
# Loop model based on Raymond's code, with the model predicting the effect of
# removing humans. 
###############################################################################



#------------------------------------------------------------------------------
# read in the interactions table
#------------------------------------------------------------------------------

# adjust the path here to match where you put the data file
setwd("~/Documents/Australia/R_code/loop_code/")
x = read.table('hotgrouped.csv',
             sep=',',
             header=T)


#------------------------------------------------------------------------------
# validation data: what can we use to ground-truth our models?
#------------------------------------------------------------------------------

# these are responses to removal of humans 
model_validation = matrix(c('human',-1,'fox',1,'camel',1),
                        nrow = 3,
                        ncol = 2,
                        byrow = T)

colnames(model_validation) = c('response_node','response_value')

model_validation
# response_node response_value
# [1,] "human"       "-1"           
# [2,] "fox"         "1"           
# [3,] "camel"      "1"         


#------------------------------------------------------------------------------
# initialise some variables and settings
#------------------------------------------------------------------------------

nwrand = 1000  # number of randomisations with different interaction strengths
max_nwrand = 20*nwrand  # try a maximum of this many realisations 

node_names = unique(union(x$To,x$From))  # get the list of nodes
n = length(node_names)  # number of nodes

stable_count_h = 0  # number of realisations that were stable

h_idx = grep('(human)',node_names,ignore.case=T)  # just humans

# summary of predictions for each element in the model
rsummary_h = matrix(0,nrow=3,ncol=length(node_names))
rownames(rsummary_h) = c('Negative','Zero','Positive')
colnames(rsummary_h) = node_names



#------------------------------------------------------------------------------
# do the actual computations

cat(sprintf('Main loop started at: %s\n',date()))
flush.console()

# initialise the interactions matrix A (signs) for trophic interactions   
Atro=matrix(0,nrow=n,ncol=n)
colnames(Atro)=node_names
rownames(Atro)=node_names

# initialise the interactions matrix B (strengths) for trophic interactions   
Btro=matrix(0,nrow=n,ncol=n)
colnames(Btro)=node_names
rownames(Btro)=node_names

# set conversion efficiency (how much a predator gets from a prey)

e <- 0.1

# computations expect convention A[i,j] is the effect on i row from j column

for (k in (1:dim(x)[1])) {
  
  this_from=as.character(x[k,]$From) # the "from" element of the interaction
  this_to=as.character(x[k,]$To) # the "to" element of the interaction
  if (grepl('predator-prey',x[k,]$Type,ignore.case=T)) {
    Atro[this_from,this_to]=0.1
    Atro[this_to,this_from]=-1
  } else {
    stop(sprintf('unrecognised link type %s (%s to %s)',x[k,]$Type,this_from,this_to))
  }
}

this_valid_count=0

# now loop, each time generating new (random) weights in Btro
for (twi in 1:max_nwrand) {
  
  makeSymm <- function(m) {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }
  
  #Btro <- matrix(rbeta(n^2, shape1 = 1, shape2 = 4), nrow=n) #not symmetric
  Btro <- makeSymm(matrix(rbeta(n^2,
                                shape1 = 1,
                                shape2 = 4),
                          nrow=n))*.99+0.01 #  symmetric
  colnames(Btro)=node_names
  rownames(Btro)=node_names
  

  Ctro=Btro*Atro
  
  for (m in (1:length(node_names))) {
    sp <- node_names[m] # the "from" element of the interaction

    if (sp %in% x$From) {
      Ctro[sp, sp] <- -.001 
    } else { Ctro[sp, sp] <- -1}
  }

  # calculate the response, and check the validation criteria
  adjwA=-solve(Ctro) #negative inverse of wA
  adjwA[abs(adjwA)<1e-07]=0 # set very small responses to zero
  
  # check that response fits validation data
  this_valid=1
  # check all validation criteria: if any fail, then this model is not valid
  for (k in 1:dim(model_validation)[1]) {
    temp_response=-adjwA[,h_idx]
    if (is.matrix(temp_response)) {  #had to change this--previous code didn't work if only one column (if only one taxon is removed)
      # sum response over the individual responses
      temp_response=rowSums(temp_response)
    } else {temp_response=temp_response}
    temp_response=sign(temp_response)
    if (temp_response[model_validation[k,'response_node']] != model_validation[k,'response_value']) {
      this_valid=0 #not valid
      break
    }
  }

  if (!this_valid) next # not valid, so discard this realisation

  this_valid_count=this_valid_count+1
  
  
  # ok, this realisation is valid: is it also stable?
  if (!all(Re(eigen(Ctro,only.values=T)$values)<0)) {
    # not stable
    next
  }  
  
  # so if we got this far, this realisation is stable
  stable_count_h=stable_count_h+1

  
  # find the predicted response to cat/myxo suppression
  temp=-adjwA[,h_idx] 
  if (is.matrix(temp))   temp=rowSums(temp)    # sum response over the individual responses (i.e. total response is the sum of responses to cat and myxo suppression)
                                                #had to change this as above in cases where there is only one taxa removed (human)
  temp=sign(temp) # only interested in signs of responses
  
  # keep tabs on summary results
  rsummary_h[1,]=rsummary_h[1,]+as.double(temp==-1)
  rsummary_h[2,]=rsummary_h[2,]+as.double(temp==0)
  rsummary_h[3,]=rsummary_h[3,]+as.double(temp==1)
  
  # if (this_valid_count==nwrand) {
  #   # we are trying up to max_nwrand times, but if we've achieved our target of nwrand valid configurations, bail out
  #   break
  # }
} #end twi loop

cat(sprintf('\n'))
flush.console()


#------------------------------------------------------------------------------
# show some outputs

# outcomes of human suppression

dddd <- gsub("-", "", as.character(Sys.Date()))

v<-paste("Number of stable models:", stable_count_h,
         "Number of valid models:", this_valid_count,
         "Number of models tried:", max_nwrand,
         "\nConsumer Limitation:", max(diag(Ctro)),
         "Basal Limitation:", min(diag(Ctro)))

pdf(paste("australia_loop_outcome",dddd,max(diag(Ctro)), ".pdf", sep=''),  width=11, height=8.5)
par(mfrow=c(2,1),fig=c(0.02,0.98,0.2,0.8), mar=c(4,4,4,7), xpd=T)
barplot(rsummary_h/stable_count_h,
        cex.names=0.7,
        las=2,
        legend=T,
        xpd=T,
        args.legend=list(x = 'topright', inset = c(-.1,0)),
        ylab="Proportion")


title(main = v, cex.main= .9, line = 2)
print(v)

dev.off()

