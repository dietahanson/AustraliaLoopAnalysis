###############################################################################
# Loop model based on Raymond's code, with the model predicting the effect of
# removing humans. 
###############################################################################


setwd("~/Documents/Australia/R_code")


#------------------------------------------------------------------------------
# read in the interactions data
# adjust the path here to match where you put the data file

x=read.table('hotgrouped.csv',sep=',',header=T)


#------------------------------------------------------------------------------
# initialise variables and settings

nwrand=1000 # number of randomisations with different interaction strengths
max_nwrand=20*nwrand #try a maximum of this many realisations 

# validation data: what can we use to ground-truth our models?
# these are responses to removal of humans 
model_validation=matrix(c('human',-1,'fox',1,'camel',1),
                        nrow=3,ncol=2,byrow=TRUE)

colnames(model_validation)=c('response_node','response_value')

model_validation
# response_node response_value
# [1,] "human"       "-1"           
# [2,] "fox"         "1"           
# [3,] "camel"      "1"         


# get the list of unique names within these interactions
node_names=unique(union(x$To,x$From))

# initialise some variables that will keep track of various things as we go
#  suffixes of _h indicate that this variable applies to the human model
#  suffixes of _ph indicate that this variable applies to the post human model

stable_count_h=0 # number of realisations that were stable

# summary of predictions for each element in the model

rsummary_h=matrix(0,nrow=3,ncol=length(node_names))
rownames(rsummary_h)=c('Negative','Zero','Positive')
colnames(rsummary_h)=node_names

h_idx=grep('(human)',node_names,ignore.case=T) 
ph_idx=grep('[^(human)]',node_names,ignore.case=T) # all node names except humans

ph_node_names=node_names[ph_idx]

stable_count_ph=0 # number of models that passed the stability test post humans 

ph_results_cols=grep('(human|fox|camel)',node_names,ignore.case=T) # the columns in the A matrix that give the responses to perturbations of our eradication target species

# summary of predictions for each element in the model
rsummary_ph=matrix(0,nrow=3,ncol=length(ph_node_names))
rownames(rsummary_ph)=c('Negative','Zero','Positive')
colnames(rsummary_ph)=ph_node_names

# # summary of predictions for each element in the model, but only for the subset of simulations in which the target species were actually suppressed
rsummary_ph_subset=matrix(0,nrow=3,ncol=length(ph_node_names)) #summary for successful runs
rownames(rsummary_ph_subset)=c('Negative','Zero','Positive')
colnames(rsummary_ph_subset)=ph_node_names



n=length(node_names) #number of nodes
n_ph=length(ph_node_names) #number of nodes once humans removed


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
  if (grepl('competition',x[k,]$Type,ignore.case=T)) {
    Atro[this_from,this_to]=-1
    Atro[this_to,this_from]=-1
  } else if (grepl('habitat',x[k,]$Type,ignore.case=T) || grepl('positive',x[k,]$Type,ignore.case=T)) {
    Atro[this_to,this_from]=1
  } else if (grepl('limiting',x[k,]$Type,ignore.case=T)) {
    Atro[this_to,this_from]=-1
  } else if (grepl('negative',x[k,]$Type,ignore.case=T)) {
    Atro[this_to,this_from]=-1
  } else if (grepl('predator-prey',x[k,]$Type,ignore.case=T)) {
    Atro[this_from,this_to]=1
    Atro[this_to,this_from]=-1
  } else if (grepl('scavenging',x[k,]$Type,ignore.case=T)) {
    Atro[this_from,this_to]=1
  } else {
    stop(sprintf('unrecognised link type %s (%s to %s)',x[k,]$Type,this_from,this_to))
  }
}

this_valid_count=0

# now loop many times, each time generating new (random) weights in Btro
for (twi in 1:max_nwrand) {
  
  for (p in (1:dim(x)[1])) {
    
    this_from=as.character(x[p,]$From) # the "from" element of the interaction
    this_to=as.character(x[p,]$To) # the "to" element of the interaction
    
    
    # randomisation of weights in Btro
    
    if (grepl('predator-prey',x[p,]$Type,ignore.case=T)) {
      Btro[this_to,this_from] = rbeta(1, shape1 = 1, shape2 = 4)
      Btro[this_from,this_to] = e * Btro[this_to,this_from]
       } else {
        Btro[this_from, this_to] = 0
        Btro[this_to, this_from] = 0
      }
    }
  

  comtro=Btro*sign(Atro)
  
  
  # make diagonals -1 for basal species and -0.001 for consumers 
  
  for (m in (1:length(node_names))) {
    
    sp <- node_names[m] # the "from" element of the interaction
    
    if (sp %in% x$From) {
      comtro[sp, sp] <- -0.001
    } else { comtro[sp, sp] <- -1}
  }
  
  
  # calculate the response, and check the validation criteria
  adjwA=-solve(wA) #negative inverse of wA
  adjwA[abs(adjwA)<1e-07]=0 # set very small responses to zero
  
  # check that response fits validation data
  this_valid=1
  # check all validation criteria: if any fail, then this model is not valid
  for (k in 1:dim(model_validation)[1]) {
    temp_response=-adjwA[,h_idx]
    if (is.matrix(temp_response)) {  #had to change this--previous code didn't work if only one column (if only one taxon is removed)
      # sum response over the individual responses
      temp_response=rowSums(temp_response)
    }
    temp_response=sign(temp_response)
    if (temp_response[model_validation[k,'response_node']] != model_validation[k,'response_value']) {
      this_valid=0 #not valid
      break
    }
  }

  if (!this_valid) next # not valid, so discard this realisation

  # ok, this realisation is valid: is it also stable?
  if (!all(Re(eigen(wA,only.values=T)$values)<0)) {
    # not stable
    next
  }  
  
  # so if we got this far, this realisation is stable
  stable_count_h=stable_count_h+1
  this_valid_count=this_valid_count+1
  
  # find the predicted response to cat/myxo suppression
  temp=-adjwA[,h_idx] 
  if (is.matrix(temp))   temp=rowSums(temp)    # sum response over the individual responses (i.e. total response is the sum of responses to cat and myxo suppression)
                                                #had to change this as above in cases where there is only one taxa removed (human)
  temp=sign(temp) # only interested in signs of responses
  
  # keep tabs on summary results
  rsummary_h[1,]=rsummary_h[1,]+as.double(temp==-1)
  rsummary_h[2,]=rsummary_h[2,]+as.double(temp==0)
  rsummary_h[3,]=rsummary_h[3,]+as.double(temp==1)
  
  
  # # now simulate the eradication project
  # 
  # wA=wA[ph_idx,ph_idx] #drop the columns and rows associated with cats and myxo
  # 
  # #check stability of model with cats and myxo and their associated links removed
  # this_pcm_stable=1
  # if (!all(Re(eigen(wA,only.values=T)$values)<0)) {
  #   this_pcm_stable=0
  # }  
  # 
  # if (this_pcm_stable) {
  #   stable_count_pcm=stable_count_pcm+1
  #   
  #   adjwA=-solve(wA) #negative inverse of wA
  #   adjwA[abs(adjwA)<1e-07]=0 # set very small responses to zero
  #   
  #   # predicted response to eradication project
  #   temp=-adjwA[,pcm_results_cols]
  #   if (dim(temp)[2]>1) temp=rowSums(temp)
  #   temp=sign(temp)
  #   
  #   rsummary_pcm[1,]=rsummary_pcm[1,]+as.double(temp==-1)
  #   rsummary_pcm[2,]=rsummary_pcm[2,]+as.double(temp==0)
  #   rsummary_pcm[3,]=rsummary_pcm[3,]+as.double(temp==1)
  #   
  #   # also keep separate count of results for models in which target species were actually suppressed
  #   if (all(temp[pcm_results_cols]<0)) {
  #     rsummary_pcm_subset[1,]=rsummary_pcm_subset[1,]+as.double(temp==-1)
  #     rsummary_pcm_subset[2,]=rsummary_pcm_subset[2,]+as.double(temp==0)
  #     rsummary_pcm_subset[3,]=rsummary_pcm_subset[3,]+as.double(temp==1)
  #   }          
  # } # if this_pcm_stable
  # if (this_valid_count==nwrand) {
  #   # we are trying up to max_nwrand times, but if we've achieved our target of nwrand valid configurations, bail out
  #   break
  # }
} #end twi loop

cat(sprintf('\n'))
flush.console()


#------------------------------------------------------------------------------
# show some outputs

# outcomes of cat/myxo suppression
dev.new()
dddd <- gsub("-", "", as.character(Sys.Date()))

pdf(paste("australia_loop_outcome",dddd,".pdf", sep=''),  width=11, height=8.5)
par(mfrow=c(2,1),fig=c(0.02,0.98,0.2,0.8))
barplot(rsummary_h/stable_count_h,cex.names=0.7,las=2,legend=T,xpd=T,args.legend=list(horiz=T,y=1.25),ylab="Proportion")
title("Human suppression",line=3)
paste(cat("Number of valid and stable models:", stable_count_h, "Number of models tried:", max_nwrand))

dev.off()

# # outcomes of eradication project
# # all simulation runs
# dev.new()
# par(mfrow=c(2,1),fig=c(0.02,0.98,0.2,0.8))
# barplot(rsummary_pcm/stable_count_pcm,cex.names=0.7,las=2,legend=T,xpd=T,args.legend=list(horiz=T,y=1.25),ylab="Proportion")
# title("Eradication project predictions",line=3)
# 
# # also show only those in which the target species were actually suppressed
# dev.new()
# par(mfrow=c(2,1),fig=c(0.02,0.98,0.2,0.8))
# barplot(rsummary_pcm_subset/unique(colSums(rsummary_pcm_subset)),cex.names=0.7,las=2,legend=T,xpd=T,args.legend=list(horiz=T,y=1.25),ylab="Proportion")
# title("Eradication project predictions (subset)",line=3)



# ################################################################################
# Example 2
# # In this example, we take "uncertain" links into account by iterating over all possible combinations of uncertain links.
# 
# #------------------------------------------------------------------------------
# # read in the interactions data
# # adjust the path here to match where you put the data file
# 
# x=read.table('macquarie_interactions.csv',sep=',',header=T)
# 
# 
# #------------------------------------------------------------------------------
# # initialise variables and settings
# 
# nwrand=1000 # number of randomisations with different random weights on interaction strengths
# max_nwrand=20*nwrand #try a maximum of this many times to get the nwrand results
# 
# # for some variables, we can potentially store them as integers in order to save memory
# # but we need to be careful: R uses 32-bit integers, so for very large nwrand, this might not be suitable
# use_int_storage=.Machine$integer.max > nwrand
# 
# # validation data: what can we use to ground-truth our models?
# # these are responses to suppression of cats and myxoma (i.e. rabbits increased, tall tussock decreased, and cats decreased)
# model_validation=matrix(c('Rabbits',1,'Tall tussock vegetation',-1,'Cats',-1),nrow=3,ncol=2,byrow=TRUE)
# colnames(model_validation)=c('response_node','response_value')
# 
# #> model_validation
# #     response_node             response_value
# #[1,] "Rabbits"                 "1"           
# #[2,] "Tall tussock vegetation" "-1"          
# #[3,] "Cats"                    "-1"          
# 
# 
# # get the list of unique names within these interactions
# node_names=unique(union(x$To,x$From))
# 
# # find the set of 'unknown' links to iterate over
# unk_idx=grep('unknown',x$Importance,ignore.case=T)
# 
# # for the purposes of demonstration, we will only include the first 5 unknown links
# #  comment this line out to include all unknown links (but be prepared: the code will take a long time to run)
# unk_idx=unk_idx[1:5]
# 
# nruns=2^length(unk_idx)
# 
# 
# # initialise some variables that will keep track of various things as we go
# #  suffixes of _cm indicate that this variable applies to the cat/myxo suppression simulations
# #  suffixes of _pcm indicate that this variable applied to the eradication project simulations (i.e. post cats and myxo)
# 
# stable_count_cm=vector("numeric",nruns) #per configuration, number of realisations that were stable
# valid_count_cm=vector("numeric",nruns) #per configuration, number of realisations that were valid
# 
# if (use_int_storage) {
#   rslog_cm=array(as.integer(0),dim=c(nruns,length(node_names),3)) 
# } else {
#   rslog_cm=array(0,dim=c(nruns,length(node_names),3))
# }
# 
# # summary of predictions for each element in the model
# rsummary_cm=matrix(0,nrow=3,ncol=length(node_names))
# rownames(rsummary_cm)=c('Negative','Zero','Positive')
# colnames(rsummary_cm)=node_names
# 
# cm_idx=grep('(cats|myxoma)',node_names,ignore.case=T)
# pcm_idx=grep('[^(cats|myxoma)]',node_names,ignore.case=T) #which elements in the model are relevant for the eradication project (not cats or myxo)
# 
# pcm_node_names=node_names[pcm_idx]
# 
# stable_count_pcm=vector("numeric",nruns) # per configuration, number of realisations that passed the additional stability test (after the removal of cats and myxo from the model)
# 
# pcm_results_cols=grep('(rabbits|rats|mice)',node_names,ignore.case=T) # the columns in the A matrix that give the responses to perturbations of our eradication target species
# 
# # summary of predictions for each element in the model
# rsummary_pcm=matrix(0,nrow=3,ncol=length(pcm_node_names))
# rownames(rsummary_pcm)=c('Negative','Zero','Positive')
# colnames(rsummary_pcm)=pcm_node_names
# 
# # summary of predictions for each element in the model, but only for the subset of simulations in which the target species were actually suppressed
# rsummary_pcm_subset=matrix(0,nrow=3,ncol=length(pcm_node_names)) #summary for successful runs
# rownames(rsummary_pcm_subset)=c('Negative','Zero','Positive')
# colnames(rsummary_pcm_subset)=pcm_node_names
# 
# # for each configuration (run), and for each node, the count of outcome types (either plus, zero, or minus)
# if (use_int_storage) {
#   rslog_pcm=array(as.integer(0),dim=c(nruns,length(pcm_node_names),3)) 
# } else {
#   rslog_pcm=array(0,dim=c(nruns,length(pcm_node_names),3)) 
# }
# 
# # and the same for the subset of eradication runs
# if (use_int_storage) {
#   rslog_pcm_subset=array(as.integer(0),dim=c(nruns,length(pcm_node_names),3)) 
# } else {
#   rslog_pcm_subset=array(0,dim=c(nruns,length(pcm_node_names),3)) 
# }
# 
# # keep count of how many times each possible outcome occurred with eradication
# #  outcomes here are in terms of suppression of target species, ranging from none actually suppressed
# #  to all successfully suppressed (i.e. 8 possible combinations for 3 target species)
# if (use_int_storage) {
#   rslog_pcmall3_outcome_count=matrix(as.integer(0),nrow=nruns,ncol=8)
# } else {
#   rslog_pcmall3_outcome_count=matrix(0,nrow=nruns,ncol=8)
# }
# 
# diagval_max=-0.25 # maximum allowable value for self-limitation links
# 
# n=length(node_names) #number of nodes
# n_pcm=length(pcm_node_names) #number of nodes once cats and myxo removed
# 
# # print out a summary of model
# cat(sprintf('Summary of all non-minor interactions in the model:\n'))
# for (k in (1:dim(x)[1])) {
#   if (grepl('minor',x[k,]$Importance,ignore.case=T)) {
#     next
#   }
#   this_importance = as.character(x[k,]$Importance)
#   if (grepl('unknown',this_importance,ignore.case=T)) {
#     if (! k %in% unk_idx) {
#       this_importance = sprintf('%s ... but ignoring for demo purposes',this_importance)
#     }
#   }   
#   cat(sprintf('Link: %s to %s (%s, %s)\n',x[k,]$From, x[k,]$To, x[k,]$Type, this_importance))
# }
# 
# cat(sprintf('\n'))
# cat(sprintf('\n%d unknown links (%d model runs)\n',length(unk_idx),nruns))
# flush.console()
# 
# #------------------------------------------------------------------------------
# # some supporting functions
# 
# integer.base.b <-
#   # adapted from code posted by Spencer Graves to the R mailing list
#   function(x, b=2){
#     xi <- as.integer(x)
#     if(any(is.na(xi) | ((x-xi)!=0)))
#       print(list(ERROR="x not integer", x=x))
#     N <- length(x)
#     xMax <- max(x)
#     ndigits <- max(1,(floor(logb(xMax, base=2))+1))
#     Base.b <- array(NA, dim=c(N, ndigits))
#     for(i in 1:ndigits){
#       Base.b[, ndigits-i+1] <- (x %% b)
#       x <- (x %/% b)
#     }
#     Base.b
#   }
# 
# bitget=function(x,bit) {
#   # returns the value of bit in x (in binary representation)
#   if (bit<1) stop("positive bit required")
#   x=integer.base.b(x,2)
#   x=t(apply(x,1,rev))
#   if (bit>dim(x)[2]) 0 else x[,bit]
# }
# 
# 
# #------------------------------------------------------------------------------
# # do the actual computations
# 
# cat(sprintf('Main loop started at: %s\n  Progress: ',date()))
# flush.console()
# 
# for (ri in (1:nruns)) {
#   # show progress ...
#   temp=floor(ri/nruns*100)
#   if (temp>floor((ri-1)/nruns*100)) {
#     if (temp>1) {
#       cat(sprintf('\b\b\b\b'))
#     }
#     cat(sprintf(' %d%%',temp))
#     flush.console()
#   }
#   
#   # which links to we want included in this iteration?
#   for (k in (1:length(unk_idx))) {
#     if (bitget(ri-1,k)) {
#       #force this unknown link to be included in this run
#       x[unk_idx[k],]$Importance='Major'
#     } else {
#       x[unk_idx[k],]$Importance='Minor'; # force it to be excluded from this run
#     }
#   }
#   
#   # initialise the interactions matrix A   
#   A=matrix(0,nrow=n,ncol=n)
#   colnames(A)=node_names
#   rownames(A)=node_names
#   
#   # computations expect convention A[i,j] is the effect on i from j
#   
#   for (k in (1:dim(x)[1])) {
#     if (grepl('minor',x[k,]$Importance,ignore.case=T)) {
#       # ignore 'minor' links
#       next
#     }
#     this_from=as.character(x[k,]$From) # the "from" element of the interaction
#     this_to=as.character(x[k,]$To) # the "to" element of the interaction
#     if (grepl('competition',x[k,]$Type,ignore.case=T)) {
#       A[this_from,this_to]=-1
#       A[this_to,this_from]=-1
#     } else if (grepl('habitat',x[k,]$Type,ignore.case=T) || grepl('positive',x[k,]$Type,ignore.case=T)) {
#       A[this_to,this_from]=1
#     } else if (grepl('limiting',x[k,]$Type,ignore.case=T)) {
#       A[this_to,this_from]=-1
#     } else if (grepl('negative',x[k,]$Type,ignore.case=T)) {
#       A[this_to,this_from]=-1
#     } else if (grepl('predator-prey',x[k,]$Type,ignore.case=T)) {
#       A[this_from,this_to]=1
#       A[this_to,this_from]=-1
#     } else if (grepl('scavenging',x[k,]$Type,ignore.case=T)) {
#       A[this_from,this_to]=1
#     } else {
#       stop(sprintf('unrecognised link type %s (%s to %s)',x[k,]$Type,this_from,this_to))
#     }
#   }
#   
#   this_valid_count=0
#   
#   # now loop many times, each time generating new (random) weights in A
#   for (twi in 1:max_nwrand) {
#     # randomisation of weights in A, in range 0.01-1
#     wA=matrix(runif(n^2),nrow=n)*0.99+0.01
#     wA=wA*sign(A)
#     # but make sure diagonals are between -1 and diagval_max
#     temp=diag(rep(1,dim(A)[1]))
#     wA=wA*(1-temp)+(-1-diagval_max)*diag(runif(n))+diagval_max*temp
#     
#     # calculate the response, and check the validation criteria
#     adjwA=-solve(wA) #negative inverse of wA
#     adjwA[abs(adjwA)<1e-07]=0 # set very small responses to zero
#     
#     # check that response fits validation data
#     this_valid=1
#     # check all validation criteria: if any fail, then this model is not valid
#     for (k in 1:dim(model_validation)[1]) {
#       temp_response=-adjwA[,cm_idx]
#       if (dim(temp_response)[2]>1) {
#         # sum response over the individual responses
#         temp_response=rowSums(temp_response)
#       }
#       temp_response=sign(temp_response)
#       if (temp_response[model_validation[k,'response_node']] != model_validation[k,'response_value']) {
#         this_valid=0 #not valid
#         break
#       }
#     }
#     
#     if (!this_valid) next # not valid, so discard this realisation
#     
#     # ok, this realisation is valid: is it also stable?
#     if (!all(Re(eigen(wA,only.values=T)$values)<0)) {
#       # not stable
#       next
#     }  
#     
#     # so if we got this far, this realisation is both valid and stable
#     stable_count_cm[ri]=stable_count_cm[ri]+1
#     this_valid_count=this_valid_count+1
#     
#     # find the predicted response to cat/myxo suppression
#     temp=-adjwA[,cm_idx]
#     if (dim(temp)[2]>1)   temp=rowSums(temp)    # sum response over the individual responses (i.e. total response is the sum of responses to cat and myxo suppression)
#     
#     temp=sign(temp) # only interested in signs of responses
#     
#     rslog_cm[ri,,1]=rslog_cm[ri,,1]+as.integer(temp==-1) 
#     #note that if rslog_cm is not integer (i.e. use_int_storage is false), then rslog_cm[ri,,1]+as.integer(temp==-1) will be of type double anyway, so we will be OK
#     rslog_cm[ri,,2]=rslog_cm[ri,,2]+as.integer(temp==0)
#     rslog_cm[ri,,3]=rslog_cm[ri,,3]+as.integer(temp==1)
#     
#     rsummary_cm[1,]=rsummary_cm[1,]+as.double(temp==-1)
#     rsummary_cm[2,]=rsummary_cm[2,]+as.double(temp==0)
#     rsummary_cm[3,]=rsummary_cm[3,]+as.double(temp==1)
#     
#     
#     # now simulate the eradication project
#     
#     wA=wA[pcm_idx,pcm_idx] #drop the columns and rows associated with cats and myxo
#     
#     #check stability of model with cats and myxo and their associated links removed
#     this_pcm_stable=1
#     if (!all(Re(eigen(wA,only.values=T)$values)<0)) {
#       this_pcm_stable=0
#     }  
#     
#     if (this_pcm_stable) {
#       stable_count_pcm[ri]=stable_count_pcm[ri]+1
#       
#       adjwA=-solve(wA) #negative inverse of wA
#       adjwA[abs(adjwA)<1e-07]=0 # set very small responses to zero
#       
#       # predicted response to eradication project
#       temp=-adjwA[,pcm_results_cols]
#       if (dim(temp)[2]>1) temp=rowSums(temp)
#       temp=sign(temp)
#       rslog_pcm[ri,,1]=rslog_pcm[ri,,1]+as.integer(temp==-1)
#       rslog_pcm[ri,,2]=rslog_pcm[ri,,2]+as.integer(temp==0)
#       rslog_pcm[ri,,3]=rslog_pcm[ri,,3]+as.integer(temp==1)
#       
#       rsummary_pcm[1,]=rsummary_pcm[1,]+as.double(temp==-1)
#       rsummary_pcm[2,]=rsummary_pcm[2,]+as.double(temp==0)
#       rsummary_pcm[3,]=rsummary_pcm[3,]+as.double(temp==1)
#       
#       
#       # now if result matched expected decreases, keep track separately
#       if (all(temp[pcm_results_cols]<0)) {
#         rslog_pcm_subset[ri,,1]=rslog_pcm_subset[ri,,1]+as.integer(temp==-1)
#         rslog_pcm_subset[ri,,2]=rslog_pcm_subset[ri,,2]+as.integer(temp==0)
#         rslog_pcm_subset[ri,,3]=rslog_pcm_subset[ri,,3]+as.integer(temp==1)
#         
#         rsummary_pcm_subset[1,]=rsummary_pcm_subset[1,]+as.double(temp==-1)
#         rsummary_pcm_subset[2,]=rsummary_pcm_subset[2,]+as.double(temp==0)
#         rsummary_pcm_subset[3,]=rsummary_pcm_subset[3,]+as.double(temp==1)
#       }
#       
#       # also keep track of combinations of outcomes
#       for (tempq in 1:8) {
#         #temp[pcm_results_cols] this is the response of the pcm target species (rabbits, rats, and mice)
#         # for this outcome combination, we want to know if only a certain subset of the target species were actually suppressed
#         if (identical(as.double(temp[pcm_results_cols]==-1),as.double(integer.base.b(tempq-1,2)))) {
#           rslog_pcmall3_outcome_count[ri,tempq]=rslog_pcmall3_outcome_count[ri,tempq]+1;
#         }
#       }
#     } # if this_pcm_stable
#     if (this_valid_count==nwrand) {
#       # we are trying up to max_nwrand times, but if we've achieved our target of nwrand valid configurations, bail out
#       break
#     }
#   } #end twi loop
#   valid_count_cm[ri]=this_valid_count
# } # end looping nruns
# cat(sprintf('\n'))
# flush.console()
# 
# #------------------------------------------------------------------------------
# # show some outputs
# 
# # first do some scaling of our logged variables so that each model configuration is equally weighted in the overall results
# #  (this might not be the case if we were unable to achieve the desired nwrand realisations for some model configs)
# 
# tempsummary_cm=matrix(0,nrow=length(node_names),ncol=3)
# tempsummary_pcm=matrix(0,nrow=length(pcm_node_names),ncol=3)
# tempsummary_pcm_subset=matrix(0,nrow=length(pcm_node_names),ncol=3)
# for (k in (1:nruns)) {
#   tempsummary_cm=tempsummary_cm+rslog_cm[k,,]/rowSums(rslog_cm[k,,])[1]
#   tempsummary_pcm=tempsummary_pcm+rslog_pcm[k,,]/rowSums(rslog_pcm[k,,])[1]
#   tempsummary_pcm_subset=tempsummary_pcm_subset+rslog_pcm_subset[k,,]/rowSums(rslog_pcm_subset[k,,])[1]
# }
# rsummary_cm=t(tempsummary_cm/nruns)
# rownames(rsummary_cm)=c('Negative','Zero','Positive')
# colnames(rsummary_cm)=node_names
# 
# rsummary_pcm=t(tempsummary_pcm/nruns)
# rownames(rsummary_pcm)=c('Negative','Zero','Positive')
# colnames(rsummary_pcm)=pcm_node_names
# rsummary_pcm_subset=t(tempsummary_pcm_subset/nruns)
# rownames(rsummary_pcm_subset)=c('Negative','Zero','Positive')
# colnames(rsummary_pcm_subset)=pcm_node_names
# 
# 
# # outcomes of cat/myxo suppression
# dev.new()
# par(mfrow=c(2,1),fig=c(0.02,0.98,0.2,0.8))
# barplot(rsummary_cm,cex.names=0.7,las=2,legend=T,xpd=T,args.legend=list(horiz=T,y=1.25),ylab="Proportion")
# title("Cat/myxoma suppression",line=3)
# 
# 
# # outcomes of eradication project
# # all simulation runs
# dev.new()
# par(mfrow=c(2,1),fig=c(0.02,0.98,0.2,0.8))
# barplot(rsummary_pcm,cex.names=0.7,las=2,legend=T,xpd=T,args.legend=list(horiz=T,y=1.25),ylab="Proportion")
# title("Eradication project predictions",line=3)
# 
# # also show only those in which the target species were actually suppressed
# dev.new()
# par(mfrow=c(2,1),fig=c(0.02,0.98,0.2,0.8))
# barplot(rsummary_pcm_subset,cex.names=0.7,las=2,legend=T,xpd=T,args.legend=list(horiz=T,y=1.25),ylab="Proportion")
# title("Eradication project predictions (subset)",line=3)
# 
# 
# 
# # what outcomes did we get, and how often?
# #  weight results equally across model configurations
# tempoutcome=matrix(0,nrow=nruns,ncol=8)
# for (k in (1:nruns)) {
#   tempoutcome[k,]=rslog_pcmall3_outcome_count[k,]/sum(rslog_pcmall3_outcome_count[k,])
# }
# tempoutcome=as.matrix(colMeans(tempoutcome)*100)
# 
# # construct vector of names of the species that were suppressed in each outcome scenario
# tempfun = function(q) paste(pcm_node_names[pcm_results_cols[which(rev(integer.base.b(q-1)==1))]],collapse="/")
# tempnames=sapply(1:8,tempfun)
# tempnames[1]="None"
# 
# rownames(tempoutcome)=tempnames
# colnames(tempoutcome)="Percentage of models"
# # print results
# cat("Frequency of outcomes: which species were successfully suppressed, how often?")
# print(tempoutcome)
# 


