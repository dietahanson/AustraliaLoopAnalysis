###############################################################################
# Loop model, with the model predicting the effect of
# removing humans. 
###############################################################################



#------------------------------------------------------------------------------
# read in the interactions table
#------------------------------------------------------------------------------

# adjust the path here to match where you put the data file
setwd("~/Documents/Australia/R_code/loop_code/old_code/")
x = read.table('hotgrouped.csv',
             sep=',',
             header = T)


#------------------------------------------------------------------------------
# validation data: what can we use to ground-truth our models?
#------------------------------------------------------------------------------

# these are responses to removal of humans 
model_validation = matrix(c('human',-1,'dingo',1,'skinks',-1),
                        nrow = 3,
                        ncol = 2,
                        byrow = T)

colnames(model_validation) = c('response_node','response_value')

model_validation
# response_node response_value
# [1,] "human"       "-1"           
# [2,] "dingo"         "1"           
# [3,] "skinks"      "1"         


#------------------------------------------------------------------------------
# initialise some variables and settings
#------------------------------------------------------------------------------

nwrand = 100  # number of randomisations with different interaction strengths
max_nwrand = 20*nwrand  # try a maximum of this many realisations 

node_names = unique(union(x$To,x$From))  # get the list of nodes
n = length(node_names)  # number of nodes

this_valid_count = 0  # number of realisations that were valid
stable_count_h = 0  # number of realisations that were stable

h_idx = grep('(human)',node_names,ignore.case=T)  # just humans

# summary of predictions for each element in the model
rsummary_h = matrix(0,nrow = 3,ncol = length(node_names))
rownames(rsummary_h) = c('Negative','Zero','Positive')
colnames(rsummary_h) = node_names


# function to make community matrix symmetric
makeSymm = function(m) {
  m[upper.tri(m)] = t(m)[upper.tri(m)]
  return(m)
}



#------------------------------------------------------------------------------
# Set up the trophic matrices
#------------------------------------------------------------------------------


cat(sprintf('Main loop started at: %s\n',date()))
flush.console()

# initialise the community matrix A (signs) for trophic interactions   
Atro = matrix(0,nrow = n,ncol = n)
colnames(Atro) = node_names
rownames(Atro) = node_names

# initialise the community matrix B (strengths) for trophic interactions   
Btro = matrix(0,nrow = n,ncol = n)
colnames(Btro) = node_names
rownames(Btro) = node_names


# computations expect convention A[i,j] is the effect on i row from j column
for (k in (1:dim(x)[1])) {  # loop through all links in table x
  
  this_from = as.character(x[k,]$From)  # the "from" (predator) element 
  this_to = as.character(x[k,]$To)  # the "to" (prey) element
  
  if (!grepl('predator-prey',x[k,]$Type,ignore.case = T))  
    next  # skip if not a trophic interaction
    Atro[this_from,this_to] = 0.1  # predator conversion efficiency
    Atro[this_to,this_from] = -1
  }



#------------------------------------------------------------------------------
# Set up the non-trophic matrices
#------------------------------------------------------------------------------


# initialise the community matrix A (signs) for non-trophic interactions
Ant = matrix(0,nrow = n,ncol = n)
colnames(Ant) = node_names
rownames(Ant) = node_names

# initialise the community matrix B (strengths) for non-trophic interactions
Bnt = matrix(0,nrow = n,ncol = n)
colnames(Bnt) = node_names
rownames(Bnt) = node_names


# computations expect convention A[i,j] is the effect on i row from j column
for (z in (1:dim(x)[1])) {  

  # if (grepl("(unknown|minor)", x[m, ]$Importance, ignore.case=TRUE) ) {
  #   next  # ignore "minor" or "unknown" links
  # }
  
  if (grepl("Predator-prey", x[z,]$Type, ignore.case = T)) {
    next  # ignore trophic interactions
  }
  
  this_from = as.character(x[z, ]$From) ## the "from" (predator) element 
  this_to = as.character(x[z, ]$To) ## the "to" (prey) element 
  
  if (grepl("competition", x[z, ]$Type, ignore.case=T)) {
    Ant[this_from, this_to] = -1
    Ant[this_to, this_from] = -1
  } else if (grepl("habitat", x[z, ]$Type, ignore.case=T) ||
             grepl("positive", x[z, ]$Type, ignore.case=T)) {
    Ant[this_to, this_from] = 1
  } else if (grepl("limiting", x[z, ]$Type, ignore.case=T)) {
    Ant[this_to, this_from] = -1
  } else if (grepl("negative", x[z, ]$Type, ignore.case=T)) {
    Ant[this_to, this_from] = -1
  } else if (grepl("scavenging", x[z, ]$Type, ignore.case=T)) {
    Ant[this_from, this_to] = 1
  } else {
    stop(sprintf("unrecognised link type %s (%s to %s)",
                 x[z, ]$Type, this_from, this_to))
  }
}



#------------------------------------------------------------------------------
# now loop max_nwrand, each time generating new (random) weights in B matrices
#------------------------------------------------------------------------------


for (twi in 1:max_nwrand) {
  
  
  # populate trophic B matrix with beta-dist values between .01 and 1
  Btro = makeSymm(matrix(rbeta(n^2,
                                shape1 = 1,
                                shape2 = 4),
                          nrow = n))*.99+0.01  
  colnames(Btro) = node_names
  rownames(Btro) = node_names
  

  Ctro=Btro*Atro  # add signs to strengths for the trophic matrix

  # populate diagonal elements of trophic matrix 
  for (m in (1:length(node_names))) {
    
    sp = node_names[m]  

    if (sp %in% x$From) {  # if the species is a predator (in the "From" col)
      Ctro[sp, sp] = -.1 
    } else { Ctro[sp, sp] = -1}  # if the species is basal
  }
  
  
  
  # populate non-trophic B matrix with uniform values between .01 and 1
  Bnt = matrix(runif(n^2), nrow=n)*0.99+0.01 
  colnames(Bnt) = node_names
  rownames(Bnt) = node_names
  
  Cnt=Bnt*Ant  # add signs to strengths for the non-trophic matrix
  diag(Cnt) = 0  # diagonal should be zero for non-trophic interactions
  
  
  # sum the trophic and non-trophic matrices to get the total matrix
  Ctot = Ctro + Cnt

  
  
  # calculate the response, and check the validation criteria
  adjwA = -solve(Ctot)  # negative inverse of total matrix
  adjwA[abs(adjwA)<1e-07] = 0  # set very small responses to zero
  
  
  # check that response fits validation data
  this_valid = 1
  
  # check all validation criteria: if any fail, then this model is not valid
  for (k in 1:dim(model_validation)[1]) {
    
    temp_response = -adjwA[,h_idx]  # look at human column
    temp_response = sign(temp_response)  # just take sign
    
    if (temp_response[model_validation[k,'response_node']] 
        != model_validation[k,'response_value']) {
      this_valid = 0  # not valid
      break
    }
  }

  if (!this_valid) next  # not valid, so discard this realisation

  this_valid_count = this_valid_count+1  # valid, so add to count
  
  
  # this realisation is valid: is it also stable (eigenvalues are negative/real)
  if (!all(Re(eigen(Ctot,only.values = T)$values)<0)) {  # not stable
    next
  }  
  
  stable_count_h = stable_count_h+1  # stable, so add to count

  
  # find the predicted response to human decline
  temp = -adjwA[,h_idx] 
  temp = sign(temp) # only interested in signs of responses
  
  # keep tabs on summary results
  rsummary_h[1,] = rsummary_h[1,]+as.double(temp == -1)
  rsummary_h[2,] = rsummary_h[2,]+as.double(temp == 0)
  rsummary_h[3,] = rsummary_h[3,]+as.double(temp == 1)

  if (stable_count_h == nwrand) {
    # try up to max_nwrand times, but stop when nwrand achieved
    break
  }
}  # end twi loop

cat(sprintf('\n'))
flush.console()


#------------------------------------------------------------------------------
# plot the results
#------------------------------------------------------------------------------

dddd = gsub("-", "", as.character(Sys.Date()))

v = paste("Number of valid models:", this_valid_count,
          "Number of stable models:", stable_count_h,
          "Number of models tried:", twi,
          "\nConsumer Limitation:", max(diag(Ctot)),
          "Basal Limitation:", min(diag(Ctot)))

# pdf(paste("australia_loop_outcome",
#           dddd,
#           # max(diag(Ctot)),
#           ".pdf",
#           sep = ''),
#     width = 11, height = 8.5)
# 
# par(mfrow = c(2,1),
#     fig = c(0.02,0.98,0.2,0.8),
#     mar = c(4,4,4,7),
#     xpd = T)
# 
# barplot(rsummary_h/stable_count_h,
#         cex.names = 0.7,
#         las = 2,
#         legend = T,
#         xpd = T,
#         args.legend = list(x = 'topright', inset = c(-.1,0)),
#         ylab = "Proportion")
# 
# title(main = v, cex.main= .9, line = 2)

print(v)


barplot(rsummary_h/stable_count_h,
        cex.names = 0.7,
        las = 2,
        legend = T,
        xpd = T,
        args.legend = list(x = 'topright', inset = c(-.1,0), cex = 0.7),
        ylab = "Proportion")

title(main = v, cex.main = .9, line = 2)

