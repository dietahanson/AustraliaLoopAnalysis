## Code to produce loop models and to choose the best model.
## The only three required inputs are "networktable": a table of all 
## interactions, with 3 columns; from = first species in interaction, 
## to = second species in interaction, and type = category of interaction 
## (predatorprey, habitat, etc), "models": a table of models, with each row 
## being a model and each column being an interaction type to include in the 
## model. Interaction types are included with a 1 and excluded with a 0, and 
## "outcomes": a table of the known outcomes which will be used to validate the 
## models. The first column contains the species for which you have known 
## outcomes, and the second column should contain 1 or -1 to specify the 
## direction of the known outcome. The first entry of this table must be the 
## node which is being pressed. 



##------------------------------------------------------------------------------
## All required functions and packages
##------------------------------------------------------------------------------


library(data.table)
library(svMisc)

# Function to make a matrix symmetrical (used for trophic matrices)
makeSymm = function(m) {
  m[upper.tri(m)] = t(m)[upper.tri(m)]
  return(m)
}



# Function to populate a trophic matrix with interaction strengths. Requires 
# inputs n = number of nodes, node_names, A = sign matrix, x = master table of 
# interactions
trophic = function(n, node_names, A, x) {    
  
  # first fill the matrix with random beta values
  Bt = makeSymm(matrix(rbeta(n^2,  # make symmetric so ij is a function of ji
                             shape1 = 1, shape2 = 4),  # parameters of beta dist
                       nrow = n))*.99 + 0.01 # make all values between 0.1 and 1  
  colnames(Bt) = node_names
  rownames(Bt) = node_names
  
  # then add signs to strengths
  Ct=Bt*A  
  
  # last, populate diagonal elements 
  for (m in (1:length(node_names))) {
    
    sp = node_names[m] 
    
    # only predators will be in the from column (for predatorprey interactions)
    if (sp %in% x[grepl("predatorprey", x$type, ignore.case=TRUE),]$from) {  
      Ct[sp, sp] = -.1  # if a predator  
    } else { Ct[sp, sp] = -1}  # if basal
  }
  return(Ct)
}



# Function to populate a non-trophic matrix with interaction strengths. Requires 
# inputs n = number of nodes, node_names, A = sign matrix, x = master table of 
# interactions
nontrophic = function(n, node_names, A, x) {  
  
  # first fill the matrix with random uniform values
  Bnt = matrix(runif(n^2), nrow=n)*0.99+0.01  # values between 0.1 and 1
  colnames(Bnt) = node_names
  rownames(Bnt) = node_names
  
  # then add signs to strengths
  Cnt=Bnt*A  
  
  # last, populate diagonal elements
  diag(Cnt) = 0  # diagonal will be populated from the trophic matrix
  return(Cnt)
}



# Function to make a unique identifier for each species pair in edges table
pairnamer = function (df) {return(paste(unlist(sort(df)), collapse = ""))
}



# Function to generate the total community matrix, which is the sum of the 
# various interaction type community matrices required by the model. Requires 
# need = a list of interaction types (predatorprey, habitat, etc.) required by
# the model, n = number of nodes, node_names, am = a list of sign matrices for 
# each interaction type, named by interaction type, and x = master table of 
# interactions
sampler = function(need, n, node_names, am, x) {
  
  sm = list()  # list to contain the interaction type matrices    
  
  # for each interaction type required by the model, populate a matrix with 
  # random values (will use the trophic function if a predatorprey interaction,
  # otherwise will use the nontrophic function) and put into the sm list
  for (g in 1:length(need)) {  
    
    if (grepl("predatorprey", need[g], ignore.case = T)) {
      tm = trophic(n, node_names, am[[need[g]]], x)  # appropriate sign matrix
      sm[[g]] = tm
    } else {
      ntm = nontrophic(n, node_names, am[[need[g]]])
      sm[[g]] = ntm
    }
  } 
  
  # then add all type matrices
  popmod = Reduce("+", sm)  
  
  # now create a table "edges" that has columns; from = first species in 
  # interaction, to = second species in interaction, type = "N" for negative or
  # "P" for positive, and Pair = unique pair name where ij and ji are both 
  # called ij. Note this is not exactly the same format as the master 
  # interactions table x as it will only contain the interactions called by the 
  # particular model, the type is only negative or positive, and Pairs is added.
  # Also, each two-way interaction is expanded to two rows, so one predatorprey 
  # interaction in x will here be one positive interaction from prey to predator
  # and one negative interaction for predator to prey
  edges = data.frame(as.table(popmod), stringsAsFactors = F)
  edges$Pair = apply(edges[,c("Var1", "Var2")],1,pairnamer) 
  edges = edges[abs(edges$Freq)>0,]  # we only need non-zero interactions
  setnames(edges,
           old = c("Var1", "Var2", "Freq", "Pair"),
           new = c("from", "to", "type", "Pair"))
  for (i in 1:nrow(edges)) {
    if (edges$from[i]==edges$to[i]) {edges$type[i]="N"} else {
      if (edges$type[i] < 0) {edges$type[i] = "P"} else {
        edges$type[i] = "N"}
      }
    }
  

  
  # create a list containing the community matrix and the edges table
  res = list(cmat = popmod, edges = edges)
  
}


# Function to check the stability of a simulated community matrix W. If all 
# eigenvalues have a negative real part, matrix is stable
stable.community = function(W) {
  all(Re(eigen(W,symmetric = F,only.values = T)$values)<0)
}



# Function to return sign of s, with values of magnitude less than epsilon
# rounded down to zero
signum = function(s,epsilon = 1.0E-5) {
  (s > epsilon) - (s < -epsilon)
}



# Function to generate a function to check a perturbation outcome against some 
# expectation. Must supply a vector of named elements that specify the relative 
# magnitude of the press perturbation, and a vector of named elements that 
# specify the signs of the change in the monitored nodes
press.validate = function(edges,perturb,monitor,epsilon = 1.0E-5) {
  
  # first make a function that calls the name of a node and retrieves its index
  # in the edges "from" levels
  index = function(name) {   
    k = match(name,levels(edges$from))
    if(any(is.na(k)))
      warning("Unknown nodes:",paste(name[is.na(k)],collapse = " "))
    k
  }
  
  # get the indices of perturbed and monitored nodes
  k.perturb = index(names(perturb))
  k.monitor = index(names(monitor))
  
  # make a vector the length of the edges table and make the the pertubed node 
  # element the negative of the perturbation
  S.press = double(nlevels(edges$from))
  S.press[k.perturb] = -perturb
  monitor = sign(monitor) 
  
  # return a function to check the perturbation condition, and a logical 
  # response asking if the outcome matches the expectation of monitor
  function(W) {
    s = tryCatch(solve(W,S.press),error=function(e) NULL)
    !is.null(s) && all(signum(s[k.monitor],epsilon)==monitor)
  }
}


# Function to generate a function to determine the impact of a 
# press perturbation
press.impact = function(edges,perturb,monitor) {
  
  # get the index of a given node ("name") in the levels of the "from" column
  index = function(name) {
    k = match(name,levels(edges$from))
    if(any(is.na(k)))  # give a warning if not found
      warning("Unknown nodes:",paste(name[is.na(k)],collapse=" "))
    k
  }
  
  # get the indices of perturbed nodes and set their position in a vector of 
  # nodes to -1 (or some x greater than -1 if there is more than perturbed 
  # node and the relative magnitude is x greater than some other perturbed node)
  k.perturb = index(names(perturb))
  S.press = double(nlevels(edges$from))
  S.press[k.perturb] = -perturb
  
  if(length(monitor)==0) {
    impact = function(W) solve(W,S.press)  # look at all impacts
  } else {
    k.monitor = index(names(monitor))
    impact = function(W) solve(W,S.press)[k.monitor]  # only impacts on monitor
  }
  
  ## Return function to compute impact
  impact
}


#------------------------------------------------------------------------------
# Get data and set variables
#------------------------------------------------------------------------------

# on local
setwd("~/Documents/Australia/R_code/loop_code/")  # adjust as needed
networktable = "~/Documents/Australia/R_code/loop_code/hotgrouped.csv"
models = "~/Documents/Australia/R_code/loop_code/models.csv"
outcomes = "~/Documents/Australia/R_code/loop_code/outcomes.csv"

# on vulpes
# networktable = "/home/dieta/Australia/loop_code/hotgrouped.csv"
# models = "/home/dieta/Australia/loop_code/models.csv"
# outcomes = "/home/dieta/Australia/loop_code/outcomes.csv"



x = read.csv(networktable, stringsAsFactors = F)  # interactions table
x[] = lapply(x, tolower)  # change everything to lowercase

y = read.csv(models, stringsAsFactors = F)  # models table
y[] = lapply(y, tolower)  # change everything to lowercase

z = read.csv(outcomes, stringsAsFactors = F, header = F)  # outcomes table
z[] = lapply(z, tolower) # change everything to lowercase

# check that there are no nodes or interactions appearing in y or z that are 
# not in x

if (sum(!z$V1 %in% union(x$from, x$to))>0) {
  print(paste("Error:",
              z[!z$V1 %in% union(x$from, x$to),]$V1, 
              "in outcomes table is not found in network table"))}

if (sum(!colnames(y) %in% x$type)>0) {
  print(paste("Error:",
               colnames(y)[!colnames(y) %in% x$type], 
               "in models table is not found in network table"))}

if (sum(!unique(x$type) %in% colnames(y))>0) {
  print(paste("Error:", 
              unique(x$type)[!unique(x$type) %in% colnames(y)], 
              "in network table is not found in models table"))}



node_names = unique(union(x$to,x$from))  # get the list of nodes
n = length(node_names)  # number of nodes

t = unique(x$type)  # how many different interaction types there are


# set how many stable and valid realizations you want to achieve
n.samples = 500
n.max = 4000*n.samples  # try a maximum of this many times (to prevent runaways)

# which model is going to be used? (This should be the model with the highest 
# posterior probability from the model_select script)

m = 1

# set the perturbed node

perturb = setNames(as.numeric(z$V2[1]), z$V1[1])
print(paste((z$V1[1]),
            "is being",
            if (sign(as.numeric(z$V2[1]))<0) {
              "decreased"
              } else {
                "increased"
                }))

# set the monitored nodes, whose sign represents the expected outcome
monitor = setNames(as.numeric(z$V2), z$V1)

# date
dddd <- gsub("-", "", as.character(Sys.Date()))



#------------------------------------------------------------------------------
# Make a list of all sign matrices needed. Output a list "am". 
#------------------------------------------------------------------------------


am = list()  # list of sign matrices

for (f in 1:length(t)) {
  
  
  # initialise the sign matrix A   
  A = matrix(0,nrow = n,ncol = n)
  colnames(A) = node_names
  rownames(A) = node_names
  
  if (grepl("predatorprey", t[f], ignore.case = T)) {
    
    for (k in (1:dim(x)[1])) {  # loop through all links in table x
      
      this_from = as.character(x[k,]$from)  
      this_to = as.character(x[k,]$to)  
      
      if (!grepl("predatorprey",x[k,]$type,ignore.case = T))  
        next  # skip if not a trophic interaction
      A[this_from,this_to] = 0.1  # predator conversion efficiency
      A[this_to,this_from] = -1
      
      am[[f]] = A
      names(am)[[f]] = t[f]
    }
    
  } else if (grepl("simplefire", t[f], ignore.case = T)) {
    
    for (k in (1:dim(x)[1])) {  # loop through all links in table x
      
      this_from = as.character(x[k,]$from) 
      this_to = as.character(x[k,]$to)  
      
      if (!grepl("simplefire",x[k,]$type,ignore.case = T))  
        next  # skip if not a simple fire interaction
      A[this_from,this_to] = 1  
      A[this_to,this_from] = -1
      
      am[[f]] = A
      names(am)[[f]] = t[f]
    }
    
  } else if (grepl("\\bpositive\\b", t[f], ignore.case = T)) {
    
    for (k in (1:dim(x)[1])) {  # loop through all links in table x
      
      this_from = as.character(x[k,]$from) 
      this_to = as.character(x[k,]$to)
      
      if (!grepl("\\bpositive\\b",x[k,]$type,ignore.case = T))  
        next  # skip if not a positive interaction
      A[this_to,this_from] = 1
      
      am[[f]] = A
      names(am)[[f]] = t[f]
    }
    
  } else if (grepl("complexfire", t[f], ignore.case = T)) {
    
    for (k in (1:dim(x)[1])) {  # loop through all links in table x
      
      this_from = as.character(x[k,]$from)   
      this_to = as.character(x[k,]$to) 
      
      if (!grepl("complexfire",x[k,]$type,ignore.case = T))  
        next  # skip if not a complex fire interaction
      A[this_from, this_to] = 1
      A[this_to, this_from] = -1
      
      am[[f]] = A
      names(am)[[f]] = t[f]
    }
    
  } else if (grepl("complexpositive", t[f], ignore.case = T)) {
    
    for (k in (1:dim(x)[1])) {  # loop through all links in table x
      
      this_from = as.character(x[k,]$from)   
      this_to = as.character(x[k,]$to) 
      
      if (!grepl("complexpositive",x[k,]$type,ignore.case = T))  
        next  # skip if not a complex positive interaction
      A[this_to, this_from] = 1
      
      am[[f]] = A
      names(am)[[f]] = t[f]
    }
  } else if (grepl("complexnegative", t[f], ignore.case = T)) {
    
    for (k in (1:dim(x)[1])) {  # loop through all links in table x
      
      this_from = as.character(x[k,]$from)   
      this_to = as.character(x[k,]$to) 
      
      if (!grepl("complexnegative",x[k,]$type,ignore.case = T))  
        next  # skip if not a complex negative interaction
      A[this_to, this_from] = -1
      
      am[[f]] = A
      names(am)[[f]] = t[f]
    }
  }
}


#------------------------------------------------------------------------------
# Run the analysis which randomly selected a model, populates a matrix with
# random interaction strengths, runs a loop analysis, and counts if the 
# realization was stable and valid (matched the expected outcome). Whichever 
# model was stable and valid most often has the highest posterior probability
#------------------------------------------------------------------------------

print(paste("Main loop started at:", date()))
print("Progress")  # just making new lines so progress doesn't cover up above
print("Progress")  # just making new lines so progress doesn't cover up above


## Outcomes
impacts = matrix(0,n.samples,n)
valid = 0
tried = 0
accepted = 0

# summary of predictions for each element in the model
rsummary = matrix(0,nrow = 3,ncol = length(node_names))
rownames(rsummary) = c('Negative','Zero','Positive')
colnames(rsummary) = node_names

while((accepted < n.samples) && (tried < n.max)) {
  
  
  # get a list of the types of interactions the model needs
  need = as.numeric(y[m,])
  names(need) = colnames(y)
  need = names(need[need > 0])
  

  # simulate a matrix W for the model
  H = sampler(need, n, node_names, am, x)
  W = H$cmat
  tried = tried + 1  # count this towards attempted realizations
  
  # generate a function to see if the impact of the press perturbation is the
  # same as the expected outcomes "monitor"
  press = press.validate(H$edges,perturb = perturb, monitor = monitor)
  
  # call the function to check press condition 
  if(!press(W)) next  # skip if not valid
  
  valid = valid + 1
  
  # check stability and exit if not stable
  if(!stable.community(W)) next
  
  accepted = accepted + 1  # count if valid and stable
  
  progress(accepted, max.value = n.samples)  # show progress in the console
  
  # generate a function to determine the outcome of the press perturbation
  impact <- press.impact(H$edges,perturb = perturb,monitor = NULL)
  
  ## Monitor impact post press
  imp <- signum(impact(W))  # get the sign of the impact
  
  # keep tabs on summary results
  rsummary[1,] = rsummary[1,]+as.double(imp == -1)
  rsummary[2,] = rsummary[2,]+as.double(imp == 0)
  rsummary[3,] = rsummary[3,]+as.double(imp == 1)
}

cat(paste("Main loop finished at:", date(),"\n"))
colnames(impacts) <- node_names

################################################################################
## make some outputs
################################################################################

# print some stats for the log file

print(paste("accepted:",accepted))
print(paste("valid:", valid))
print(paste("tried:",tried))

barplot(rsummary/accepted,
        cex.names = 0.7,
        las = 2,
        legend = T,
        xpd = T,
        args.legend = list(x = 'topright', inset = c(-.1,0), cex = 0.7),
        ylab = "Proportion")

title(main = v, cex.main = .9, line = 2)


