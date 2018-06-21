##functions that are needed as input into bayestest_dieta

library("data.table")

#-------------------------------------------------------------------------------
# functions
#-------------------------------------------------------------------------------

# function to make a matrix symmetrical (used for trophic matrices)
makeSymm = function(m) {
  m[upper.tri(m)] = t(m)[upper.tri(m)]
  return(m)
}


# function to populate a trophic matrix with interaction strengths
trophic = function(n, node_names, A, x) {  
  
  Bt = makeSymm(matrix(rbeta(n^2,
                             shape1 = 1,
                             shape2 = 4),
                       nrow = n))*.99+0.01  
  colnames(Bt) = node_names
  rownames(Bt) = node_names
  
  
  Ct=Bt*A  # add signs to strengths
  
  # populate diagonal elements 
  for (m in (1:length(node_names))) {
    
    sp = node_names[m]  
    
    if (sp %in% x[x$Type=="Predatorprey",]$From) {   # if a predator 
      Ct[sp, sp] = -.1 
    } else { Ct[sp, sp] = -1}  # if the species is basal
  }
  return(Ct)
}


# function to populate a non-trophic matrix with interaction strengths
nontrophic = function(n, node_names, A, x) {  
  
  Bnt = matrix(runif(n^2), nrow=n)*0.99+0.01
  colnames(Bnt) = node_names
  rownames(Bnt) = node_names
  
  
  Cnt=Bnt*A  # add signs to strengths
  
  # populate diagonal elements of trophic matrix 
  diag(Cnt) = 0  # diagonal should be zero for non-trophic interactions
  return(Cnt)
}


# function to make a unique identifier for each species pair in edges df
pairnamer = function (df) {return(paste(unlist(sort(df)), collapse = ""))
}


##function to generate random interaction strengths

sampler = function(need, n, node_names, am, x) {
  
  sm = list()    
  
  for (g in 1:length(need)) {
    
    if (grepl("predatorprey", need[g], ignore.case = T)) {
      tm = trophic(n, node_names, am[[need[g]]], x)
      sm[[g]] = tm
    } else {
      ntm = nontrophic(n, node_names, am[[need[g]]])
      sm[[g]] = ntm
    }
  } 
  
  popmod = Reduce("+", sm)  # add all sub matrices
  
  edges = data.frame(as.table(popmod), stringsAsFactors = F)
  edges$Pair = apply(edges[,c("Var1", "Var2")],1,pairnamer) 
  edges = edges[abs(edges$Freq)>0,]
  setnames(edges,
           old = c("Var1", "Var2", "Freq", "Pair"),
           new = c("From", "To", "Type", "Pair"))
  edges$Group = 0
  edges$Type = sapply(edges$Type,
                      function(x) if (x < 0) { x = "N"} else {x = "P"})
  
  res = list(cmat = popmod, edges = edges)
  
}



#------------------------------------------------------------------------------
# get data and initialize settings
#------------------------------------------------------------------------------


setwd("~/Documents/Australia/R_code/loop_code/")  # adjust as needed

x = read.table("hotgrouped_fake.csv",  # interactions table 
               sep = ",",
               header = T,
               stringsAsFactors = F)

y = read.csv("models.csv", header = T, stringsAsFactors = F)


node_names = unique(union(x$To,x$From))  # get the list of nodes
n = length(node_names)  # number of nodes

t = unique(x$Type)  # how many different interaction types there are






#------------------------------------------------------------------------------
# loop through all possible interaction types and make sign matrices for 
# each. Output a list "am". This is too long! probably some better way....
#------------------------------------------------------------------------------


am = list()  # list of sign matrices

for (f in 1:length(t)) {
  
  
  # initialise the sign matrix A   
  A = matrix(0,nrow = n,ncol = n)
  colnames(A) = node_names
  rownames(A) = node_names
  
  if (grepl("predatorprey", t[f], ignore.case = T)) {
    
    for (k in (1:dim(x)[1])) {  # loop through all links in table x
      
      this_from = as.character(x[k,]$From)  
      this_to = as.character(x[k,]$To)  
      
      if (!grepl("predatorprey",x[k,]$Type,ignore.case = T))  
        next  # skip if not a trophic interaction
      A[this_from,this_to] = 0.1  # predator conversion efficiency
      A[this_to,this_from] = -1
      
      am[[f]] = A
      names(am)[[f]] = t[f]
    }
    
  } else if (grepl("habitat", t[f], ignore.case = T)) {
    
    for (k in (1:dim(x)[1])) {  # loop through all links in table x
      
      this_from = as.character(x[k,]$From) 
      this_to = as.character(x[k,]$To)  
      
      if (!grepl("habitat",x[k,]$Type,ignore.case = T))  
        next  # skip if not a habitat interaction
      A[this_to,this_from] = 1
      
      am[[f]] = A
      names(am)[[f]] = t[f]
    }
    
  } else if (grepl("positive", t[f], ignore.case = T)) {
    
    for (k in (1:dim(x)[1])) {  # loop through all links in table x
      
      this_from = as.character(x[k,]$From) 
      this_to = as.character(x[k,]$To)
      
      if (!grepl("positive",x[k,]$Type,ignore.case = T))  
        next  # skip if not a positive interaction
      A[this_to,this_from] = 1
      
      am[[f]] = A
      names(am)[[f]] = t[f]
    }
    
  } else if (grepl("competition", t[f], ignore.case = T)) {
    
    for (k in (1:dim(x)[1])) {  # loop through all links in table x
      
      this_from = as.character(x[k,]$From)   
      this_to = as.character(x[k,]$To) 
      
      if (!grepl("competition",x[k,]$Type,ignore.case = T))  
        next  # skip if not a competitive interaction
      A[this_from, this_to] = -1
      A[this_to, this_from] = -1
      
      am[[f]] = A
      names(am)[[f]] = t[f]
    }
  }
}



  
  
  
  
  
  
  