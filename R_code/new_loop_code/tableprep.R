###############################################################################
## Code to translate interaction tables into QPress form
###############################################################################


setwd("~/Documents/Australia/R_code/new_loop_code")


it = read.table('hotgroupedtest.csv',
                     sep=',',
                     header = T,
                stringsAsFactors = F)

it$Group = 0

for (k in (1:dim(it)[1])) {  # loop through all links in table it
  
  this_from = as.character(it[k,]$From)  # the "from" (predator) element 
  this_to = as.character(it[k,]$To)  # the "to" (prey) element
  
  if (grepl("predator-prey",it[k,]$Class,ignore.case = T))  {
  
    it$Type[k] = "N"
    it$Pair[k] = k
    it = rbind(it, data.frame(From = this_to,
                              To = this_from,
                              Class = "Predator-prey",
                              Group = 0,
                              Type = "P",
                              Pair = k))
    
  } else if  (grepl("competition", it[k, ]$Class, ignore.case=TRUE)) {
    
    it$Type[k] = "N"
    it$Pair[k] = k
    it = rbind(it, data.frame(From = this_to,
                              To = this_from,
                              Class = "competition",
                              Group = 0,
                              Type = "N",
                              Pair = k))
  
    
  } else if (grepl("habitat", it[k, ]$Class, ignore.case=TRUE) ||
             grepl("positive", it[k, ]$Class, ignore.case=TRUE)) {
    
    it$Type[k] = "P"
    it$Pair[k] = k
    
    if (this_from)
    
  } else {
    stop(sprintf("unrecognised link type %s (%s to %s)", it[k, ]$Class, this_from, this_to))
  }   

}
