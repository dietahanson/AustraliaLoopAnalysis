###############################################
###Put Stefani's groups into her network
###############################################

library(reshape2)


hotnet <- read.csv("hotnet.csv",  #get network
                   header = T,
                   stringsAsFactors = F) 

hotnet$weight <- NULL #take away weights

stefgroup <- read.csv("stefgroups.csv",  #get groups 
                      header = T,
                      stringsAsFactors = F) 

stefgroup<-melt(stefgroup, #melt common names and latin names
                id.vars = "Group",
                value.name = "name")
                                                                          


hotnet$From <- as.factor(stefgroup$Group[match(hotnet$predator, #match predators
                                               stefgroup$name)])

cat(sum(is.na(hotnet$From)), #print how many did not match
    "predators with missing matches:",
    hotnet[is.na(hotnet$From)])

hotnet$To <- as.factor(stefgroup$Group[match(hotnet$prey,   #match prey
                                           stefgroup$name)])

cat(sum(is.na(hotnet$To)), #print how many did not match
    "prey with missing matches:",
    hotnet[is.na(hotnet$To)])

hotnet$Type <- "Predator-prey" #add interaction type

#write.csv(hotnet, file = "hotgrouped.csv", row.names = F)


################################################
###Make a 
################################################



