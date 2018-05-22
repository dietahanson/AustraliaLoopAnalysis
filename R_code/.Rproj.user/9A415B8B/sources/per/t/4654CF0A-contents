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
###Make an nw1 input file for n_w
################################################

nwfile <- data.frame(hotnet$prey, hotnet$predator)
nwfile$hotnet.prey <- gsub(" ", "_", nwfile$hotnet.prey, fixed = TRUE)
nwfile$hotnet.predator <- gsub(" ", "_", nwfile$hotnet.predator, fixed = TRUE)


nnodes <- length(unique(union(nwfile$hotnet.prey,
                              nwfile$hotnet.predator)))

cat(nnodes,"\n", file="martu.nw1")
write.table(nwfile,
            file = "martu.nw1",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F, 
            append = T)





