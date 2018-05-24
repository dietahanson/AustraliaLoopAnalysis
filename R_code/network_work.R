###############################################
###Add groups to network
###############################################


stefgroup <- read.csv("stefgroups.csv",  #get stefanis groups 
                      header = T,
                      stringsAsFactors = F) 

tg <- read.delim("tg_martu.txt", #get trophic groups from n_w
                 header = T, 
                 stringsAsFactors = F,
                 sep = "\t")

tg$Species <- gsub("_", " ", tg$Species)  #fix species names
tg$Species <- gsub("human", "Homo sapiens", tg$Species)
tg$Species <- gsub("\\*", "", tg$Species)

aic <- read.delim("aic_martu.txt", #get aic groups from n_w
                 header = T, 
                 stringsAsFactors = F,
                 sep = "\t")

aic$Species <- gsub("_", " ", aic$Species) #fix species names
aic$Species <- gsub("human", "Homo sapiens", aic$Species)
aic$Species <- gsub("\\*", "", aic$Species)

groups <- data.frame(common = stefgroup$common, #make table of species names
                     latin = stefgroup$latin,
                     stefgroup = stefgroup$stefgroup)

groups$tg <- tg$X.TroG[match(groups$latin, tg$Species)] #add trophic groups

groups$aic <- aic$X.AicG[match(groups$latin, aic$Species)] #add aic groups

groups[is.na(groups$tg),] #returns species with no group match in tg
groups[is.na(groups$aic),] #returns species with no group match in aic




###################################################################
#####Make input for Loop Analysis
###################################################################



hotnet <- read.csv("hotnet.csv",  #get network
                   header = T,
                   stringsAsFactors = F) 

hotnet$weight <- NULL #take away weights



##Change which groups you want to use (stefani's, tg, or aic) by changing the
##groups column that is selected

hotnet$From <- as.factor(groups$tg[match(hotnet$predator, #match predators
                                               groups$latin)])

cat(sum(is.na(hotnet$From)), #print how many did not match
    "predators with missing matches:")
hotnet[is.na(hotnet$From),]

hotnet$To <- as.factor(groups$tg[match(hotnet$prey,   #match prey
                                           groups$latin)])

cat(sum(is.na(hotnet$To)), #print how many did not match
    "prey with missing matches:")
    hotnet[is.na(hotnet$To),]

hotnet$Type <- "Predator-prey" #add interaction type

#write.csv(hotnet, file = "hotgrouped.csv", row.names = F)





