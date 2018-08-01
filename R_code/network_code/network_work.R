###############################################
###Add groups to network
###############################################
library(splitstackshape)

stefgroup <- read.csv("stefgroups.csv",  # get stefanis groups 
                      header = T,
                      stringsAsFactors = F) 

tg <- read.delim("tg_martu.txt",  # get trophic groups from n_w
                 header = T, 
                 stringsAsFactors = F,
                 sep = "\t")

tg$Species <- gsub("_", " ", tg$Species)   # fix species names
tg$Species <- gsub("human", "Homo sapiens", tg$Species)
tg$Species <- gsub("\\*", "", tg$Species)

aic <- read.delim("aic_martu.txt",  # get aic groups from n_w
                 header = T, 
                 stringsAsFactors = F,
                 sep = "\t")

aic$Species <- gsub("_", " ", aic$Species)  # fix species names
aic$Species <- gsub("human", "Homo sapiens", aic$Species)
aic$Species <- gsub("\\*", "", aic$Species)

agg <- read.delim("agg_martu.txt",  # get aggregated groups from n_w
                  header = T,
                  stringsAsFactors = F, 
                  sep = "\t")


agg$X.Species <- gsub("_", " ", agg$X.Species)  # fix species names
agg$X.Species <- gsub("human", "Homo sapiens", agg$X.Species)
agg$X.Species <- gsub("\\*", "", agg$X.Species)
agg$X.Species <- gsub("sp..", "sp.", agg$X.Species, fixed = T)

agg <- data.frame(species = agg$X.Species,
                  group = agg$X.)

agglong <- as.data.frame(cSplit(agg, "species",
                                ".",
                                direction = "long",
                                stripWhite = F))

agglong$species <- gsub("sp$", "sp.", agglong$species)



groups <- data.frame(common = stefgroup$common,  # make table of species names
                     latin = stefgroup$latin,
                     stefgroup = stefgroup$stefgroup)

groups$tg <- tg$X.TroG[match(groups$latin,
                             tg$Species)]  # add trophic groups

groups$aic <- aic$X.AicG[match(groups$latin,
                               aic$Species)]  # add aic groups

groups$agg <- agglong$group[match(groups$latin,
                                  agglong$species)]  # add agg groups

groups[is.na(groups$tg),]  # returns species with no group match in tg
groups[is.na(groups$aic),]  # returns species with no group match in aic
groups[is.na(groups$agg),]  # returns species with no group match in agg

# write.csv(groups,
#           file = "compare_groups.csv",
#           row.names = F)

###################################################################
#####Make input for Loop Analysis
###################################################################



hotnet <- read.csv("hotnet.csv",   # get network
                   header = T,
                   stringsAsFactors = F)


# Change which groups you want to use (stefani's, tg, or aic) by changing the
# groups column that is selected

hotnet$from <- as.factor(groups$stefgroup[match(hotnet$predator,  # match predators
                                               groups$latin)])

cat(sum(is.na(hotnet$from)),  # print how many did not match
    "predators with missing matches:")
hotnet[is.na(hotnet$from),]

hotnet$to <- as.factor(groups$stefgroup[match(hotnet$prey,   # match prey
                                           groups$latin)])

cat(sum(is.na(hotnet$to)),  # print how many did not match
    "prey with missing matches:")
    hotnet[is.na(hotnet$to),]



# now we have to take away duplicates because more than one member of a group
# ate something from another group, so the link from group to group would be 
# included multiple times

hotnetsmall <- hotnet[!duplicated(hotnet[c("from","to", "type")]), 
                      c("from", "to", "type")]

#write.csv(hotnetsmall, file = "hotgrouped.csv", row.names = F)





