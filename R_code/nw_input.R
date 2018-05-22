################################################
###Make an nw1 input file for n_w
################################################

hotnet<- read.csv("hotnet.csv", stringsAsFactors = F)

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
