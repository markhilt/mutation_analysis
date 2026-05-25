# Plot mutations as their frequency of reads.

require(ggplot2)
require(stringr)
require(reshape2)
require(dplyr)

setwd("Y:/Research/Germline_project/mutations_analysis/plots")

import_data <-function(tabfile){
  dat <- read.delim(tabfile, check.names=TRUE)
  df <- select(dat, -CHROM,-POS,-ID,-REF,-ALT,-QUAL)  # Extract interesting fields and reformat
  
  # Calculate the frequencies in each element
  calcfreq <- function(x){
    vals <- unlist(strsplit(x,","))
    ref <- as.numeric(vals[1])
    alt <- as.numeric(vals[2])
    tot <- ref+alt
    freq <- alt/tot
    return(freq)
  }
  freq <- as.data.frame(apply(df, c(1,2), calcfreq))
  # Keep track of mean frequency per mutation - this is the plotting order
  freq$mutation <- rowMeans(freq)
  freq <- freq[order(-freq$mutation),] 
  freq$mutation <- paste0(rep("mut",nrow(freq)),1:nrow(freq))
  dat2 <- melt(freq,id.vars="mutation",variable.name="Sample",value.name="Frequency")
  return(dat2)
}

# Import dataset. The file has GATK VariantsToTable format.
f <- "Y:/Research/Germline_project/mutations_analysis/Stipe-gills_mutations.curated.txt"
dat <- import_data(f)

# Reformat data frame
dat$Tissue <- ifelse(grepl("S$", dat$Sample, ignore.case = T), "Stipe", "Gill")
dat$Mushroom <- str_remove(str_remove(dat$Sample, "[LS]$"),"X")
# mutation_order arbitrarily determined. Adapt as needed.
mutation_order <- c("mut1","mut4","mut5","mut6","mut7","mut8","mut9","mut19",
"mut10","mut15","mut12","mut13","mut14","mut16","mut20","mut24",
"mut21","mut23","mut25","mut27",
"mut28","mut29","mut30","mut11","mut22","mut26","mut17","mut18","mut2","mut3")

p <- ggplot(data=dat, aes(x=mutation, y=Frequency, fill=mutation)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_x_discrete(limits = mutation_order) +
  ggtitle("Mutations in Stipe-Gills") +
  theme(axis.text.x=element_blank(),
        legend.position = "none") +
  facet_grid(vars(Tissue),vars(Mushroom))

ggsave(file="stipe-gills.ring06.svg", plot=p, width=10, height=6)
