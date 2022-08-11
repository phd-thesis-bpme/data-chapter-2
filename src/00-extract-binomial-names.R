####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 00-extract-binomial-names.R
# Created March 2022
# Last Updated August 2022

####### Import Libraries and External Files #######

library(magrittr)

####### Read Data #################################

#' This file comes directly from NA-POPS: analysis/01-combine-data.R
#' which is not publicly available. Will be loaded in as "project_counts"
load("data/raw/counts.rda")

codes <- read.csv("data/raw/IBP-Alpha-Codes20.csv")
families <- read.csv("data/raw/NACC_list_species.csv")
birdtree <- read.csv("data/raw/birdtree_taxonomy.csv")

####### Main Code #################################

# Extract all binomial names for species that NA-POPS has data for
codes_red <- codes[which(codes$SPEC %in% unique(project_counts$Species)), ]
families$species <- trimws(families$species, which = "right") # clean rogue whitespace
binomial <- merge(codes_red, families[, c("species", "family")], 
                      by.x = "SCINAME", by.y = "species")[, c("SPEC", "COMMONNAME", "SCINAME", "family")]
names(binomial) <- c("Code", "English", "Scientific", "Family")

#' Attempt to match up taxonomies with Birdtree taxonomy
#' First, perform a naive merge with common name, and determine which species were dropped
names(birdtree)[1] <- "Scientific_BT"
names(birdtree)[3] <- "English_BT"
birdtree_IOC <- birdtree[which(birdtree$Taxonomy == "IOC27"), ]
birdtree <- birdtree[which(birdtree$Taxonomy == "BL3"), ]
binomial_matched_common <- merge(binomial, birdtree[, c("Scientific_BT", "English_BT", "Taxonomy")],
                        by.x = "English", by.y = "English_BT")
dropped <- setdiff(binomial$English, binomial_matched_common$English)

#' Now, merge a subsetted dataframe that contains only the dropped species, 
#' with the birdtree dataframe, this time by scientific name. You will end up
#' with two English names (likely differing by only one letter, ugh). From here,
#' merge back in the BT scientific names by matching with "English_BT". Finally,
#' bind these new rows back into "binomial_matched_common" from before". Take the
#' difference between this new dataframe and the original "binomial" dataframe
binomial_dropped <- binomial[which(binomial$English %in% dropped), ]
binomial_matched_sci <- merge(binomial_dropped, birdtree[, c("Scientific_BT", "English_BT", "Taxonomy")],
                              by.x = "Scientific", by.y = "Scientific_BT") %>%
  merge(birdtree[,c("Scientific_BT", "English_BT")],
        by = "English_BT")

binomial_matched <- rbind(binomial_matched_common,
                          binomial_matched_sci[, c("English", "Code", "Scientific", "Family", "Scientific_BT", "Taxonomy")])
dropped <- setdiff(binomial$English, binomial_matched$English)

#' Now, check for any that are in IOC27 taxonomy
binomial_dropped <- binomial[which(binomial$English %in% dropped), ]
binomial_matched_ioc <- merge(binomial_dropped, birdtree_IOC[, c("Scientific_BT", "English_BT", "Taxonomy")],
                                 by.x = "English", by.y = "English_BT")
if (nrow(binomial_matched_ioc > 0))
{
  binomial_matched <- rbind(binomial_matched,
                            binomial_matched_ioc[, c("English", "Code", "Scientific", "Family", "Scientific_BT", "Taxonomy")])
}
dropped <- setdiff(binomial$English, binomial_matched$English)

binomial_dropped <- binomial[which(binomial$English %in% dropped), ]
binomial_matched_sci_ioc <- merge(binomial_dropped, birdtree_IOC[, c("Scientific_BT", "English_BT", "Taxonomy")],
                              by.x = "Scientific", by.y = "Scientific_BT") %>%
  merge(birdtree_IOC[,c("Scientific_BT", "English_BT")],
        by = "English_BT")
if (nrow(binomial_matched_sci_ioc > 0))
{
  binomial_matched <- rbind(binomial_matched,
                            binomial_matched_sci_ioc[, c("English", "Code", "Scientific", "Family", "Scientific_BT", "Taxonomy")])
}
dropped <- setdiff(binomial$English, binomial_matched$English)

#' We will now just have to do the remaining species by hand. Previous tests of
#' fuzzy matching was okay for some, but then missed some others, so just best to
#' do it manually unfortunately. This section here will likely have to be updated
#' from time to time as more species get added to NA-POPS. I'll also put an explanation
#' beside each species explaning why it must be done manually.
binomial_remaining <- binomial[which(binomial$English %in% dropped), 
                               c("English", "Code", "Scientific", "Family")]
binomial_remaining$Scientific_BT <- NA
binomial_remaining$Taxonomy <- NA
print(binomial_remaining$English)

#' LeConte's Sparrow couldn't fuzzily match well with "Le Conte's Sparrow, without
#' including another wrong match
binomial_remaining[which(binomial_remaining$English == "LeConte's Sparrow"),
                   "Scientific_BT"] <- 
  birdtree[which(birdtree$English_BT == "Le Conte's Sparrow"), "Scientific_BT"]
binomial_remaining[which(binomial_remaining$English == "LeConte's Sparrow"),
                   "Taxonomy"] <- 
  birdtree[which(birdtree$English_BT == "Le Conte's Sparrow"), "Taxonomy"]

#' Eastern Whip-poor-will goes with BirdTree's Whip-poor-will. Fuzzy matching coulnd't
#' handle
binomial_remaining[which(binomial_remaining$English == "Eastern Whip-poor-will"),
                   "Scientific_BT"] <- 
  birdtree[which(birdtree$English_BT == "Whip-poor-will"), "Scientific_BT"]
binomial_remaining[which(binomial_remaining$English == "Eastern Whip-poor-will"),
                   "Taxonomy"] <- 
  birdtree[which(birdtree$English_BT == "Whip-poor-will"), "Taxonomy"]

#' Gray Hawk is "Grey Hawk" in BirdTree, and different scientific name, easy enough
binomial_remaining[which(binomial_remaining$English == "Gray Hawk"),
                   "Scientific_BT"] <- 
  birdtree[which(birdtree$English_BT == "Grey Hawk"), "Scientific_BT"]
binomial_remaining[which(binomial_remaining$English == "Gray Hawk"),
                   "Taxonomy"] <- 
  birdtree[which(birdtree$English_BT == "Grey Hawk"), "Taxonomy"]

#' Eastern Yellow Wagtail should just be "Yellow Wagtail", and diff sci name
binomial_remaining[which(binomial_remaining$English == "Eastern Yellow Wagtail"),
                   "Scientific_BT"] <- 
  birdtree[which(birdtree$English_BT == "Yellow Wagtail"), "Scientific_BT"]
binomial_remaining[which(binomial_remaining$English == "Eastern Yellow Wagtail"),
                   "Taxonomy"] <- 
  birdtree[which(birdtree$English_BT == "Yellow Wagtail"), "Taxonomy"]

#' Thick-billed Longspur got changed from McCown's Longspur
binomial_remaining[which(binomial_remaining$English == "Thick-billed Longspur"),
                   "Scientific_BT"] <- 
  birdtree[which(birdtree$English_BT == "McCown's Longspur"), "Scientific_BT"]
binomial_remaining[which(binomial_remaining$English == "Thick-billed Longspur"),
                   "Taxonomy"] <- 
  birdtree[which(birdtree$English_BT == "McCown's Longspur"), "Taxonomy"]

#' Black-throated Gray Warbler is Black-throated Grey Warbler in BirdTree
#' with a different sci name
binomial_remaining[which(binomial_remaining$English == "Black-throated Gray Warbler"),
                   "Scientific_BT"] <- 
  birdtree[which(birdtree$English_BT == "Black-throated Grey Warbler"), "Scientific_BT"]
binomial_remaining[which(binomial_remaining$English == "Black-throated Gray Warbler"),
                   "Taxonomy"] <- 
  birdtree[which(birdtree$English_BT == "Black-throated Grey Warbler"), "Taxonomy"]

#' Woodhouse's Scrub-Jay, Bell's Sparrow, Sagebrush Sparrow, Pacific Wren 
#' didn't exist at the time of Jetz et al. 2012 so just drop.
#' Northern Shrike apparently was just lumped with "Great Gray Shrike" in birdtree?
#' I don't quite understand this one, but I'm just going to drop it as well.
binomial_remaining <- binomial_remaining[-which(is.na(binomial_remaining$Scientific_BT)), ]

# Glue these on to the rest of the matches, and you've got your matched table
binomial_matched <- rbind(binomial_matched, binomial_remaining)

# Now add Bicknell's Thrush as "new data"
binomial_matched <- rbind(binomial_matched,
                          data.frame(English = "Bicknell's Thrush",
                                     Code = "BITH",
                                     Scientific = birdtree[which(birdtree$English_BT == "Bicknell's Thrush"), "Scientific_BT"],
                                     Family = "Turdidae",
                                     Scientific_BT = birdtree[which(birdtree$English_BT == "Bicknell's Thrush"), "Scientific_BT"],
                                     Taxonomy = birdtree[which(birdtree$English_BT == "Bicknell's Thrush"), "Taxonomy"]))
####### Output ####################################

write.table(binomial_matched,
            file = "data/generated/binomial_names.csv",
            sep = ",",
            row.names = FALSE)
