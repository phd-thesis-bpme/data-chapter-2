####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 00-extract-binomial-names.R
# Created March 2022
# Last Updated August 2022

####### Import Libraries and External Files #######

library(magrittr)
library(fedmatch)
library(data.table)

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

#' Now the remaining dropped species will be a bit trickier. For the most part, 
#' it will likely be the case that the names are just spelled slightly different,
#' AND the scientific name differs (aren't those NOT suppose to change? Ugh...).
#' So let's try fuzzy matching by Common name first
#' SPECIES THAT WILL NOT/SHOULD NOT WORK:
#' * Sagebrush Sparrow/Sage Sparrow (SBSP used to be SAGS with BESP)
#' * Eastern Yellow Wagtail (Just yellow wagtail)
#' * Pacific Wren (may have to drop)
#' * Woodhouse's Scrub Jay (investigate)
#' * Eastern Whip-poor-will (just Whip-poor-will)
#' * Thick-billed Longspur (used to be Smith's Longspur)
#' * Northern Shrike (????)
#' * Bell's sparrow (same explanation as first on this list)


# fuzzy_match <- merge_plus(data1 = binomial_dropped,
#                           data2 = birdtree[, c("Scientific_BT", "English_BT")],
#                           by.x = "English", by.y = "English_BT",
#                           match_type = "fuzzy",
#                           unique_key_1 = "Scientific",
#                           unique_key_2 = "Scientific_BT",
#                           fuzzy_settings = build_fuzzy_settings(maxDist = .08))
#binomial[which(binomial$Scientific == "Ixoreus naevius"), "Scientific"] <- "Zoothera naevia"


####### Output ####################################

write.table(binomial,
            file = "data/generated/binomial_names.csv",
            sep = ",",
            row.names = FALSE)
