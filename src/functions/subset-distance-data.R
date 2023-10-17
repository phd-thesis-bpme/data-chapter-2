subset_distance_data <- function(distance_stan_data = NULL,
                                 sps = NULL)
{
  library(dplyr)
  orig_data_df <- data.frame(species = distance_stan_data$sp_list,
                             sp_n = distance_stan_data$species)

  sp_df <- orig_data_df %>%
    group_by(species,sp_n) %>%
    dplyr::summarise(n_counts = dplyr::n())

  sub_data <- data.frame(mig_strat = distance_stan_data$mig_strat,
                         habitat = distance_stan_data$habitat,
                         pitch = distance_stan_data$pitch,
                         mass = distance_stan_data$mass,
                         sp_n = 1:length(distance_stan_data$mig_strat))
  sp_df <- dplyr::left_join(sp_df,sub_data)

  sp_sel <- sp_df[which(sp_df$species %in% sps),]
  sp_sel[,"new_sp_n"] <- 1:length(sps) # new species indicators

  #identify which of the original data objects should be retained
  data_sel <- which(orig_data_df$sp_n %in% sp_sel$sp_n)
  sp_data_repl <- orig_data_df %>%
    left_join(.,sp_sel,
              by = "sp_n")

  distance_stan_data2 <- distance_stan_data; rm(distance_stan_data); gc();
  distance_stan_data2$bands_per_sample <- distance_stan_data2$bands_per_sample[data_sel]
  distance_stan_data2$abund_per_band <- distance_stan_data2$abund_per_band[data_sel,]
  distance_stan_data2$max_dist <- distance_stan_data2$max_dist[data_sel,]
  distance_stan_data2$species <- sp_data_repl[data_sel,"new_sp_n"]
  distance_stan_data2$mig_strat <- as.integer(unlist(sp_sel[,"mig_strat"]))
  distance_stan_data2$habitat <- as.integer(unlist(sp_sel[,"habitat"]))
  distance_stan_data2$mass <- as.numeric(unlist(sp_sel[,"mass"]))
  distance_stan_data2$pitch <- as.numeric(unlist(sp_sel[,"pitch"]))

  distance_stan_data2$n_species <- length(sps)
  distance_stan_data2$n_samples <- length(data_sel)
  distance_stan_data2$sp_list <- distance_stan_data2$sp_list[which(distance_stan_data2$sp_list %in% sps)]

  return(distance_stan_data2)
}