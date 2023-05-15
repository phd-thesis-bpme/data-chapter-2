analysis: removal distance

removal: prepare-removal removal-prior
	setsid R CMD BATCH ./src/06-removal-model.R &

removal-prior: data/generated/removal_stan_data_pred.rda
	setsid R CMD BATCH ./src/03-removal-prior-check.R &

prepare-removal: data/generated/corr_matrix_predict.rda data/generated/binomial_names.csv data/raw/traits.csv data/raw/time_count_matrix.rda data/raw/time_design.rda
	setsid R CMD BATCH ./src/02-prepare-removal-data.R &

data/generated/corr_matrix_predict.rda: data/raw/all_species.nex data/generated/binomial_names.csv
	setsid R CMD BATCH ./src/01-prepare-phylogenies.R &

data/generated/binomial_names.csv: data/raw/counts.rda data/raw/IBP-Alpha-Codes20.csv data/raw/NACC_list_species.csv data/raw/birdtree_taxonomy.csv
	setsid R CMD BATCH ./src/00-extract-binomial-names.R &

distance: data/generated/distance_stan_data_pred.rda distance-prior
	setsid R CMD BATCH ./src/07-distance-model.R &

distance-prior: data/generated/distance_stan_data_pred.rda
	setsid R CMD BATCH ./src/05-distance-prior-check.R &

data/generated/distance_stan_data_pred.rda: data/raw/dist_count_matrix.rda data/raw/dist_design.rda data/generated/corr_matrix_predict.rda data/generated/binomial_names.csv data/raw/traits.csv
	setsid R CMD BATCH ./src/04-prepare-distance-data.R &

