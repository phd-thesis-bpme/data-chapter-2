distance: distance-prior distance-model

distance-prior:
	setsid R CMD BATCH ./src/06-distance-prior-check.R &

distance-model:
	setsid R CMD BATCH ./src/07-distance-model.R &
