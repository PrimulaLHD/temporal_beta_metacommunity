## generating species competition scenario
library(mcomsimr)
eq_int <- species_int_mat(
  species = 70,
  intra = 1,
  min_inter = 1,
  max_inter = 1,
  comp_scaler = 0.05,
  plot = TRUE
)

set.seed(123)
stab_int <- species_int_mat(
  species = 70,
  intra = 1,
  min_inter = 0,
  max_inter = 0.5,
  comp_scaler = 0.05,
  plot = TRUE
)
