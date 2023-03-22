##generate species traits
library(mcomsimr)

set.seed(123)
narrow_optima <- env_traits(
  species = 70,
  max_r = 5,
  min_env = 0,
  max_env = 1,
  env_niche_breadth = 0.5,
  plot = TRUE,
  optima_spacing = "random"
)

set.seed(123)
flat_optima <- env_traits(
  species = 70,
  max_r = 5,
  min_env = 0,
  max_env = 1,
  env_niche_breadth = 10,
  plot = TRUE,
  optima_spacing = "random"
)
