## generate dispersal gradients in a log space
dispersal_gradient <- exp(seq(log(1e-5), log(1), length.out = 16))[1:15]
