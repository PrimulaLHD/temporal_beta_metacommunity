beta_temporal_metacomm <- function(m1, m2, method = "Jac") {
  # m1 and m2 are the metacommunity matrix in time 1 and time 2, respectively
  # or species-habitat network in time 1, and in time 2, respectively
  if (all(dim(m1) != dim(m2))) {
    print("Error, m1 and m2 must have equal dimensions!")
  }
  
  if (all(dim(m1) == dim(m2))) {
    # share items
    D <- m1 - m2
    a <- sum(m1[m1 == m2]) #  number of shared items between m1 and m2
    b <- length(D[D == 1]) # number of unique items in m1
    c <- length(D[D == -1]) # number of unique items in m2
    
    
    # part 1, local components
    m1_0 <- bipartite::empty(m1)
    m2_0 <- bipartite::empty(m2)
    
    shared_sp1 <- colnames(m1_0) %in% colnames(m2_0)
    shared_patch1 <- rownames(m1_0) %in% rownames(m2_0)
    
    shared_sp2 <- colnames(m2_0) %in% colnames(m1_0)
    shared_patch2 <- rownames(m2_0) %in% rownames(m1_0)
    
    m1_00 <- (m1_0[shared_patch1, shared_sp1])
    m2_00 <- (m2_0[shared_patch2, shared_sp2])
    
    D0 <- m1_00 - m2_00
    b0 <-
      length(D0[D0 == 1]) # unique numbers of item in m1, local-extinction driven part
    c0 <-
      length(D0[D0 == -1]) # unique numbers of item in m2, local-colonization driven part
    
    # part 2, regional species changes
    
    R_mat_loss <- D[rowSums(m2) > 0, colSums(m2) == 0]
    b_R <-
      length(R_mat_loss[R_mat_loss == 1]) # regional species loss (extinction) driven part
    
    
    R_mat_gain <- D[rowSums(m1) > 0, colSums(m1) == 0]
    c_R <-
      length(R_mat_gain[R_mat_gain == -1]) # regional species gain (colonization) driven part
    
    # part 3, landscape modification
    
    L_mat_loss <- D[rowSums(m2) == 0, colSums(m2) > 0]
    b_L <- length(L_mat_loss[L_mat_loss == 1])
    
    L_mat_gain <- D[rowSums(m1) == 0, colSums(m1) > 0]
    c_L <- length(L_mat_gain[L_mat_gain == -1])
    
    # part 4, both landscape modification and regional species change driven part
    
    RL_mat_loss <- D[rowSums(m2) == 0, colSums(m2) == 0]
    b_RL <- length(RL_mat_loss[RL_mat_loss == 1])
    
    RL_mat_gain <- D[rowSums(m1) == 0, colSums(m1) == 0]
    c_RL <- length(RL_mat_gain[RL_mat_gain == -1])
    
    if (method == "Jac") {
      # Jaccard dissimilarity
      beta_Local <- (b0 + c0) / (a + b + c)
      beta_Regional <- (b_R + c_R) / (a + b + c)
      beta_Landscape <- (b_L + c_L) / (a + b + c)
      beta_RL <- (b_RL + c_RL) / (a + b + c)
      beta_temporal <-
        ((b + c) / (a + b + c)) # total temporal beta diversity of link types
      
      ## further partitioning
      beta_extinction <- (b0 + b_R + b_L + b_RL) / (a + b + c)
      beta_colonization <- (c0 + c_R + c_L + c_RL) / (a + b + c)
      beta_Local_colonization <- c0 / (a + b + c)
      beta_Local_extinction <- b0  / (a + b + c)
      beta_Regional_colonization <- c_R / (a + b + c)
      beta_Regional_extinction <- b_R / (a + b + c)
      beta_Landscape_gain <- c_L / (a + b + c)
      beta_Landscape_loss <- b_L / (a + b + c)
      beta_RL_colonization <- c_RL / (a + b + c)
      beta_RL_extinction <- b_RL / (a + b + c)
      
    }
    
    if (method == "Sor") {
      #Sørensen
      beta_Local <- (b0 + c0) / (2 * a + b + c)
      beta_Regional <- (b_R + c_R) / (2 * a + b + c)
      beta_Landscape <- (b_L + c_L) / (2 * a + b + c)
      beta_RL <- (b_RL + c_RL) / (2 * a + b + c)
      beta_temporal <-
        ((b + c) / (2 * a + b + c)) # total temporal beta diversity of link types
      
      beta_extinction <- (b0 + b_R + b_L + b_RL) / (2 * a + b + c)
      beta_colonization <- (c0 + c_R + c_L + c_RL) / (2 * a + b + c)
      
      ## further partitioning
      beta_extinction <- (b0 + b_R + b_L + b_RL) / (2 * a + b + c)
      beta_colonization <- (c0 + c_R + c_L + c_RL) / (2 * a + b + c)
      beta_Local_colonization <- c0 / (2 * a + b + c)
      beta_Local_extinction <- b0  / (2 * a + b + c)
      beta_Regional_colonization <- c_R / (2 * a + b + c)
      beta_Regional_extinction <- b_R / (2 * a + b + c)
      beta_Landscape_gain <- c_L / (2 * a + b + c)
      beta_Landscape_loss <- b_L / (2 * a + b + c)
      beta_RL_colonization <- c_RL / (2 * a + b + c)
      beta_RL_extinction <- b_RL / (2 * a + b + c)
      
      
    }
    
  }
  
  res <-
    data.frame(
      Component = c(
        "beta_Local",
        "beta_Regional",
        "beta_Landscape",
        "beta_RL",
        "beta_temporal",
        "beta_extinction",
        "beta_colonization",
        "beta_Local_colonization",
        "beta_Local_extinction",
        "beta_Regional_colonization",
        "beta_Regional_extinction",
        "beta_Landscape_gain",
        "beta_Landscape_loss",
        "beta_RL_colonization",
        "beta_RL_extinction"
      ),
      beta = c(
        beta_Local,
        beta_Regional,
        beta_Landscape,
        beta_RL,
        beta_temporal,
        beta_extinction,
        beta_colonization,
        beta_Local_colonization,
        beta_Local_extinction,
        beta_Regional_colonization,
        beta_Regional_extinction,
        beta_Landscape_gain,
        beta_Landscape_loss,
        beta_RL_colonization,
        beta_RL_extinction
      )
    )
  
  return(res)
  
}

## run example (Not run)
# m1 <- matrix(c(1,0,0,
#                0,1,1,
#                0,1,0), nrow = 3, ncol = 3, byrow = T)
# colnames(m1) <- c("A","B","C")
# rownames(m1) <- c("CC","DD","EE")
# m2 <- matrix(c(0,0,0,
#                0,1,1,
#                1,1,0), nrow = 3, ncol = 3, byrow = T)
# colnames(m2) <- c("A","B","C")
# rownames(m2) <- c("CC","DD","EE")
# beta_temporal_metacomm(m1 =m1, m2= m2)


