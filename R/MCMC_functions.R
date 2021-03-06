# prepare ECs and USA pseudo-bulk counts:
prepare_PB_counts = function(one_cluster_ids_kept, sce, clusters, n_samples,
                             EC_list = NULL){
  # select cell-types
  sce = sce[,clusters == one_cluster_ids_kept]
  
  # compute pseudo-bulk SUA counts for each sample:
  S = U = A = matrix(0, nrow = nrow(sce), ncol = n_samples)
  for(i in seq_len(n_samples)){
    sel_cells = sce$sample_id == levels(sce$sample_id)[i]
    sce_one_sample = sce[, sel_cells]
    S[,i] = rowSums(assays(sce_one_sample)$spliced)
    U[,i] = rowSums(assays(sce_one_sample)$unspliced)
    A[,i] = rowSums(assays(sce_one_sample)$ambiguous)
  }
  #S = Matrix(data=S, sparse = TRUE)
  #U = Matrix(data=U, sparse = TRUE)
  #A = Matrix(data=A, sparse = TRUE)
  
  if(!is.null(EC_list)){
    EC_counts = EC_list[[1]]
    # get cells with selected cell types in sce
    cells_sel = colnames(sce)
    # select cells of cell type cl
    EC_counts = lapply(EC_counts, function(x){
      sel = rownames(x) %in% cells_sel
      x[sel,]
    })
    # compute pseudo-bulk counts aggregating counts across cells:
    counts = lapply(EC_counts, colSums)
    # make counts a sparse Vector:
    #counts = lapply(counts, function(x){
    #  sel = x > 0.5
    #  sparseVector(x[sel], which(sel), length(x))
    #})
  }else{
    counts = NULL
  }
  
  list(counts, S, U, A)
}

find_mode <- function(x, adjust, ...) {
  dx <- density(x, adjust = adjust, ...)
  dx$x[which.max(dx$y)]
}

# compute gene-level p-value:
compute_pval = function(A, B, K = 3, N){
  R = nrow(A)
  A = A[sample.int(R, R),] # n indicates the nr of elements of the chain (exluded burn-in)
  
  gamma = A - B
  
  CV  = cov(gamma) # cov is 20ish times faster than posterior mode (very marginal cost).
  mode = apply(gamma, 2, find_mode, adjust = 10)
  
  p = K-1
  
  p_value = vapply(seq_len(K), function(k){
    sel  = seq_len(K)[-k]
    # Normal (classical Wald test)
    stat = t(mode[sel]) %*% ginv(CV[sel, sel], tol = 0) %*% mode[sel]
    1-pchisq(stat, df = K-1)
  }, FUN.VALUE = numeric(1))
  
  # USA posterior mean and SD of both conditions, A and B.
  mode_A_USA = colSums(A) # find.mode (mode) or sum (mean)
  mode_A_USA = mode_A_USA/sum(mode_A_USA)
  sd_A_USA = sqrt(diag(var(A)))
  mode_B_USA = colSums(B) # find.mode (mode) or sum (mean)
  mode_B_USA = mode_B_USA/sum(mode_B_USA)
  sd_B_USA = sqrt(diag(var(B)))
  
  # US posterior mean and SD of both conditions, A and B.
  A[,1] = A[,1] + 0.5 * A[,3]
  A[,2] = A[,2] + 0.5 * A[,3]
  A = A[, seq_len(2)]
  
  B[,1] = B[,1] + 0.5 * B[,3]
  B[,2] = B[,2] + 0.5 * B[,3]
  B = B[, seq_len(2)]
  
  mode_A = colSums(A) # find.mode (mode) or sum (mean)
  mode_A = mode_A/sum(mode_A)
  sd_A = sqrt(diag(var(A)))
  mode_B = colSums(B) # find.mode (mode) or sum (mean)
  mode_B = mode_B/sum(mode_B)
  sd_B = sqrt(diag(var(B)))
  
  c( mean(p_value), # p-value for Diff. Reg.
    mode_A, mode_B, sd_A, sd_B, # posterior mean and SD for US pi
    mode_A_USA, mode_B_USA, sd_A_USA, sd_B_USA) # posterior mean and SD for USA pi
}

# convergence diagnostic:
my_heidel_diag = function(x, R, by., pvalue = 0.01){
  start.vec <- seq(from = 1, to = R/2, by = by.)
  S0 <- my_spectrum0_ar(window(x, start = R/2), R/2+1)
  
  converged <- FALSE
  for(i in seq(along = start.vec)){
    x <- window(x, start = start.vec[i])
    n <- R + 1 - start.vec[i] # niter(x)
    B <- cumsum(x) - sum(x) * seq_len(n)/n
    Bsq <- (B * B)/(n * S0)
    I <- sum(Bsq)/n
    p = my_pcramer_(I)
    if(converged <- !is.na(I) && p < 1 - pvalue){
      break
    }
  }
  
  if( !converged || is.na(I) ) {
    nstart <- NA
  }else {
    nstart <- start.vec[i]
  }
  return(c(converged, nstart, 1 - p))
}

# additional function for my_heidel_diag
my_spectrum0_ar = function(x, R){
  lm.out <- lm(x ~ seq_len(R) )
  if(identical(all.equal(sd(residuals(lm.out)), 0), TRUE)) {
    v0 <- 0
  }else{
    ar.out <- ar(x, aic = TRUE)
    v0 <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
  }
  #  return(list(spec = v0, order = order))
  v0
}

# additional function for my_heidel_diag
my_pcramer_ = function(q, eps = 1e-05){
  log.eps <- log(eps)
  y = vapply(seq(0, 3, by = 1), function(k){
    z <- gamma(k + 0.5) * sqrt(4 * k + 1)/(gamma(k + 1) * 
                                             pi^(3/2) * sqrt(q))
    u <- (4 * k + 1)^2/(16 * q)
    ifelse(u > -log.eps, 0, z * exp(-u) * besselK(x = u, 
                                                  nu = 1/4))
  }, FUN.VALUE = numeric(1))
  return(sum(y))
}

