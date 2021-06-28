process_spec_diff <- function(this_species, hex_data, ebird_specs, write_inds = T, 
                              verbose = F, overwrite = F, timeout = Inf, maxK = 20,
                              output_path = "intermediate_data/ind_results/",
                              warn_path = "intermediate_data/warnings/") {
  
  withTimeout({
    this_name_clean <- ebird_specs$name_clean[ebird_specs$SCIENTIFIC.NAME == this_species]
    
    flog.appender(appender.file(paste0(warn_path, this_name_clean, "_warning.log")))
    this_output_file <- paste0(output_path, "res_", this_name_clean, ".csv")
    if (file.exists(this_output_file) && !overwrite) {
      return(NULL)
    }
    
    target_df <- data.frame(
      species = this_species,
      common_name = this_name_clean,
      difference = NA,
      diff_SE = NA,
      k = NA,
      inat_medPred = NA,
      inat_SE = NA,
      inat_medResp = NA,
      inat_GAM_pval = NA,
      ebird_medPred = NA,
      ebird_SE = NA,
      ebird_medResp = NA,
      ebird_GAM_pval = NA
    )
    
    if (verbose) {
      cat("Processing", this_species, "-", this_name_clean, "\n")
    }
    
    this_dat_ebird <- hex_data %>% 
      filter(species == this_species, dataset == "ebird_n") %>% 
      mutate(lon = scale(lon), lat = scale(lat))
    this_dat_inat <- hex_data %>% 
      filter(species == this_species, dataset == "inat_n") %>% 
      mutate(lon = scale(lon), lat = scale(lat))
    
    temp_k <- min(floor(sqrt(nrow(this_dat_ebird))), maxK)
    
    this_fit_ebird <- gam(cbind(success, total - success) ~ te(lon, lat, k = temp_k), 
                          family = quasibinomial(link = "logit"), 
                          optimizer = c("outer", "bfgs"),
                          control = list(maxit = 100000),
                          data = this_dat_ebird)
    this_fit_inat <- gam(cbind(success, total - success) ~ te(lon, lat, k = temp_k), 
                         family = quasibinomial(link = "logit"),
                         optimizer = c("outer", "bfgs"),
                         control = list(maxit = 100000),
                         data = this_dat_inat)

    check_ebird = k.check(this_fit_ebird)
    check_inat = k.check(this_fit_inat)
    
    target_df$k = temp_k
    target_df$ebird_GAM_pval = check_ebird[4]
    target_df$inat_GAM_pval = check_inat[4]
    
    #FYI check_ebird[2] is edf, check_ebird[3] is kindex
    
    
    rmvn <- function(n,mu,sig) { ## MVN random deviates
      L <- mroot(sig);m <- ncol(L);
      t(mu + L%*%matrix(rnorm(m*n),m,n)) 
    }
    
    Xp_ebird <- predict(this_fit_ebird ,type="lpmatrix") 
    Xp_inat  <- predict(this_fit_inat ,type="lpmatrix") 
    
    br_ebd <- rmvn(10000, coef(this_fit_ebird),this_fit_ebird$Vp) ## 1000 replicate param. vectors
    br_inat <- rmvn(10000, coef(this_fit_inat),this_fit_inat$Vp) ## 1000 replicate param. vectors
    
    res1 <- rep(0,10000)
    res2 <- rep(0,10000)
    res3 <- rep(0,10000)
    for (i in 1:10000) { 
      pr_ebd <- Xp_ebird %*% br_ebd[i,] ## replicate predictions
      pr_inat <- Xp_inat %*% br_inat[i,] ## replicate predictions
      res1[i] <- median(pr_ebd) ## median eBird prediction
      res2[i] <- median(pr_inat) ## median iNat prediction
      res3[i] <- median(pr_inat - pr_ebd) ## median difference
    }

    
    saveRDS(list(
      ebird = Xp_ebird, inat = Xp_inat, 
      ebird_coef = coef(this_fit_ebird),
      inat_coef = coef(this_fit_inat),
      species = this_species, name_clean = this_name_clean
    ), paste0("intermediate_data/lpmtx/lp_matrix_", this_name_clean, ".RDS"))
    
    
    target_df$ebird_medPred <- median(res1)
    target_df$inat_medPred <- median(res2)
    target_df$difference <- median(res3)
    
    target_df$ebird_medResp <- nimble::expit(median(res1))
    target_df$inat_medResp <- nimble::expit(median(res2))
    
    target_df$ebird_SE <- sd(res1)
    target_df$inat_SE <- sd(res2)
    target_df$diff_SE <- sd(res3)
  
    if (write_inds) write_csv(target_df, this_output_file)
    target_df
    
  }, timeout = timeout)
}
