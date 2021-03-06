# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Rcpp_MCMC <- function(n_samples, n_genes, n_groups, n_genes_keep, keep_genes_id, numeric_groups, sample_ids_per_group, n_samples_per_group, N_MCMC, burn_in, PI_gene, PI_SU, list_X_unique, list_EC_gene_id, list_EC_SU_id, counts, MCMC_bar_pi_1, MCMC_bar_pi_2, MCMC_bar_pi_3, chol, delta_SU, prior_TF, mean_log_delta, sd_log_delta, sd_prior_non_informative) {
    .Call(`_DifferentialRegulation_Rcpp_MCMC`, n_samples, n_genes, n_groups, n_genes_keep, keep_genes_id, numeric_groups, sample_ids_per_group, n_samples_per_group, N_MCMC, burn_in, PI_gene, PI_SU, list_X_unique, list_EC_gene_id, list_EC_SU_id, counts, MCMC_bar_pi_1, MCMC_bar_pi_2, MCMC_bar_pi_3, chol, delta_SU, prior_TF, mean_log_delta, sd_log_delta, sd_prior_non_informative)
}

Rcpp_MCMC_sce <- function(n_samples, n_genes, n_groups, numeric_groups, sample_ids_per_group, n_samples_per_group, N_MCMC, burn_in, PI_SU, X_sce, MCMC_bar_pi_1, MCMC_bar_pi_2, MCMC_bar_pi_3, chol, delta_SU, prior_TF, mean_log_delta, sd_log_delta, sd_prior_non_informative) {
    .Call(`_DifferentialRegulation_Rcpp_MCMC_sce`, n_samples, n_genes, n_groups, numeric_groups, sample_ids_per_group, n_samples_per_group, N_MCMC, burn_in, PI_SU, X_sce, MCMC_bar_pi_1, MCMC_bar_pi_2, MCMC_bar_pi_3, chol, delta_SU, prior_TF, mean_log_delta, sd_log_delta, sd_prior_non_informative)
}

