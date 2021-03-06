#include <cassert>
#include <ctime>
#include <iostream>
#include <random>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

// Armadillo
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
//using namespace arma;

// "> (btw, in general it's more efficient to access matrices as columns rather than rows)."

void prior_informative(double& prior,
                       Rcpp::NumericVector const& x,
                       double const& mean, 
                       double const& sd,
                       double const& sd_prior_non_informative,
                       arma::mat mat){
  double sum_exp_x = sum(exp(x));
  double log_disp = log(sum_exp_x);
  
  mat.zeros(); 
  // can I initialize a matrix like this with all 0's ???
  for (unsigned int k = 0; k < 3; ++k) {
    mat(k,k) = 1;
    mat(2,k) = exp(x[k])/sum_exp_x;
  }
  
  double sign;
  arma::log_det(prior, sign, mat);
  
  prior += R::dnorm(log_disp, mean, sd, true) + R::dnorm(x[0], mean - log(3), sd_prior_non_informative, true) + R::dnorm(x[1], mean - log(3), sd_prior_non_informative, true);
}


double
  ll_alpha(Rcpp::NumericMatrix const& pi, 
           Rcpp::NumericVector const& alpha,
           unsigned int const& N)
  {
    double ll = 0.0;
    double sum_alpha = 0.0;
    
    for (unsigned int k = 0; k < 3; ++k) {
      sum_alpha += alpha[k];
      
      ll -= N * R::lgammafn(alpha[k]);
      
      for (unsigned int i = 0; i < N; ++i) {
        ll += log(pi(i,k)) * ( alpha[k] - 1);
      }
    }
    ll += N * R::lgammafn(sum_alpha);
    
    if( ISNAN(ll) ){
      Rcout << "na ll is" << std::endl << ll << std::endl;
      ll = -std::numeric_limits<double>::infinity();
      Rcout << "ll becomes" << std::endl << ll << std::endl;
    }
    
    return ll;
  }


void covRcpp_bis(Rcpp::NumericMatrix& Y,
                 Rcpp::NumericMatrix& cov,
                 double const& c_diag,
                 double const& c_prop, 
                 unsigned int const& df) {
  // Centering the matrix:
  Y(_, 0) = Y(_, 0) - mean(Y(_, 0));
  Y(_, 1) = Y(_, 1) - mean(Y(_, 1));
  Y(_, 2) = Y(_, 2) - mean(Y(_, 2));
  
  // COV only updates the bottom right corner element!
  
  // Computing the covariance matrix
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j <= i; ++j) { // I REMOVE THE FIRST 100 ROWS OF THE MCMC AND START FROM THE 101-st row.
      cov(i,j) = c_prop * Rcpp::sum(Y(Rcpp::_, i)*Y(Rcpp::_, j))/df;
      cov(j,i) = cov(i,j);
    }
    cov(i,i) += c_diag;
  }
}

void my_rmvnorm_final(Rcpp::NumericVector& alpha_prop,
                      arma::mat const& chol_R, 
                      Rcpp::NumericVector const& mean,
                      arma::vec& sample_arma){
  sample_arma = Rcpp::as<arma::vec>(Rcpp::rnorm(3));
  sample_arma = trans(sample_arma.t() * chol_R);
  
  alpha_prop = Rcpp::wrap(sample_arma);
  alpha_prop += mean;
}


// [[Rcpp::export]]
Rcpp::NumericVector Rcpp_MCMC( unsigned int const& n_samples, // N samples
                               unsigned int const& n_genes, // N genes
                               unsigned int const& n_groups, // N groups
                               unsigned int const& n_genes_keep, // TF vector indicating genes to be analyzed (SU differential testing)
                               Rcpp::IntegerVector const& keep_genes_id, // group id for every sample (must start from 0)
                               Rcpp::IntegerVector const& numeric_groups, // group id for every sample (must start from 0)
                               Rcpp::ListOf<Rcpp::IntegerVector> const& sample_ids_per_group, // each list = vector with ids of samples 
                               Rcpp::IntegerVector const& n_samples_per_group,
                               unsigned int const& N_MCMC, // MCMC iter
                               unsigned int const& burn_in, // burn-in
                               Rcpp::NumericMatrix& PI_gene, // prob of each gene (for every sample)
                               Rcpp::ListOf<Rcpp::NumericMatrix>& PI_SU, // prob of each gene (for every sample)
                               Rcpp::ListOf<Rcpp::NumericMatrix> const& list_X_unique, // SU uniquely mapping counts
                               Rcpp::ListOf<Rcpp::ListOf<Rcpp::IntegerVector>> const& list_EC_gene_id, // SU uniquely mapping counts
                               Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector>> const& list_EC_SU_id, // TRUE -> S; FALSE -> U
                               Rcpp::ListOf<Rcpp::IntegerVector> const& counts, // EC counts (integers)
                               Rcpp::ListOf<Rcpp::NumericMatrix>& MCMC_bar_pi_1,
                               Rcpp::ListOf<Rcpp::NumericMatrix>& MCMC_bar_pi_2,
                               Rcpp::ListOf<Rcpp::NumericMatrix>& MCMC_bar_pi_3,
                               Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericMatrix>>& chol,
                               Rcpp::ListOf<Rcpp::NumericMatrix>& delta_SU,
                               bool const& prior_TF,
                               double const& mean_log_delta,
                               double const& sd_log_delta,
                               double const& sd_prior_non_informative
) {
  // Obtain environment containing chol function:
  Rcpp::Environment base("package:base"); 
  // Make chol function callable from C++precision$prior
  Rcpp::Function R_chol = base["chol"];    
  
  for (unsigned int gr=0 ; gr < n_groups; gr++) {
    Rcpp::NumericMatrix tmp_1(N_MCMC, n_genes_keep), tmp_2(N_MCMC, n_genes_keep), tmp_3(N_MCMC, n_genes_keep);
    MCMC_bar_pi_1[gr] = tmp_1;
    MCMC_bar_pi_2[gr] = tmp_2;
    MCMC_bar_pi_3[gr] = tmp_3;
  }
  double tmp_1, tmp_2, tmp_3;
  
  // X is a temporary matrix to update the (s,u, a) reads of each sample:
  Rcpp::NumericMatrix X(n_genes, 3), prior_delta_SU(n_genes_keep, n_groups), cv_alpha(3, 3), pi_SU_sample(n_genes, 3), Y(N_MCMC, 3);
  
  // keep a temporary vectorized form of pi_S to speed-up memory access
  Rcpp::NumericVector pi_gene_sample(n_genes), X_gene(n_genes), prop_delta_SU(3), original_delta(3), exp_prop_delta_SU(3), MCMC_ll(N_MCMC), pi_SU_tmp(3);
  
  unsigned int EC_len, iter, sample, ec, i, gene_id, gr, cond_01, id, df, gene_to_keep_id;
  double prob_tot, dir_sample_SUA_sum, prior_prop_delta_SU, alpha;
  //ll_delta_diff
  
  bool cond, cond_SUA;
  SEXP tmp_chol; //tmp object to store cholesky matrix (computed from R function).
  
  // constant to add to the diagonal of the ARW matrix.
  double c_diag = 0.001;// the proportionality constant for the ARW matrix.
  double c_prop = 0.2; //0.3;
  
  arma::vec sample_arma(3); // arma for ARW proposal
  arma::mat arma_mat_prior(3,3);
  
  
  // initalize prior for delta_SU:
  for (gr=0 ; gr < n_groups; gr++) {
    for (gene_id=0 ; gene_id < n_genes_keep; gene_id++) {
      gene_to_keep_id = keep_genes_id[gene_id];
      
      original_delta = delta_SU[gr]( gene_to_keep_id,_);
      
      if(gene_id < n_genes_keep){
        // Initialize the prior:
        if( prior_TF ){ // if prior has been specified
          prior_informative(prior_prop_delta_SU, log(original_delta),  mean_log_delta, sd_log_delta, sd_prior_non_informative, arma_mat_prior);
          prior_delta_SU(gene_id, gr) = prior_prop_delta_SU;
        }else{ // if prior has NOT been specified: TO DO: compute prior once only at the beginning!
          prior_delta_SU(gene_id, gr) = R::dnorm( log(original_delta[0]), 0, sd_prior_non_informative, true) + R::dnorm( log(original_delta[1]), 0, sd_prior_non_informative, true) + R::dnorm( log(original_delta[2]), 0, sd_prior_non_informative, true);
        }
      }
    }
    
  }
  
  for (iter=0 ; iter < N_MCMC; iter++) {
    
    //Rcout << "iter" << std::endl << iter << std::endl;
    
    Rcpp::checkUserInterrupt(); // interrupt code from R (benchmark cost).
    
    //Rcout << "starting iter: ------------------------------------> " << std::endl << iter << std::endl;
    
    for (sample=0 ; sample < n_samples; sample++) {
      // initialize X with uniquely mapping counts:
      X.column(0) = list_X_unique[sample].column(0);
      X.column(1) = list_X_unique[sample].column(1);
      X.column(2) = list_X_unique[sample].column(2);
      
      // REMOVE latent variables sampling:
      pi_SU_sample = PI_SU[sample];
      pi_gene_sample = PI_gene.column(sample);
      List list_EC_gene_id_sample = list_EC_gene_id[sample];
      List list_EC_SU_id_sample = list_EC_SU_id[sample];
      
      for (ec=0 ; ec < list_EC_gene_id_sample.size(); ec++) {
        Rcpp::IntegerVector ec_index = list_EC_gene_id_sample[ec];
        EC_len = ec_index.length();
        
        Rcpp::NumericVector pi_j(EC_len);
        Rcpp::IntegerVector y_tmp(EC_len); //, index(EC_len);
        
        Rcpp::IntegerVector ec_SU = list_EC_SU_id_sample[ec];
        
        for (i=0 ; i < EC_len; i++) {
          id  = ec_index[i]; // gene id
          pi_j[i] = pi_gene_sample[ id ] * pi_SU_sample( id, ec_SU[i] );
          // replace pi_gene by X_gene as in BANDITS (maybe more stable, simplifies convergence?)
          // pi_j[i] = X_gene[ id ] * pi_SU_sample( id, ec_SU[i] );
        }
        prob_tot = std::accumulate(pi_j.begin(), pi_j.end(), 0.0);
        
        if( prob_tot > 0.0){
          for (i=0 ; i < EC_len; i++) {
            pi_j[i] /= prob_tot;
          }
          
          rmultinom( counts[sample][ ec ], pi_j.begin(), EC_len, y_tmp.begin());
          
          for (i=0 ; i < EC_len; i++) {
            X( ec_index[i], ec_SU[i] ) += y_tmp[i];
          }
        }
      }
      X_gene = Rcpp::rowSums(X);
      
      for (gene_id=0 ; gene_id < n_genes; gene_id++) {
        original_delta = delta_SU[ numeric_groups[sample] ]( gene_id ,_);
        
        // sample pi_S (prob SUA read):
        dir_sample_SUA_sum = 0.0;
        for (unsigned int k = 0; k < 3; ++k) {
          pi_SU_tmp[k] = as<double>(Rcpp::rgamma(1, X(gene_id, k) + original_delta[k], 1));
          dir_sample_SUA_sum += pi_SU_tmp[k];
        }
        
        cond_SUA = true;
        for (unsigned int k = 0; k < 3; ++k) {
          pi_SU_tmp[k] = pi_SU_tmp[k]/dir_sample_SUA_sum;
          
          if( ( Rcpp::NumericVector::is_na(pi_SU_tmp[k]) ) || ( pi_SU_tmp[k] < pow(10, -100) )  ) {
            cond_SUA = false;
          }
        }
        
        // Update pi_SU_tmp
        if(cond_SUA){
          pi_SU_sample.row(gene_id) = pi_SU_tmp;
        }
        
        // sample pi_gene (prob gene):
        pi_gene_sample[gene_id] = as<double>(Rcpp::rgamma(1, X_gene[gene_id] + 1, 1));
        // normalizing pi_gene_sample is actually not necessary: pi_j is later normalized anyway...unless it leads to underflow or overflow problems.
      }
      
      PI_SU[sample] = pi_SU_sample;
    }
    
    double ll_prop, ll_original, ll = 0;
    // hyper-parameters for PI_SU Metropolis sampling:
    for (gene_id=0 ; gene_id < n_genes_keep; gene_id++) {
      gene_to_keep_id = keep_genes_id[gene_id];
      
      for (gr=0 ; gr < n_groups; gr++) {
        // this requires 2 samples per group (to create a matrix):
        // gene SU pi:
        NumericMatrix PI_SU_gene(n_samples_per_group[gr], 3);
        for(sample = 0; sample < n_samples_per_group[ gr ]; sample ++){
          id = sample_ids_per_group[gr][sample];
          PI_SU_gene.row(sample) = PI_SU[ id ].row(gene_to_keep_id);
          // put group first -> select sub-matrix of PI_SU
        }
        
        original_delta = delta_SU[gr](gene_to_keep_id,_);
        
        if( true ){ //}iter < 100 ){ // simple R.W. for the initial 200 iter.
          prop_delta_SU = Rcpp::rnorm(3, 0, 0.2);
          prop_delta_SU += log(original_delta);
        }
        else{ // ARW for the following iterations:
          if( (iter == 100) || (iter == burn_in) ){ // update covar matrix after 200 iter and when reaching the burn-in:
            Y.column(0) = MCMC_bar_pi_1[gr].column(gene_id);
            Y.column(1) = MCMC_bar_pi_2[gr].column(gene_id);
            Y.column(2) = MCMC_bar_pi_3[gr].column(gene_id);
            
            if( iter == burn_in ){ //subset to consider results until "iter-1"
              Y = Y(Range(100, iter-1), _); // for burn-in, also remove initial 100 iterations
              df = iter - 101;
            }else{
              Y = Y(Range(50, iter-1), _);
              df = iter - 51;
            }
            
            covRcpp_bis(Y, cv_alpha, c_diag, c_prop, df);
            
            tmp_chol = R_chol( cv_alpha );
            chol[gr][gene_id] = Rcpp::as<NumericMatrix>(tmp_chol);
          }
          
          // modify here to compute the cholesky decomposition above.
          my_rmvnorm_final(prop_delta_SU, Rcpp::as<arma::mat>( chol[gr][gene_id] ), log(original_delta), sample_arma);
        }
        
        exp_prop_delta_SU = exp(prop_delta_SU);
        // condition for proposed values:
        cond = true;
        
        // check if there are NA's in the proposed values:
        if( NumericVector::is_na(prop_delta_SU[0]) || NumericVector::is_na(prop_delta_SU[1]) || NumericVector::is_na(prop_delta_SU[2]) ){
          cond = false;
          Rcout << "1st control failed: " << std::endl << prop_delta_SU << std::endl;
        }
        
        // check that:
        // 1) All exp(alpha_prop) are strictly > 0; 2) the sum of exp(alpha_prop) is < 10^10
        if(cond){
          if( !is_true( all( exp_prop_delta_SU > 0 ) ) || (sum(exp_prop_delta_SU) >  pow(10, 9) ) ){
            cond = false;
            Rcout << "2nd control failed: " << std::endl << exp_prop_delta_SU << std::endl;
          }
        }
        
        if(cond){
          ll_prop = ll_alpha(PI_SU_gene, exp_prop_delta_SU, n_samples_per_group[ gr ]);
          
          ll_original = ll_alpha(PI_SU_gene, original_delta, n_samples_per_group[ gr ]);
          
          if( prior_TF ){ // if prior has been specified
            prior_informative(prior_prop_delta_SU, prop_delta_SU, mean_log_delta, sd_log_delta, sd_prior_non_informative, arma_mat_prior);
          }else{ // if prior has NOT been specified: TO DO: compute prior once only at the beginning!
            prior_prop_delta_SU = R::dnorm(prop_delta_SU[0], 0, sd_prior_non_informative, true) + R::dnorm(prop_delta_SU[1], 0, sd_prior_non_informative, true) + R::dnorm(prop_delta_SU[2], 0, sd_prior_non_informative, true);
          }
          
          // ACCEPT/REJECT proposal for group A:
          alpha = std::min( ll_prop - ll_original + prior_prop_delta_SU - prior_delta_SU(gene_id, gr), 0.0);
          
          if( !NumericVector::is_na(alpha) ){ // check no NA's in alpha.
            cond_01 = Rcpp::rbinom(1, 1, exp(alpha))(0);
            
            if(cond_01 == 1){ // # Update parameter and prior
              ll_original = ll_prop;
              original_delta = clone(exp_prop_delta_SU);
              delta_SU[gr](gene_to_keep_id,_) = original_delta;
              prior_delta_SU(gene_id, gr) = prior_prop_delta_SU;
            }
          }
        }
        
        // only until burn_in, store values for both alpha and beta, used to compute the 
        MCMC_bar_pi_1[gr](iter,gene_id) = original_delta[0];
        MCMC_bar_pi_2[gr](iter,gene_id) = original_delta[1];
        MCMC_bar_pi_3[gr](iter,gene_id) = original_delta[2];
      }
      ll += ll_original;
    }
    MCMC_ll[iter] = ll + sum(prior_delta_SU);
  }
  
  return MCMC_ll;
}
