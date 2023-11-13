#include <cassert>
#include <ctime>
#include <iostream>
#include <random>
#include <vector>

// [[Rcpp::plugins(cpp17)]]

// Armadillo
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
//using namespace arma;

// "> (btw, in general it's more efficient to access matrices as columns rather than rows)."

void prior_informative_EC_US(double& prior,
                             Rcpp::NumericVector const& x,
                             Rcpp::NumericVector const& prior_log_disp, 
                             Rcpp::NumericVector const& prior_log_S){
  double log_disp = log(sum(exp(x)));
  
  prior = x[1] - log_disp;
  prior += ((R::dnorm(log_disp, prior_log_disp[0], prior_log_disp[1], true) + R::dnorm(x[0], prior_log_S[0], prior_log_S[1], true)));
}

double ll_alpha_EC_US(Rcpp::NumericMatrix const& pi, 
                      Rcpp::NumericVector const& alpha,
                      unsigned int const& N){
  double ll = 0.0;
  double sum_alpha = 0.0;
  
  for (unsigned int k = 0; k < 2; ++k){
    sum_alpha += alpha[k];
    
    ll -= N * R::lgammafn(alpha[k]);
    
    for (unsigned int i = 0; i < N; ++i){
      ll += log(pi(i,k)) * ( alpha[k] - 1);
    }
  }
  ll += N * R::lgammafn(sum_alpha);
  
  if( ISNAN(ll) ){
    //Rcout << "na ll is" << std::endl << ll << std::endl;
    ll = -std::numeric_limits<double>::infinity();
    //Rcout << "ll becomes" << std::endl << ll << std::endl;
  }
  
  return ll;
}


void covRcpp_bis_EC_US(Rcpp::NumericMatrix& Y,
                       Rcpp::NumericMatrix& cov,
                       double const& c_diag,
                       double const& c_prop, 
                       unsigned int const& df){
  // Centering the matrix:
  Y(_, 0) = Y(_, 0) - mean(Y(_, 0));
  Y(_, 1) = Y(_, 1) - mean(Y(_, 1));
  // COV only updates the bottom right corner element!
  
  // Computing the covariance matrix
  for (unsigned int i = 0; i < 2; ++i){
    for (unsigned int j = 0; j <= i; ++j){ // I REMOVE THE FIRST 100 ROWS OF THE MCMC AND START FROM THE 101-st row.
      cov(i,j) = c_prop * Rcpp::sum(Y(Rcpp::_, i)*Y(Rcpp::_, j))/df;
      cov(j,i) = cov(i,j);
    }
    cov(i,i) += c_diag;
  }
}

void my_rmvnorm_final_EC_US(Rcpp::NumericVector& alpha_prop,
                            arma::mat const& chol_R, 
                            Rcpp::NumericVector const& mean,
                            arma::vec& sample_arma){
  //sample_arma = arma::randn(2);
  sample_arma = Rcpp::as<arma::vec>(Rcpp::rnorm(2));
  sample_arma = trans(sample_arma.t() * chol_R);
  
  alpha_prop = Rcpp::wrap(sample_arma);
  alpha_prop += mean;
}


// [[Rcpp::export]]
Rcpp::NumericVector Rcpp_MCMC_EC_US( unsigned int const& n_samples, // N samples
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
                                     Rcpp::NumericMatrix const& list_X_unique, // SU uniquely mapping counts
                                     Rcpp::ListOf<Rcpp::ListOf<Rcpp::IntegerVector>> const& list_EC_gene_id, // SU uniquely mapping counts
                                     Rcpp::ListOf<Rcpp::IntegerVector> const& counts, // EC counts (integers)
                                     Rcpp::ListOf<Rcpp::NumericMatrix>& MCMC_bar_pi_1,
                                     Rcpp::ListOf<Rcpp::NumericMatrix>& MCMC_bar_pi_2,
                                     Rcpp::ListOf<Rcpp::NumericMatrix>& delta_SU,
                                     Rcpp::NumericVector const& prior_log_disp, 
                                     Rcpp::NumericVector const& prior_log_S, 
                                     Rcpp::IntegerVector const& sample_EC,
                                     Rcpp::IntegerVector const& sample_SU_TF,
                                     Rcpp::NumericVector const& eff_len_S,
                                     Rcpp::NumericVector const& eff_len_U,
                                     double const& c_prop){
  // List containing arma ARW matrices:
  List chol_1(n_genes_keep);
  List chol_2(n_genes_keep);
  
  // Obtain environment containing chol function:
  //Rcpp::Environment base("package:base"); 
  // Make chol function callable from C++
  //Rcpp::Function R_chol = base["chol"];    
  
  for (unsigned int gr=0 ; gr < n_groups; gr++){
    Rcpp::NumericMatrix tmp_1(N_MCMC, n_genes_keep), tmp_2(N_MCMC, n_genes_keep);
    MCMC_bar_pi_1[gr] = tmp_1;
    MCMC_bar_pi_2[gr] = tmp_2;
  }
  double tmp_1, tmp_2, tmp;
  
  // X is a temporary matrix to update the (s,u, a) reads of each sample:
  Rcpp::NumericMatrix X_all_samples(2*n_genes, n_samples), prior_delta_SU(n_genes_keep, n_groups), cv_alpha(2, 2), pi_SU_sample(n_genes, 2), Y(N_MCMC, 2), Y_sel;
  
  // keep a temporary vectorized form of pi_S to speed-up memory access
  Rcpp::NumericVector X(2*n_genes), pi_gene_SU_sample(2*n_genes), prop_delta_SU(2), original_delta(2), exp_prop_delta_SU(2), MCMC_ll(N_MCMC), pi_SU_tmp(2);
  
  unsigned int EC_len, iter, sample, ec, i, gene_id, gr, cond_01, id, df, gene_to_keep_id;
  double prob_tot, dir_sample_SUA_sum, prior_prop_delta_SU, alpha;
  //ll_delta_diff
  
  bool cond, cond_SUA;
  
  // constant to add to the diagonal of the ARW matrix.
  double c_diag = 0.001;// the proportionality constant for the ARW matrix.
  //double c_prop = 0.3;
  
  arma::vec sample_arma(2); // arma for ARW proposal
  
  // initalize prior for delta_SU:
  for (gr=0 ; gr < n_groups; gr++){
    for (gene_id=0 ; gene_id < n_genes_keep; gene_id++){
      original_delta = delta_SU[gr]( gene_id,_);
      
      // Initialize the prior:
      //prior_prop_delta_SU = R::dnorm(log(original_delta[0]), mean_log_delta[0], sd_log_delta[0], true) + R::dnorm(log(original_delta[1]), mean_log_delta[1], sd_log_delta[1], true);
      prior_informative_EC_US(prior_prop_delta_SU, log(original_delta), prior_log_disp, prior_log_S);
      prior_delta_SU(gene_id, gr) = prior_prop_delta_SU;
    }
  }
  
  for (iter=0 ; iter < N_MCMC; iter++){
    
    //Rcout << "iter" << std::endl << iter << std::endl;
    
    Rcpp::checkUserInterrupt(); // interrupt code from R (benchmark cost).
    
    //Rcout << "starting iter: ------------------------------------> " << std::endl << iter << std::endl;
    
    // sample loop STARTS; 3 Gibbs samplers: latent variables, pi_gene, and pi_SU 
    for (sample=0 ; sample < n_samples; sample++){
      
      pi_SU_sample = PI_SU[sample];
      
      // sample ECs every XX iterations:
      if( sample_EC[iter] == 1 ){
        //Rcout << "EC iter" << std::endl << iter << std::endl;
        
        // initialize X with uniquely mapping counts:
        for (gene_id=0 ; gene_id < 2*n_genes; gene_id++){
          X[gene_id] = list_X_unique(gene_id, sample);
        }
        
        // pi_gene_SU_sample = pi_gene * pi_US/eff_len
        pi_gene_SU_sample = PI_gene.column(sample);
        
        List list_EC_gene_id_sample = list_EC_gene_id[sample];
        
        // latent variable sampling:
        for (ec=0 ; ec < list_EC_gene_id_sample.size(); ec++){
          Rcpp::IntegerVector ec_index = list_EC_gene_id_sample[ec];
          EC_len = ec_index.length();
          
          Rcpp::NumericVector pi_j(EC_len);
          Rcpp::IntegerVector y_tmp(EC_len); //, index(EC_len);
          
          // TODO: try to get elements in a vector form: pi_gene_SU_sample[ ec_index ] ?
          for (i=0 ; i < EC_len; i++){
            //id  = ec_index[i]; // gene id
            pi_j[i] = pi_gene_SU_sample[ ec_index[i] ];
          }
          prob_tot = sum(pi_j); 
          //std::accumulate(pi_j.begin(), pi_j.end(), 0.0);
          
          if( prob_tot > 0.0){
            for (i=0 ; i < EC_len; i++){
              pi_j[i] /= prob_tot;
            }
            
            rmultinom( counts[sample][ ec ], pi_j.begin(), EC_len, y_tmp.begin());
            
            for (i=0 ; i < EC_len; i++){
              //id  = ec_index[i]; // gene id
              X[ ec_index[i] ] += y_tmp[i];
            }
          }
        }
        // store the latent counts, of the sample-th sample, in a matrix:
        X_all_samples.column(sample) = clone(X);
      }else{
        // take the latent counts, of the sample-th sample, from the matrix:
        X = X_all_samples.column(sample);
      }
      
      // pi_SU only for selected genes (n_genes_keep):
      for (gene_id=0 ; gene_id < n_genes_keep; gene_id++){
        gene_to_keep_id = keep_genes_id[gene_id];
        
        original_delta = delta_SU[ numeric_groups[sample] ]( gene_id ,_);
        
        pi_SU_tmp[0] = as<double>(Rcpp::rgamma(1, X[gene_to_keep_id] + original_delta[0], 1));
        pi_SU_tmp[1] = as<double>(Rcpp::rgamma(1, X[n_genes + gene_to_keep_id] + original_delta[1], 1));
        dir_sample_SUA_sum = pi_SU_tmp[0] + pi_SU_tmp[1];
        
        pi_SU_tmp[0] /= dir_sample_SUA_sum;
        pi_SU_tmp[1] /= dir_sample_SUA_sum;
        
        cond_SUA = true;
        if( ( Rcpp::NumericVector::is_na(pi_SU_tmp[0]) ) || ( pi_SU_tmp[0] < pow(10, -100) ) ||  ( Rcpp::NumericVector::is_na(pi_SU_tmp[1]) ) || ( pi_SU_tmp[1] < pow(10, -100) )){
          cond_SUA = false;
        }
        
        if(cond_SUA){
          pi_SU_sample.row(gene_to_keep_id) = pi_SU_tmp;
        }
      }
      PI_SU[sample] = clone(pi_SU_sample);
      
      // sample pi_gene every XX iterations:
      if( sample_EC[iter + 1] == 1 ){
        for (gene_id=0 ; gene_id < n_genes; gene_id++){
          
          // pi_SU for NON-selected genes:
          if( sample_SU_TF[gene_id] == 1 ){
            pi_SU_tmp[0] = as<double>(Rcpp::rgamma(1, X[gene_id] + 1, 1));
            pi_SU_tmp[1] = as<double>(Rcpp::rgamma(1, X[n_genes + gene_id] + 1, 1));
            dir_sample_SUA_sum = pi_SU_tmp[0] + pi_SU_tmp[1];
            
            pi_SU_tmp[0] /= dir_sample_SUA_sum;
            pi_SU_tmp[1] /= dir_sample_SUA_sum;
            
            cond_SUA = true;
            if( ( Rcpp::NumericVector::is_na(pi_SU_tmp[0]) ) || ( pi_SU_tmp[0] < pow(10, -100) ) ||  ( Rcpp::NumericVector::is_na(pi_SU_tmp[1]) ) || ( pi_SU_tmp[1] < pow(10, -100) )){
              cond_SUA = false;
            }
            
            if(cond_SUA){
              pi_SU_sample.row(gene_id) = pi_SU_tmp;
            }
          }
          
          // pi_gene for all genes (n_genes):
          // sample pi_gene (prob gene):
          id = n_genes + gene_id;
          tmp = as<double>(Rcpp::rgamma(1, X[gene_id] + X[id] + 1, 1));
          // normalizing pi_gene_sample is actually not necessary: pi_j is later normalized anyway...unless it leads to underflow or overflow problems.
          
          // IMPORTANT: UPDATE pi_gene_SU_sample BEFORE EC computation!
          // otherwise ECs computed on old parameters
          pi_gene_SU_sample[gene_id] = tmp * pi_SU_sample( gene_id, 0 )/eff_len_S[gene_id];
          pi_gene_SU_sample[id] = tmp * pi_SU_sample( gene_id, 1 )/eff_len_U[gene_id];
        }
        
        PI_gene.column(sample) = clone(pi_gene_SU_sample);
      }
    }
    // sample loop ENDS.
    
    // hyper-parameters (delta_SU) Metropolis sampling:
    double ll_prop, ll_original, ll = 0;
    
    // prepare Y_sel matrix for ARW matrix calculation:
    // update matrix multiple times before burn-in:
    if( iter == 200 ){
      df = iter - 101;
      Rcpp::NumericMatrix Y_sel(df, 2);
    }
    if( iter == 300 ){
      df = iter - 201;
      Rcpp::NumericMatrix Y_sel(df, 2);
    }
    
    for (gene_id=0 ; gene_id < n_genes_keep; gene_id++){
      gene_to_keep_id = keep_genes_id[gene_id];
      
      for (gr=0 ; gr < n_groups; gr++){
        // this requires 2 samples per group (to create a matrix):
        // gene SU pi:
        NumericMatrix PI_SU_gene(n_samples_per_group[gr], 2);
        for(sample = 0; sample < n_samples_per_group[ gr ]; sample ++){
          id = sample_ids_per_group[gr][sample];
          PI_SU_gene.row(sample) = PI_SU[ id ].row(gene_to_keep_id);
          // put group first -> select sub-matrix of PI_SU
        }
        
        original_delta = delta_SU[gr](gene_id,_);
        
        if( iter < 200 ){ // simple R.W. for the initial 200 iter.
          prop_delta_SU = Rcpp::rnorm(2, 0, 0.2);
          prop_delta_SU += log(original_delta);
        }
        else{ // ARW for the following iterations:
          // update ARW matrices at 100, 200, 300, 400, 500 and burn_in iterations:
          if( (iter == 200) || (iter == 300) || (iter == 400) || (iter == 500) || (iter == burn_in) ){
            // update covar matrix after 200 iter and when reaching the burn-in:
            Y.column(0) = log(MCMC_bar_pi_1[gr].column(gene_id));
            Y.column(1) = log(MCMC_bar_pi_2[gr].column(gene_id));
            
            if( iter == 200 ){
              Y_sel = Y(Range(100, iter-1), _); // for burn-in, also remove initial 100 iterations
            }else{
              Y_sel = Y(Range(200, iter-1), _); // for burn-in, also remove initial 100 iterations
            }
            
            covRcpp_bis_EC_US(Y_sel, cv_alpha, c_diag, c_prop, df);
            
            // in BANDITS we used arma::chol
            if(gr == 0){
              chol_1[gene_id] = arma::chol(Rcpp::as<arma::mat>(cv_alpha));
              //R_chol( cv_alpha );
            }else{
              chol_2[gene_id] = arma::chol(Rcpp::as<arma::mat>(cv_alpha));
              //R_chol( cv_alpha );
            }
          }
          
          if(gr == 0){
            my_rmvnorm_final_EC_US(prop_delta_SU, chol_1[gene_id], log(original_delta), sample_arma);
          }else{
            my_rmvnorm_final_EC_US(prop_delta_SU, chol_2[gene_id], log(original_delta), sample_arma);
          }
        }
        
        exp_prop_delta_SU = exp(prop_delta_SU);
        // condition for proposed values:
        cond = true;
        
        // check if there are NA's in the proposed values:
        if( NumericVector::is_na(prop_delta_SU[0]) || NumericVector::is_na(prop_delta_SU[1]) ){
          cond = false;
          //Rcout << "1st control failed: " << std::endl << prop_delta_SU << std::endl;
        }
        
        // check that:
        // 1) All exp(alpha_prop) are strictly > 0; 2) the sum of exp(alpha_prop) is < 10^10
        if(cond){
          if( !is_true( all( exp_prop_delta_SU > 0 ) ) || (sum(exp_prop_delta_SU) >  pow(10, 9) ) ){
            cond = false;
            //Rcout << "2nd control failed: " << std::endl << exp_prop_delta_SU << std::endl;
          }
        }
        
        if(cond){
          ll_prop = ll_alpha_EC_US(PI_SU_gene, exp_prop_delta_SU, n_samples_per_group[ gr ]);
          ll_original = ll_alpha_EC_US(PI_SU_gene, original_delta, n_samples_per_group[ gr ]);
          
          //prior_prop_delta_SU = R::dnorm(prop_delta_SU[0], mean_log_delta[0], sd_log_delta[0], true) + R::dnorm(prop_delta_SU[1], mean_log_delta[1], sd_log_delta[1], true);
          prior_informative_EC_US(prior_prop_delta_SU, prop_delta_SU, prior_log_disp, prior_log_S);
          
          // ACCEPT/REJECT proposal for group A:
          alpha = std::min( ll_prop - ll_original + prior_prop_delta_SU - prior_delta_SU(gene_id, gr), 0.0);
          
          if( !NumericVector::is_na(alpha) ){ // check no NA's in alpha.
            cond_01 = Rcpp::rbinom(1, 1, exp(alpha))(0);
            
            if(cond_01 == 1){ // # Update parameter and prior
              ll_original = ll_prop;
              original_delta = clone(exp_prop_delta_SU);
              delta_SU[gr](gene_id,_) = original_delta;
              prior_delta_SU(gene_id, gr) = prior_prop_delta_SU;
            }
          }
        }
        
        // only until burn_in, store values for both alpha and beta, used to compute the 
        MCMC_bar_pi_1[gr](iter,gene_id) = original_delta[0];
        MCMC_bar_pi_2[gr](iter,gene_id) = original_delta[1];
        
        ll += ll_original;
      }
    }
    MCMC_ll[iter] = ll + sum(prior_delta_SU);
  }
  
  // create 2 tmp vectors, to speed-up access to the matrix:
  Rcpp::NumericVector tmp_vec_1(N_MCMC), tmp_vec_2(N_MCMC);
  for (gene_id=0 ; gene_id < n_genes_keep; gene_id++){
    gene_to_keep_id = keep_genes_id[gene_id];
    for (unsigned int gr=0 ; gr < n_groups; gr++){
      tmp_vec_1 = MCMC_bar_pi_1[gr].column(gene_id);
      tmp_vec_2 = MCMC_bar_pi_2[gr].column(gene_id);
      for (iter=0 ; iter < N_MCMC; iter++){
        original_delta[0] = tmp_vec_1[iter];
        original_delta[1] = tmp_vec_2[iter];
        
        // divide by eff length:
        //original_delta[0] = original_delta[0]/eff_len_S[gene_to_keep_id];
        //original_delta[1] = original_delta[1]/eff_len_U[gene_to_keep_id];
        // rescale so that pi's add to 1:
        tmp = original_delta[0] + original_delta[1];
        
        // store normalized values:
        tmp_vec_1[iter] = original_delta[0]/tmp;
        tmp_vec_2[iter] = original_delta[1]/tmp;
      }
      MCMC_bar_pi_1[gr].column(gene_id) = tmp_vec_1;
      MCMC_bar_pi_2[gr].column(gene_id) = tmp_vec_2;
    }
  }
  
  return MCMC_ll;
}