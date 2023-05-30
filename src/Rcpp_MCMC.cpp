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
                       Rcpp::NumericVector const& mean, 
                       Rcpp::NumericVector const& sd){
  double log_disp = log(sum(exp(x)));
  // log-determinant:
  prior = x[2] - log_disp;
  prior += R::dnorm(log_disp, mean[0], sd[0], true) + R::dnorm(x[0], mean[1], sd[1], true) + R::dnorm(x[1], mean[2], sd[2], true);
}

double
  ll_alpha(Rcpp::NumericMatrix const& pi, 
           Rcpp::NumericVector const& alpha,
           unsigned int const& N)
  {
    double ll = 0.0;
    double sum_alpha = 0.0;
    
    for (unsigned int k = 0; k < 3; ++k){
      sum_alpha += alpha[k];
      
      ll -= N * R::lgammafn(alpha[k]);
      
      for (unsigned int i = 0; i < N; ++i){
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
                 unsigned int const& df){
  // Centering the matrix:
  Y(_, 0) = Y(_, 0) - mean(Y(_, 0));
  Y(_, 1) = Y(_, 1) - mean(Y(_, 1));
  Y(_, 2) = Y(_, 2) - mean(Y(_, 2));
  
  // COV only updates the bottom right corner element!
  
  // Computing the covariance matrix
  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j <= i; ++j){ // I REMOVE THE FIRST 100 ROWS OF THE MCMC AND START FROM THE 101-st row.
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
  //sample_arma = arma::randn(3); 
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
                               Rcpp::NumericMatrix& PI_gene_times_SU, // prob of each gene (for every sample)
                               Rcpp::ListOf<Rcpp::NumericMatrix>& PI_SU, // prob of each gene (for every sample)
                               Rcpp::NumericMatrix const& list_X_unique, // SU uniquely mapping counts
                               Rcpp::ListOf<Rcpp::ListOf<Rcpp::IntegerVector>> const& list_EC_gene_id, // SU uniquely mapping counts
                               Rcpp::ListOf<Rcpp::IntegerVector> const& counts, // EC counts (integers)
                               Rcpp::ListOf<Rcpp::NumericMatrix>& MCMC_bar_pi_1,
                               Rcpp::ListOf<Rcpp::NumericMatrix>& MCMC_bar_pi_2,
                               Rcpp::ListOf<Rcpp::NumericMatrix>& MCMC_bar_pi_3,
                               Rcpp::ListOf<Rcpp::NumericMatrix>& delta_SU,
                               Rcpp::NumericVector const& mean_log_delta, 
                               Rcpp::NumericVector const& sd_log_delta,
                               Rcpp::IntegerVector const& sample_EC,
                               Rcpp::ListOf<Rcpp::NumericMatrix>& X_list,
                               Rcpp::IntegerVector const& sample_SU_TF,
                               double const& c_prop){
  // List containing arma ARW matrices:
  List chol_1(n_genes_keep);
  List chol_2(n_genes_keep);
  
  // Obtain environment containing chol function:
  //Rcpp::Environment base("package:base"); 
  // Make chol function callable from C++precision$prior
  //Rcpp::Function R_chol = base["chol"];    
  
  for (unsigned int gr=0 ; gr < n_groups; gr++){
    Rcpp::NumericMatrix tmp_1(N_MCMC, n_genes_keep), tmp_2(N_MCMC, n_genes_keep), tmp_3(N_MCMC, n_genes_keep);
    MCMC_bar_pi_1[gr] = tmp_1;
    MCMC_bar_pi_2[gr] = tmp_2;
    MCMC_bar_pi_3[gr] = tmp_3;
  }
  double tmp_1, tmp_2, tmp_3;
  
  // X is a temporary matrix to update the (s,u, a) reads of each sample:
  Rcpp::NumericMatrix X_all_samples(3*n_genes, n_samples), prior_delta_SU(n_genes_keep, n_groups), cv_alpha(3, 3), pi_SU_sample(n_genes, 3), Y(N_MCMC, 3), Y_sel;
  
  // keep a temporary vectorized form of pi_S to speed-up memory access
  Rcpp::NumericVector X(3*n_genes), pi_gene_SU_sample(3*n_genes), prop_delta_SU(3), original_delta(3), exp_prop_delta_SU(3), MCMC_ll(N_MCMC), pi_SU_tmp(3);
  
  unsigned int EC_len, iter, sample, ec, i, gene_id, gr, cond_01, id, df, gene_to_keep_id, id_U, id_A;
  double prob_tot, dir_sample_SUA_sum, prior_prop_delta_SU, alpha, tmp, param;
  //ll_delta_diff
  
  bool cond, cond_SUA;
  
  // constant to add to the diagonal of the ARW matrix.
  double c_diag = 0.001;// the proportionality constant for the ARW matrix.
  
  arma::vec sample_arma(3); // arma for ARW proposal
  arma::mat arma_mat_prior(3,3);
  
  // initalize prior for delta_SU:
  for (gr=0 ; gr < n_groups; gr++){
    for (gene_id=0 ; gene_id < n_genes_keep; gene_id++){
      original_delta = delta_SU[gr]( gene_id,_);
      
      // Initialize the prior:
      prior_informative(prior_prop_delta_SU, log(original_delta),  mean_log_delta, sd_log_delta);
      prior_delta_SU(gene_id, gr) = prior_prop_delta_SU;
    }
  }
  
  for (iter=0 ; iter < N_MCMC; iter++){
    
    //Rcout << "iter" << std::endl << iter << std::endl;
    
    Rcpp::checkUserInterrupt(); // interrupt code from R (benchmark cost).
    
    //Rcout << "starting iter: ------------------------------------> " << std::endl << iter << std::endl;
    
    for (sample=0 ; sample < n_samples; sample++){
      pi_SU_sample = PI_SU[sample];
      
      if( sample_EC[iter] == 1 ){
        // initialize X with uniquely mapping counts:
        X = list_X_unique.column(sample);
        
        pi_gene_SU_sample = PI_gene_times_SU.column(sample);
        List list_EC_gene_id_sample = list_EC_gene_id[sample];
        
        for (ec=0 ; ec < list_EC_gene_id_sample.size(); ec++){
          Rcpp::IntegerVector ec_index = list_EC_gene_id_sample[ec];
          EC_len = ec_index.length();
          
          Rcpp::NumericVector pi_j(EC_len);
          Rcpp::IntegerVector y_tmp(EC_len); //, index(EC_len);
          
          for (i=0 ; i < EC_len; i++){
            pi_j[i] = pi_gene_SU_sample[ ec_index[i] ];
          }
          prob_tot = std::accumulate(pi_j.begin(), pi_j.end(), 0.0);
          
          if( prob_tot > 0.0){
            for (i=0 ; i < EC_len; i++){
              pi_j[i] /= prob_tot;
            }
            
            rmultinom( counts[sample][ ec ], pi_j.begin(), EC_len, y_tmp.begin());
            
            for (i=0 ; i < EC_len; i++){
              X[ ec_index[i] ] += y_tmp[i];
            }
          }
        }
        X_all_samples.column(sample) = clone(X);
      }else{
        // take the latent counts, of the sample-th sample, from the matrix:
        X = X_all_samples.column(sample);
      }
      
      // pi_SU only for selected genes (n_genes_keep):
      for (gene_id=0 ; gene_id < n_genes_keep; gene_id++){
        gene_to_keep_id = keep_genes_id[gene_id];
        
        original_delta = delta_SU[ numeric_groups[sample] ]( gene_id ,_);
        
        pi_SU_tmp[0] = as<double>(Rcpp::rgamma(1, X(gene_to_keep_id) + original_delta[0], 1));
        pi_SU_tmp[1] = as<double>(Rcpp::rgamma(1, X(n_genes + gene_to_keep_id) + original_delta[1], 1));
        pi_SU_tmp[2] = as<double>(Rcpp::rgamma(1, X(2*n_genes + gene_to_keep_id) + original_delta[2], 1));
        dir_sample_SUA_sum = pi_SU_tmp[0] + pi_SU_tmp[1] + pi_SU_tmp[2];
        
        pi_SU_tmp[0] /= dir_sample_SUA_sum;
        pi_SU_tmp[1] /= dir_sample_SUA_sum;
        pi_SU_tmp[2] /= dir_sample_SUA_sum;
        
        cond_SUA = true;
        for (unsigned int k = 0; k < 3; ++k){
          if( ( Rcpp::NumericVector::is_na(pi_SU_tmp[k]) ) || ( pi_SU_tmp[k] < pow(10, -100) )  ){
            cond_SUA = false;
          }
        }
        
        // Update pi_SU_tmp
        if(cond_SUA){
          pi_SU_sample.row(gene_to_keep_id) = pi_SU_tmp;
        }
      }
      PI_SU[sample] = pi_SU_sample;
      
      if( sample_EC[iter + 1] == 1 ){
        for (gene_id=0 ; gene_id < n_genes; gene_id++){ 
          id_U = gene_id + n_genes;
          id_A = gene_id + 2 * n_genes;
          
          // pi_SU for NON-selected genes (selected were already sampled above):
          //if( sample_SU_TF[gene_id] == 1 ){
          pi_SU_tmp[0] = as<double>(Rcpp::rgamma(1, X[gene_id] + 1, 1));
          pi_SU_tmp[1] = as<double>(Rcpp::rgamma(1, X[id_U] + 1, 1));
          pi_SU_tmp[2] = as<double>(Rcpp::rgamma(1, X[id_A] + 1, 1));
          dir_sample_SUA_sum = pi_SU_tmp[0] + pi_SU_tmp[1] + pi_SU_tmp[2];
          
          pi_SU_tmp[0] /= dir_sample_SUA_sum;
          pi_SU_tmp[1] /= dir_sample_SUA_sum;
          pi_SU_tmp[2] /= dir_sample_SUA_sum;
          
          cond_SUA = true;
          for (unsigned int k = 0; k < 3; ++k){
            if( ( Rcpp::NumericVector::is_na(pi_SU_tmp[k]) ) || ( pi_SU_tmp[k] < pow(10, -100) )  ){
              cond_SUA = false;
            }
          }
          
          if(cond_SUA){
            pi_SU_sample.row(gene_id) = pi_SU_tmp;
          }
          //}
          
          param = X(gene_id) + X(id_U) + X(id_A) + 1;
          tmp = as<double>(Rcpp::rgamma(1, param, 1));
          // normalizing pi_gene (i.e., tmp) is actually not necessary: pi_j is later normalized anyway...unless it leads to underflow or overflow problems.
          
          // IMPORTANT: UPDATE pi_gene_SU_sample BEFORE EC computation!
          // otherwise ECs computed on old parameters
          pi_gene_SU_sample(gene_id) = tmp * pi_SU_sample( gene_id, 0 );
          pi_gene_SU_sample(id_U) = tmp * pi_SU_sample( gene_id, 1 );
          pi_gene_SU_sample(id_A) = tmp * pi_SU_sample( gene_id, 2 );
        }
        PI_gene_times_SU.column(sample) = pi_gene_SU_sample;
      }
    }
    
    // prepare Y_sel matrix for ARW matrix calculation:
    // update matrix multiple times before burn-in:
    if( iter == 200 ){
      df = iter - 101;
      Rcpp::NumericMatrix Y_sel(df, 3);
    }
    if( iter == 300 ){
      df = iter - 201;
      Rcpp::NumericMatrix Y_sel(df, 3);
    }
    
    double ll_prop, ll_original, ll = 0;
    // hyper-parameters for PI_SU Metropolis sampling:
    for (gene_id=0 ; gene_id < n_genes_keep; gene_id++){
      gene_to_keep_id = keep_genes_id[gene_id];
      
      for (gr=0 ; gr < n_groups; gr++){
        // this requires 2 samples per group (to create a matrix):
        // gene SU pi:
        NumericMatrix PI_SU_gene(n_samples_per_group[gr], 3);
        for(sample = 0; sample < n_samples_per_group[ gr ]; sample ++){
          id = sample_ids_per_group[gr][sample];
          PI_SU_gene.row(sample) = PI_SU[ id ].row(gene_to_keep_id);
          // put group first -> select sub-matrix of PI_SU
        }
        
        original_delta = delta_SU[gr](gene_id,_);
        
        if( iter < 200 ){ // simple R.W. for the initial 200 iter.
          prop_delta_SU = Rcpp::rnorm(3, 0, 0.2);
          prop_delta_SU += log(original_delta);
        }
        else{ // ARW for the following iterations:
          // update ARW matrices at 100, 200, 300, 400, 500 and burn_in iterations:
          if( (iter == 200) || (iter == 300) || (iter == 400) || (iter == 500) || (iter == burn_in) ){
            Y.column(0) = log(MCMC_bar_pi_1[gr].column(gene_id));
            Y.column(1) = log(MCMC_bar_pi_2[gr].column(gene_id));
            Y.column(2) = log(MCMC_bar_pi_3[gr].column(gene_id));
            
            if( iter == 200 ){
              Y_sel = Y(Range(100, iter-1), _); // for burn-in, also remove initial 100 iterations
            }else{
              Y_sel = Y(Range(200, iter-1), _); // for burn-in, also remove initial 100 iterations
            }
            
            covRcpp_bis(Y_sel, cv_alpha, c_diag, c_prop, df);
            
            if(gr == 0){
              chol_1[gene_id] = arma::chol(Rcpp::as<arma::mat>(cv_alpha));
              //R_chol( cv_alpha );
            }else{
              chol_2[gene_id] = arma::chol(Rcpp::as<arma::mat>(cv_alpha));
              //R_chol( cv_alpha );
            }          
          }
          
          if(gr == 0){
            my_rmvnorm_final(prop_delta_SU, chol_1[gene_id], log(original_delta), sample_arma);
          }else{
            my_rmvnorm_final(prop_delta_SU, chol_2[gene_id], log(original_delta), sample_arma);
          }
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
          
          prior_informative(prior_prop_delta_SU, prop_delta_SU, mean_log_delta, sd_log_delta);
          
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
        MCMC_bar_pi_3[gr](iter,gene_id) = original_delta[2];
      }
      ll += ll_original;
    }
    MCMC_ll[iter] = ll + sum(prior_delta_SU);
  }
  
  // reparametrize from delta's to pi's
  // create 3 tmp vectors, to speed-up access to the matrix:
  Rcpp::NumericVector tmp_vec_1(N_MCMC), tmp_vec_2(N_MCMC), tmp_vec_3(N_MCMC);
  for (gene_id=0 ; gene_id < n_genes_keep; gene_id++){
    for (unsigned int gr=0 ; gr < n_groups; gr++){
      tmp_vec_1 = MCMC_bar_pi_1[gr].column(gene_id);
      tmp_vec_2 = MCMC_bar_pi_2[gr].column(gene_id);
      tmp_vec_3 = MCMC_bar_pi_3[gr].column(gene_id);
      for (iter=0 ; iter < N_MCMC; iter++){
        original_delta[0] = tmp_vec_1[iter];
        original_delta[1] = tmp_vec_2[iter];
        original_delta[2] = tmp_vec_3[iter];
        
        // rescale so that pi's add to 1:
        tmp = original_delta[0] + original_delta[1] + original_delta[2];
        
        // store normalized values:
        tmp_vec_1[iter] = original_delta[0]/tmp;
        tmp_vec_2[iter] = original_delta[1]/tmp;
        tmp_vec_3[iter] = original_delta[2]/tmp;
      }
      MCMC_bar_pi_1[gr].column(gene_id) = tmp_vec_1;
      MCMC_bar_pi_2[gr].column(gene_id) = tmp_vec_2;
      MCMC_bar_pi_3[gr].column(gene_id) = tmp_vec_3;
    }
  }
  
  return MCMC_ll;
}
