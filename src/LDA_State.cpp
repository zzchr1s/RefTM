#include<vector>
#include<unordered_set>
#include<cmath>
#include<limits.h>
#include<random>

#include<RcppArmadillo.h>


#include "lda_svi.h"

using namespace std;
using namespace Rcpp;

LDA_State::LDA_State(int n_docs,int vocab_size,int n_topics,int n_topics_ref,std::unordered_map<int,std::unordered_map<int,int>> &data,double eta_val,double alpha_val, Rcpp::NumericMatrix refLambda){
	D=n_docs;
	V=vocab_size;
	K=n_topics;
                K0 = n_topics_ref;
	dtm = data;
	refLambda= refLambda;

	t=0;

	eta = eta_val;
	alpha = alpha_val;

	gamma = arma::mat(D,K);
	lambda = arma::mat(K,V); 

	Elogbeta = arma::mat(K,V);
	expElogbeta = arma::mat(K,V);
	
	for (auto &elem:gamma){
		elem = 0;
	}
	for (auto &elem:lambda){
		elem = R::rgamma(100,0.01);
	}

	for (int i=0;i<K;i++){
		Elogbeta.row(i) = dirichlet_expectation(lambda.row(i));
		expElogbeta.row(i) = arma::exp(Elogbeta.row(i));
	}
               
               if (K0 !=0){
                          arma::mat tmp = Rcpp::as<arma::mat>(refLambda);
                          for (int i=0;i<K0;i++){
		Elogbeta.row(i) = dirichlet_expectation(tmp.row(i));
		expElogbeta.row(i) = arma::exp(Elogbeta.row(i));
	          }
               }
                
}

void LDA_State::update_minibatch(std::vector<int> documents,int maxiter,double tau_0,double kappa, Rcpp::NumericMatrix refLambda){

	t++;

	int batchsize = documents.size();
	double rho_t = std::pow((tau_0+t),-kappa);


	arma::mat sufficient_statistics(K,V);

	unordered_map<int,int> indices;

	for (int i=0;i<documents.size();i++){
		indices.insert(make_pair(documents[i],i));
	}

	for (auto doc_id : documents){
		for (double &elem : gamma.row(doc_id)){
			elem = R::rgamma(100,0.01);
		}
	}




	arma::mat Elogtheta(batchsize,K);
	arma::mat expElogtheta(batchsize,K);

	for (int i=0;i<batchsize;i++){
		Elogtheta.row(i) = dirichlet_expectation(gamma.row(documents[i]));
		expElogtheta.row(i) = arma::exp(Elogtheta.row(i));
	}

	
	double mean_abs_change=0;	
	arma::mat sufficient_stats = arma::zeros(K,V);



	for (auto doc_id : documents){
		std::vector<int> word_ids;
		std::vector<int> word_cts;
		

		unordered_map<int,int> word_count_map  = dtm[doc_id];
		for (pair<int,int> index : word_count_map){
			word_ids.push_back(index.first);
			word_cts.push_back(index.second);
		}
		
    
		
		arma::rowvec word_cts_vec(word_cts.size());
		for (int i=0;i<word_cts.size();i++){
		  word_cts_vec[i] = static_cast<double>(word_cts[i]);
		}
		
		
	  arma::rowvec gamma_d(K);
	  for (double &elem:gamma_d){
		elem = R::rgamma(100,0.01);
	  }
	  arma::rowvec Elogtheta_d = Elogtheta.row(indices[doc_id]);
    	  arma::rowvec expElogtheta_d = expElogtheta.row(indices[doc_id]);
	
	  arma::uvec word_indices(word_ids.size());
	  for (int i=0;i<word_ids.size();i++){
	    word_indices[i] = static_cast<unsigned int>(word_ids[i]);
	  }
	  
	  arma::mat expElogbeta_d = expElogbeta.cols(word_indices);

	  arma::rowvec phinorm = expElogtheta_d * expElogbeta_d + 1e-100;

		for(int i=0;i<maxiter;i++){
		  

		  arma::rowvec gamma_d_old = gamma_d;

		  gamma_d = alpha + (expElogtheta_d % ((word_cts_vec/phinorm) * expElogbeta_d.t()));

		  Elogtheta_d = dirichlet_expectation(gamma_d);

		  expElogtheta_d = arma::exp(Elogtheta_d);

		  phinorm = expElogtheta_d * expElogbeta_d + 1e-100;

		  mean_abs_change = arma::mean(arma::abs(gamma_d - gamma_d_old));
		  if (i==i && R::runif(0,1) < 0.001){
			Rcout << mean_abs_change << endl;
		  }

		}	
		
		

	gamma.row(doc_id) = gamma_d;

	sufficient_statistics.cols(word_indices) += expElogtheta_d.t()*(word_cts_vec/phinorm);

	}
	
	sufficient_statistics = sufficient_statistics % expElogbeta;
	
	//update word-topic matrix lambda
	
	lambda = lambda * (1-rho_t) + rho_t * (eta + D*sufficient_statistics/batchsize);
                
                if (K0 !=0){
                          arma::mat tmp = Rcpp::as<arma::mat>(refLambda);
                          for (int i=0;i<K0;i++){
		lambda.row(i) = tmp.row(i);
	          }
               }

	for (int i = 0 ;i< K ;i++){
		Elogbeta.row(i) = dirichlet_expectation(lambda.row(i));
		expElogbeta.row(i) = arma::exp(Elogbeta.row(i));
	}
	

}

void LDA_State::fit_model(int passes,int batchsize,int maxiter,double tau_0,double kappa, Rcpp::NumericMatrix refLambda){
	
	for (int i=0;i<passes;i++){

		std::vector<int> docs;
		for (int j=0;j<D;j++){
			docs.push_back(j);
		}

		unsigned seed = floor(R::runif(0,1)*UINT_MAX);
		shuffle(docs.begin(), docs.end(),mt19937_64(seed));//64-bit mersenne twister
	
		int batches=0;
		while (!docs.empty()){
			if (docs.size() >= batchsize){
				std::vector<int> minibatch(docs.end() - batchsize,docs.end());
				docs.erase(docs.end() - batchsize,docs.end());
				update_minibatch(minibatch,maxiter,tau_0,kappa,refLambda);
			}
			else {
				update_minibatch(docs,maxiter,tau_0,kappa,refLambda);
				docs.erase(docs.begin(),docs.end());
			}
			batches++;
		}
	}
}
