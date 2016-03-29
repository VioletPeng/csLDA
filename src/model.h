/*
 * Copyright (C) 2007 by
 * 
 * 	Xuan-Hieu Phan
 *	hieuxuan@ecei.tohoku.ac.jp or pxhieu@gmail.com
 * 	Graduate School of Information Sciences
 * 	Tohoku University
 *
 * GibbsLDA++ is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 *
 * GibbsLDA++ is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GibbsLDA++; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

/* 
 * References:
 * + The Java code of Gregor Heinrich (gregor@arbylon.net)
 *   http://www.arbylon.net/projects/LdaGibbsSampler.java
 * + "Parameter estimation for text analysis" by Gregor Heinrich
 *   http://www.arbylon.net/publications/text-est.pdf
 */

#ifndef	_MODEL_H
#define	_MODEL_H

#include "constants.h"
#include "dataset.h"
#include "utils.h"

using namespace std;

// LDA model
class model {
public:
    // help read a long line in a file
    char * longLine;
    long max_line_len;
    // fixed options
    string wordmapfile;		// file that contains word map [string -> integer id]
    string trainlogfile;        // training log file
    string alpha_suffix;
    string tau_suffix;
    string tassign_suffix;	// suffix for topic assignment file
    string theta_suffix;	// suffix for theta file
    string pi_suffix;
    string phi0_suffix;		// suffix for phi file
    string phi1_suffix;
    string psi_suffix;
    string phi0bg_suffix;
    string phi1bg_suffix;
    string others_suffix;	// suffix for file containing other parameters
    string twords0_suffix;	// suffix for file containing words-per-topics
    string twords1_suffix;

    string dir;			// model directory
    string dfile;		// data file    
    string model_name;		// model name
    int model_status;		// model status:
				// MODEL_STATUS_UNKNOWN: unknown status
				// MODEL_STATUS_EST: estimating from scratch
				// MODEL_STATUS_ESTC: continue to estimate the model from a previous one
				// MODEL_STATUS_INF: do inference

    dataset * ptrndata;	// pointer to training dataset object
    dataset * pnewdata; // pointer to new dataset object

    mapid2word id2word; // word map [int => string]
    
    // --- model parameters and variables ---    
    int M; // dataset size (i.e., number of docs)
    int V; // vocabulary size
    int K; // number of topics
    int count_cs; //number of codeswitched docs
    double alpha_init, beta_init, tau_init, epsilon_init;
    double *alpha, *beta0, *beta1, *tau, *epsilon; // LDA hyperparameters 
    double sum_alpha, sum_beta0, sum_beta1, sum_tau,sum_epsilon;
    bool alpha_est;
    int niters; // number of Gibbs sampling iterations
    int burn_iters;  // in inference, after this # of iterations, start burn-in.
    int liter; // the iteration at which the model was saved
    int savestep; // saving period
    int twords; // print out top words per each topic
    int withrawstrs;

    double * p; // temp variable for sampling
    double * p_z_l;
    int ** z; // topic assignments for words, size M x doc.size()
    int ** l; //language assignments for words, size M x doc.size()
    int ** b; //bg assgnments for words, size M x doc.size()
    int ** nw0; // cwt[i][j]: number of instances of word/term i in language_0 assigned to topic j, size V x K
    int ** nw1; // cwt[i][j]: number of instances of word/term i in language_1 assigned to topic j, size V x K
    int *nwb0; //number of instances of word/term i in language_0 assigned to background words, size V
    int *nwb1; //number of instances of word/term i in language_1 assigned to background words, size V
    int ** nd; // na[i][j]: number of words in document i assigned to topic j, size M x K
    int ** nl; //number of words in document i assigned to language l, size M x 2
    int * nb; // number of words assigned to bg size 2
    int * nw0sum; // nwsum[j]: total number of words assigned to topic j, size K
    int * nw1sum;
    int nwb0sum;
    int nwb1sum;
    int * ndsum; // nasum[i]: total number of words in document i, size M
    double ** theta; // theta: document-topic distributions, size M x K
    double ** pi; //pi: document-lanugage distributions, size M * 2
    double ** phi0; // phi: topic-word distributions, size K x V
    double ** phi1;
    double * phi_b0;
    double * phi_b1;
    double * psi;
    double ** gamma_theta;
    double ** gamma_pi;
    double ** gamma_pi_s;
    double ** gamma_phi0;
    double ** gamma_phi1;
    double ** gamma_psi;
    // for inference only
    int inf_liter;
    int newM;
    int newV;
    int ** newz;
    int ** newnw;
    int ** newnd;
    int * newnwsum;
    int * newndsum;
    double ** newtheta;
    double ** newphi;
    // --------------------------------------
    
    model() {
	set_default_values();
    }
          
    ~model();
    
    // set default values for variables
    void set_default_values();   

    // parse command line to get options
    int parse_args(int argc, char ** argv);
    
    // initialize the model
    int init(int argc, char ** argv);
    
    // load LDA model to continue estimating or to do inference
    int load_model(string model_name);

    // load LDA model to do the inference.
    int load_for_inf(string model_name);
    int load_model_alpha(string filename);
    int load_model_tau(string filename);
    int load_model_psi(string filename);
    int load_model_phi(string filename, int lin);
    int load_model_phibg(string filename, int lin);
    char* readline(FILE *input);

    // save LDA model to files
    // model_name.tassign: topic assignments for words in docs
    // model_name.theta: document-topic distributions
    // model_name.phi: topic-word distributions
    // model_name.others: containing other parameters of the model (alpha, beta, M, V, K)
    int save_model(string model_name);
    int save_model_tassign(string filename);
    int save_model_alpha(string filename);
    int save_model_tau(string filename);
    int save_model_theta(string filename);
    int save_model_pi(string filename);
    int save_model_phi(string filename, int lin);
    int save_model_phibg(string filename, int lin);
    int save_model_psi(string filename);
    int save_model_others(string filename);
    int save_model_twords(string filename, int lin);
    
    // saving inference outputs
    int save_inf_model(string model_name);
    int save_inf_model_tassign(string filename);
    int save_inf_model_newtheta(string filename);
    int save_inf_model_newphi(string filename);
    int save_inf_model_others(string filename);
    int save_inf_model_twords(string filename);
    
    // init for estimation
    int init_est();
    int init_estc();
	
    // estimate LDA model using Gibbs sampling
    void estimate();
    void sampling(int m, int n);
    void compute_theta();
    void compute_phi();
    void compute_pi();
    void compute_psi();
    // init for inference
    int init_inf();
    // inference for new (unseen) data based on the estimated LDA model
    void inference();
    int perplexity_per_doc(int m, document * doc);
    int inf_sampling(int m, int n);
    void compute_newtheta();
    void compute_newphi();
    int multinomial(double * p, int K);
    double llhood();
    double hyperpara_est(double * alpha, double ** gamma, int MM, int KK, double step, double converge, int MaxIter);
};

#endif

