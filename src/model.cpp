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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include "constants.h"
#include "strtokenizer.h"
#include "utils.h"
#include "dataset.h"
#include "model.h"

using namespace std;

model::~model() {
    if (p) {
	delete p;
    }

    if (p_z_l) {
        delete p_z_l;
    }

    if (ptrndata) {
	delete ptrndata;
    }
    
    if (pnewdata) {
	delete pnewdata;
    }

    if (z) {
	for (int m = 0; m < M; m++) {
	    if (z[m]) {
		delete z[m];
	    }
	}
	delete z;
    }

    if (l) {
        for (int m = 0; m < M; m++) {
	    if (l[m]) {
		delete l[m];
	    }
	}
	delete l;
    }

    if(b) {
	for (int m = 0; m < M; m++) {
		if (b[m]) {
		    delete b[m];
		}
	}
	delete b;
    }
    
    if (nw0) {
	for (int w = 0; w < V; w++) {
	    if (nw0[w]) {
		delete nw0[w];
	    }
	}
	delete nw0;
    }

    if (nw1) {
	for (int w = 0; w < V; w++) {
	    if (nw1[w]) {
		delete nw1[w];
	    }
	}
	delete nw1;
    }

    if (nwb0) {
	delete nwb0;
    }

    if (nwb1) {
	delete nwb1;
    }

    if (nd) {
	for (int m = 0; m < M; m++) {
	    if (nd[m]) {
		delete nd[m];
	    }
	}
	delete nd;
    }

    if (nl) {
	for (int m = 0; m < M; m++) {
	    if (nl[m]) {
		delete nl[m];
	    }
	}
	delete nl;
    } 

    if (nb) {
	    delete nb;
    }
    
    if (nw0sum) {
	delete nw0sum;
    }   

    if (nw1sum) {
        delete nw1sum;
    }

    if (ndsum) {
	delete ndsum;
    }
    
    if (theta) {
	for (int m = 0; m < M; m++) {
	    if (theta[m]) {
		delete theta[m];
	    }
	}
	delete theta;
    }
   
    if (pi) {
        for (int m = 0; m < M; m++) {
	    if (pi[m]) {
		delete pi[m];
	    }
	}
	delete pi;
    } 

    if (phi0) {
	for (int k = 0; k < K; k++) {
	    if (phi0[k]) {
		delete phi0[k];
	    }
	}
	delete phi0;
    }

    if (phi1) {
	for (int k = 0; k < K; k++) {
	    if (phi1[k]) {
		delete phi1[k];
	    }
	}
	delete phi1;
    }

    if (phi_b0) {
        delete phi_b0;
    }

    if (phi_b1) {
        delete phi_b1;
    }

    if (psi) {
	delete psi;
    }

    if (gamma_theta) {
	for (int m = 0; m < M; m++) {
	    if (gamma_theta[m]) {
		delete gamma_theta[m];
	    }
	}
	delete gamma_theta;
    }
   
    if (gamma_pi) {
        for (int m = 0; m < M; m++) {
	    if (gamma_pi[m]) {
		delete gamma_pi[m];
	    }
	}
	delete gamma_pi;
    } 

     if (gamma_pi_s) {
        for (int m = 0; m < count_cs; m++) {
	    if (gamma_pi_s[m]) {
		delete gamma_pi_s[m];
	    }
	}
	delete gamma_pi_s;
    }
    if (gamma_phi0) {
	for (int k = 0; k < K+1; k++) {
	    if (gamma_phi0[k]) {
		delete gamma_phi0[k];
	    }
	}
	delete gamma_phi0;
    }

    if (gamma_phi1) {
	for (int k = 0; k < K+1; k++) {
	    if (gamma_phi1[k]) {
		delete gamma_phi1[k];
	    }
	}
	delete gamma_phi1;
    }

    if (gamma_psi) {
	    if (gamma_psi[0])
	        delete gamma_psi[0];
	    delete gamma_psi;
    }

    if (alpha) {
        delete alpha;
    }

    if (beta0) {
        delete beta0;
    }

    if (beta1) {
	delete beta1;
    }

    if (tau) {
	delete tau;
    }

    if (epsilon) {
	delete epsilon;
    }

    // only for inference
    if (newz) {
	for (int m = 0; m < newM; m++) {
	    if (newz[m]) {
		delete newz[m];
	    }
	}
    }
    
    if (newnw) {
	for (int w = 0; w < newV; w++) {
	    if (newnw[w]) {
		delete newnw[w];
	    }
	}
    }

    if (newnd) {
	for (int m = 0; m < newM; m++) {
	    if (newnd[m]) {
		delete newnd[m];
	    }
	}
    } 
    
    if (newnwsum) {
	delete newnwsum;
    }   
    
    if (newndsum) {
	delete newndsum;
    }
    
    if (newtheta) {
	for (int m = 0; m < newM; m++) {
	    if (newtheta[m]) {
		delete newtheta[m];
	    }
	}
    }
    
    if (newphi) {
	for (int k = 0; k < K; k++) {
	    if (newphi[k]) {
		delete newphi[k];
	    }
	}
    }
}

void model::set_default_values() {
    wordmapfile = "wordmap.txt";
    trainlogfile = "trainlog.txt";
    tassign_suffix = ".tassign";
    alpha_suffix = ".alpha";
    tau_suffix = ".tau";
    theta_suffix = ".theta";
    pi_suffix = ".pi";
    phi0_suffix = ".phi0";
    phi1_suffix = ".phi1";
    phi0bg_suffix = ".phi0_bg";
    phi1bg_suffix = ".phi1_bg";
    psi_suffix = ".psi";
    others_suffix = ".others";
    twords0_suffix = ".twords0";
    twords1_suffix = ".twords1";
    
    dir = "./";
    dfile = "trndocs.dat";
    model_name = "model-final";    
    model_status = MODEL_STATUS_UNKNOWN;
    
    ptrndata = NULL;
    pnewdata = NULL;
    
    M = 0;
    V = 0;
    K = 100;
    count_cs = 0;
    alpha_est = true;
    alpha_init = 50.0 / K;
    beta_init = 0.1;
    tau_init = 0.5;
    epsilon_init = 0.5;
    niters = 2000;
    liter = 0;
    savestep = 200;    
    twords = 0;    //number of representative words for each topic
    withrawstrs = 0;
    
    alpha = NULL;
    beta0 = NULL;
    beta1 = NULL;
    tau = NULL;
    epsilon = NULL;
    p = NULL;
    p_z_l = NULL;
    z = NULL;
    l = NULL;
    b = NULL;
    nw0 = NULL;
    nw1 = NULL;
    nwb0 = NULL;
    nwb1 = NULL;
    nd = NULL;
    nl = NULL;
    nb = NULL;
    nw0sum = NULL;
    nw1sum = NULL;
    nwb0sum = 0;
    nwb1sum = 0;
    ndsum = NULL;
    theta = NULL;
    pi = NULL;
    phi0 = NULL;
    phi1 = NULL;
    phi_b0 = NULL;
    phi_b1 = NULL;
    psi = NULL;

    gamma_theta = NULL;
    gamma_pi = NULL;
    gamma_pi_s = NULL;
    gamma_phi0 = NULL;
    gamma_phi1 =NULL;
    gamma_psi = NULL;

    newM = 0;
    newV = 0;
    newz = NULL;
    newnw = NULL;
    newnd = NULL;
    newnwsum = NULL;
    newndsum = NULL;
    newtheta = NULL;
    newphi = NULL;
}

int model::parse_args(int argc, char ** argv) {
    return utils::parse_args(argc, argv, this);
}

int model::init(int argc, char ** argv) {
    // call parse_args
    if (parse_args(argc, argv)) {
	return 1;
    }
    
    if (model_status == MODEL_STATUS_EST) {
	// estimating the model from scratch
	if (init_est()) {
	    return 1;
	}
	
    } else if (model_status == MODEL_STATUS_ESTC) {
	// estimating the model from a previously estimated one
	if (init_estc()) {
	    return 1;
	}
	
    } else if (model_status == MODEL_STATUS_INF) {
	// do inference
	if (init_inf()) {
	    return 1;
	}
    }
    
    return 0;
}


int model::load_for_inf(string model_name) {
    if (load_model_alpha(dir + model_name + alpha_suffix)) {
        return 1;
    }
    printf("successfully loaded alpha!\n"); 
    
    if (load_model_tau(dir + model_name + tau_suffix)) {
        return 1;
    }
    printf("successfully loaded tau!\n"); 
    
    if (load_model_psi(dir + model_name + psi_suffix)) {
    	return 1;
    }
    printf("successfully loaded psi!\n"); 
    
    if (load_model_phi(dir + model_name + phi0_suffix, 0)) {
	return 1;
    }
    printf("successfully loaded phi0!\n"); 

    if (load_model_phi(dir + model_name + phi1_suffix, 1)) {
	return 1;
    }
    printf("successfully loaded phi1!\n"); 
    
    if (load_model_phibg(dir + model_name + phi0bg_suffix, 0)) {
	return 1;
    }
    printf("successfully loaded phi_bg0!\n"); 

    if (load_model_phibg(dir + model_name + phi1bg_suffix, 1)) {
	return 1;
    }
    printf("successfully loaded phi_bg1!\n"); 
    
    return 0;

}


int model::load_model_alpha(string filename) {
  //  printf("in load alpha!!\n");
    ifstream fin; 
   /*if (!fin.open(filename.c_str(), ifstream::in)) {
	printf("Cannot open file %d to load model!\n", filename.c_str());
	return 1;
    }*/
    fin.open(filename.c_str(), ifstream::in);
    //char buff[BUFF_SIZE_LONG];
    string line;

    while (!fin.eof()) {
    	getline(fin, line);
	if (line == "")
		continue;
//	printf("in load alpha, line = %s\n", line.c_str());
	strtokenizer strtok(line, " \t\r\n");
	int length = strtok.count_tokens();
	//printf("in loading alpha, the dimension = %d\n", length);
	sum_alpha = 0.0;
	for (int j = 0; j < length; j++) {
		alpha[j] = atof(strtok.token(j).c_str());	
		sum_alpha += alpha[j];
	}
    }

/*    char * pointer = fgets(buff, BUFF_SIZE_LONG-1, fin);
    if (!pointer) {
	    printf("Invalid word-topic assignment file, check the number of docs!\n");
	    return 1;
    }
    line = buff;*/
    fin.close(); 
    return 0;
}


int model::load_model_tau(string filename) {
    /*FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %d to load model!\n", filename.c_str());
	return 1;
    }
    
    char buff[BUFF_SIZE_LONG];
    */
    ifstream fin; 
    fin.open(filename.c_str(), ifstream::in);
    string line;

    while (!fin.eof()) {
	    getline(fin, line);
	    if (line == "")
		    continue;
	    strtokenizer strtok(line, " \t\r\n");
	    int length = strtok.count_tokens();
	    //printf("in loading tau, the dimension = %d\n", length);

	    sum_tau = 0.0;
	    for (int j = 0; j < length; j++) {
		    tau[j] = atof(strtok.token(j).c_str());	
		    sum_tau += tau[j];
	    }
    }
    //fclose(fin);
    fin.close();
    return 0;
}


int model::load_model_psi(string filename) {
/*    FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %d to load model!\n", filename.c_str());
	return 1;
    }
    
    char buff[BUFF_SIZE_LONG];
    string line;

    char * pointer = fgets(buff, BUFF_SIZE_LONG-1, fin);
    if (!pointer) {
	    printf("Invalid word-topic assignment file, check the number of docs!\n");
	    return 1;
    }
    line = buff;*/
    ifstream fin; 
    fin.open(filename.c_str(), ifstream::in);
    string line;

    while (!fin.eof()) {
	    getline(fin, line);
	    if (line == "")
		    continue;
	    strtokenizer strtok(line, " \t\r\n");
	    int length = strtok.count_tokens();
	    //printf("in loading psi, the dimension = %d\n", length);
	    for (int j = 0; j < length; j++) {
		    psi[j] = atof(strtok.token(j).c_str());	
	    }
    }
    //fclose(fin);
    fin.close();
    return 0;
}


int model::load_model_phi(string filename, int lid) {
    int i, j;
    ifstream fin; 
    fin.open(filename.c_str(), ifstream::in);
    string line; 
    i = -1;
    while (!fin.eof()) {
	    getline(fin, line);
	    ++i;
	    if (line == "")
		    continue;
	    //cout << "line = " << line << endl;

	    strtokenizer strtok(line, " \t\r\n");
	    int length = strtok.count_tokens();
	    //printf("in loading phi, the dimension = %d\n", length);	
	    //printf("test whole line! %s\n", strtok.token(length-1).c_str());
	    // V = length;
	    if (lid == 0) {
		    for (j = 0; j < length; j++) {
			    phi0[i][j] = strtod(strtok.token(j).c_str(), NULL);
	//		    printf("%.9e ", phi0[i][j]);
		    }    
	//	    printf("\n");
	    }
	    else {
		    for (j = 0; j < length; j++) {
			    phi1[i][j] = strtod(strtok.token(j).c_str(), NULL);
		    }    
	    }
    }	
    fin.close();
    //fclose(fin);
    
    return 0;
}


int model::load_model_phibg(string filename, int lid) {
    int i, j;
    
/*    FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %d to load model!\n", filename.c_str());
	return 1;
    }

    char buff[BUFF_SIZE_LONG];
    string line;

    char * pointer = fgets(buff, BUFF_SIZE_LONG-1, fin);
    if (!pointer) {
	    printf("Invalid word-topic assignment file, check the number of docs!\n");
	    return 1;
    }

    line = buff;
  */  
    ifstream fin; 
    fin.open(filename.c_str(), ifstream::in);
    string line;

    while (!fin.eof()) {
	    getline(fin, line);
	    if (line == "")
		    continue;
	    strtokenizer strtok(line, " \t\r\n");
	    int length = strtok.count_tokens();

	    //printf("in loading phi_bg, the dimension = %d\n", length);	
	    if (lid == 0) {
		    for (j = 0; j < length; j++) {
			    phi_b0[j] = strtod(strtok.token(j).c_str(), NULL);
		    }    
	    }
	    else {
		    for (j = 0; j < length; j++) {
			    phi_b1[j] = strtod(strtok.token(j).c_str(), NULL);
		    }    
	    }
    } 
    fin.close();
    //fclose(fin);
    
    return 0;
}


/*char* model::readline(FILE *input)
{
	int len;

	if(fgets(line,max_line_len,input) == NULL)
		return NULL;

	while(strrchr(line,'\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *) realloc(line,max_line_len);
		len = (int) strlen(line);
		if(fgets(line+len,max_line_len-len,input) == NULL)
			break;
	}
	return line;
}
*/

int model::load_model(string model_name) {
    int i, j;
    
    string filename = dir + model_name + tassign_suffix;
    FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %d to load model!\n", filename.c_str());
	return 1;
    }
    
    char buff[BUFF_SIZE_LONG];
    string line;

    // allocate memory for z and ptrndata
    z = new int*[M];
    ptrndata = new dataset(M);
    ptrndata->V = V;

    for (i = 0; i < M; i++) {
	char * pointer = fgets(buff, BUFF_SIZE_LONG, fin);
	if (!pointer) {
	    printf("Invalid word-topic assignment file, check the number of docs!\n");
	    return 1;
	}
	
	line = buff;
	strtokenizer strtok(line, " \t\r\n");
	int length = strtok.count_tokens();
	
	vector<int> words;
	vector<int> topics;
	for (j = 0; j < length; j++) {
	    string token = strtok.token(j);
    
	    strtokenizer tok(token, ":");
	    if (tok.count_tokens() != 2) {
		printf("Invalid word-topic assignment line!\n");
		return 1;
	    }
	    
	    words.push_back(atoi(tok.token(0).c_str()));
	    topics.push_back(atoi(tok.token(1).c_str()));
	}
	
	// allocate and add new document to the corpus
	document * pdoc = new document(words);
	ptrndata->add_doc(pdoc, i);
	
	// assign values for z
	z[i] = new int[topics.size()];
	for (j = 0; j < topics.size(); j++) {
	    z[i][j] = topics[j];
	}
    }   
    
    fclose(fin);
    
    return 0;
}

int model::save_model(string model_name) {
    if (save_model_tassign(dir + model_name + tassign_suffix)) {
	return 1;
    }
    
    if (save_model_others(dir + model_name + others_suffix)) {
	return 1;
    }

    if (save_model_alpha(dir + model_name + alpha_suffix)) {
        return 1;
    }
    
    if (save_model_tau(dir + model_name + tau_suffix)) {
        return 1;
    }
    
    if (save_model_theta(dir + model_name + theta_suffix)) {
	return 1;
    }

    if (save_model_pi(dir + model_name + pi_suffix)) {
	return 1;
    }
    
    if (save_model_phi(dir + model_name + phi0_suffix, 0)) {
	return 1;
    }

    if (save_model_phi(dir + model_name + phi1_suffix, 1)) {
	return 1;
    }
    
    if (save_model_phibg(dir + model_name + phi0bg_suffix, 0)) {
	return 1;
    }

    if (save_model_phibg(dir + model_name + phi1bg_suffix, 1)) {
	return 1;
    }

    if (save_model_psi(dir + model_name + psi_suffix)) {
    	return 1;
    }
    
    if (twords > 0) {
	if (save_model_twords(dir + model_name + twords0_suffix, 0)) {
	    return 1;
	}
	if (save_model_twords(dir + model_name + twords1_suffix, 1)) {
	    return 1;
	}
    }
    
    return 0;
}

int model::save_model_tassign(string filename) {
    int i, j;
    
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < ptrndata->M; i++) {    
	for (j = 0; j < ptrndata->docs[i]->length; j++) {
	    fprintf(fout, "%d:%d ", ptrndata->docs[i]->words[j], z[i][j]);
	}
	fprintf(fout, "\n");
    }

    fclose(fout);
    
    return 0;
}

int model::save_model_theta(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < M; i++) {
	for (int j = 0; j < K; j++) {
	    fprintf(fout, "%.9e ", theta[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_model_pi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < M; i++) {
	for (int j = 0; j < 2; j++) {
	    fprintf(fout, "%.9e ", pi[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_model_alpha(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int j = 0; j < K; j++) {
        fprintf(fout, "%.9e ", alpha[j]);
    }
    fprintf(fout, "\n");
    
    fclose(fout);
    
    return 0;
}

int model::save_model_tau(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int j = 0; j < 2; j++) {
        fprintf(fout, "%.9e ", tau[j]);
    }
    fprintf(fout, "\n");
    
    fclose(fout);
    
    return 0;
}

int model::save_model_phi(string filename, int lin) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < K; i++) {
	for (int j = 0; j < V; j++) {
	    if (lin == 0)
		    fprintf(fout, "%.9e ", phi0[i][j]);
	    else if (lin == 1)
		    fprintf(fout, "%.9e ", phi1[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_phibg(string filename, int lin) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < V; i++) {
    	if (lin == 0)
	    fprintf(fout, "%.9e ", phi_b0[i]);
	else if (lin == 1)
	    fprintf(fout, "%.9e ", phi_b1[i]);
    }
    fprintf(fout, "\n");
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_psi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int j = 0; j < 2; j++) {
        fprintf(fout, "%.9e ", psi[j]);
    }
    fprintf(fout, "\n");
    
    fclose(fout);
    
    return 0;
}

int model::save_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    fprintf(fout, "alpha_init=%.9e\n", alpha_init);
    fprintf(fout, "beta_init=%.9e\n", beta_init);
    fprintf(fout, "tau_init=%.9e\n", tau_init);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", M);
    fprintf(fout, "nwords=%d\n", V);
    fprintf(fout, "liter=%d\n", liter);
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_twords(string filename, int lin) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    if (twords > V) {
	twords = V;
    }
    mapid2word::iterator it;
    
    for (int k = 0; k < K; k++) {
	vector<pair<int, double> > words_probs;
	pair<int, double> word_prob;
	for (int w = 0; w < V; w++) {
	    word_prob.first = w;
	    if (lin == 0)
	    	word_prob.second = phi0[k][w];
	    else if (lin == 1)
		word_prob.second = phi1[k][w];
	    words_probs.push_back(word_prob);
	}
    
        // quick sort to sort word-topic probability
	utils::quicksort(words_probs, 0, words_probs.size() - 1);
	
	fprintf(fout, "Topic %dth:\n", k);
	for (int i = 0; i < twords; i++) {
	    it = id2word.find(words_probs[i].first);
	    if (it != id2word.end()) {
		fprintf(fout, "\t%s   %.9e\n", (it->second).c_str(), words_probs[i].second);
	    }
	}
    }
    
    fclose(fout);    
    
    return 0;    
}

int model::save_inf_model(string model_name) {
    if (save_inf_model_tassign(dir + model_name + tassign_suffix)) {
	return 1;
    }
    
    if (save_inf_model_others(dir + model_name + others_suffix)) {
	return 1;
    }
    
    if (save_inf_model_newtheta(dir + model_name + theta_suffix)) {
	return 1;
    }
    
    if (save_inf_model_newphi(dir + model_name + phi0_suffix)) {
	return 1;
    }

    if (twords > 0) {
	if (save_inf_model_twords(dir + model_name + twords0_suffix)) {
	    return 1;
	}
    }
    
    return 0;
}

int model::save_inf_model_tassign(string filename) {
    int i, j;
    
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < pnewdata->M; i++) {    
	for (j = 0; j < pnewdata->docs[i]->length; j++) {
	    fprintf(fout, "%d:%d ", pnewdata->docs[i]->words[j], newz[i][j]);
	}
	fprintf(fout, "\n");
    }

    fclose(fout);
    
    return 0;
}

int model::save_inf_model_newtheta(string filename) {
    int i, j;

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (i = 0; i < newM; i++) {
	for (j = 0; j < K; j++) {
	    fprintf(fout, "%f ", newtheta[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_inf_model_newphi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < K; i++) {
	for (int j = 0; j < newV; j++) {
	    fprintf(fout, "%f ", newphi[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);    
    
    return 0;
}

int model::save_inf_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    fprintf(fout, "alpha=%f\n", alpha_init);
    fprintf(fout, "beta=%f\n", beta_init);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", newM);
    fprintf(fout, "nwords=%d\n", newV);
    fprintf(fout, "liter=%d\n", inf_liter);
    
    fclose(fout);    
    
    return 0;
}

int model::save_inf_model_twords(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    if (twords > newV) {
	twords = newV;
    }
    mapid2word::iterator it;
    map<int, int>::iterator _it;
    
    for (int k = 0; k < K; k++) {
	vector<pair<int, double> > words_probs;
	pair<int, double> word_prob;
	for (int w = 0; w < newV; w++) {
	    word_prob.first = w;
	    word_prob.second = newphi[k][w];
	    words_probs.push_back(word_prob);
	}
    
        // quick sort to sort word-topic probability
	utils::quicksort(words_probs, 0, words_probs.size() - 1);
	
	fprintf(fout, "Topic %dth:\n", k);
	for (int i = 0; i < twords; i++) {
	    _it = pnewdata->_id2id.find(words_probs[i].first);
	    if (_it == pnewdata->_id2id.end()) {
		continue;
	    }
	    it = id2word.find(_it->second);
	    if (it != id2word.end()) {
		fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
	    }
	}
    }
    
    fclose(fout);    
    
    return 0;    
}


int model::init_est() {
    int m, n, w, k;

    p = new double[K];
    p_z_l = new double[K*2];
    for (k = 0; k < K; k ++) {
        p[k] = 0.0;
    }
    for (k = 0; k< 2*K; k ++) {
        p_z_l[k] = 0.0;
    }
    // + read training data
    ptrndata = new dataset;
    if (ptrndata->read_trndata(dir + dfile, dir + wordmapfile)) {
        printf("Fail to read training data!\n");
        return 1;
    }
		
    // + allocate memory and assign values for variables
    M = ptrndata->M;
    V = ptrndata->V;
    printf("k = %d, m = %d, v = %d, lang = %s\n", K, M, V, ptrndata->docs[0]->language.c_str());
    // K: from command line or default value
    // alpha, beta: from command line or default values
    // niters, savestep: from command line or default values

    nw0 = new int*[V];
    nw1 = new int*[V];
    for (w = 0; w < V; w++) {
        nw0[w] = new int[K];
        nw1[w] = new int[K];
	for (k = 0; k < K; k++) {
    	    nw0[w][k] = 0;
            nw1[w][k] = 0;
	}
    }
    nwb0 = new int[V];
    nwb1 = new int[V];
    for (w = 0; w < V; w++) {
	nwb0[w] = 0;
	nwb1[w] = 0;
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }

    nl = new int*[M];
    for (m = 0; m < M; m++) {
        nl[m] = new int[2];
        for (k = 0; k < 2; k++) {
    	    nl[m][k] = 0;
        }
    }
    nb = new int[2];
    nb[0] = nb[1] = 0;

    nw0sum = new int[K];
    nw1sum = new int[K];
    for (k = 0; k < K; k++) {
	nw0sum[k] = 0;
    	nw1sum[k] = 0;
    }
    nwb0sum = 0;
    nwb1sum = 0;
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	ndsum[m] = 0;
    }
    
    srandom(time(0)); // initialize for random number generation
    z = new int*[M];
    l = new int*[M];
    b = new int*[M];
    for (m = 0; m < ptrndata->M; m++) {
	int N = ptrndata->docs[m]->length;
	z[m] = new int[N];
	l[m] = new int[N];
	b[m] = new int[N];
        // initialize for z
        if (ptrndata->docs[m]->language == "codeS") 
		count_cs++;
	for (n = 0; n < N; n++) {
    	    int topic = (int)(((double)random() / RAND_MAX) * K);
	    z[m][n] = topic;
	    int temp_lang;
    	    if (ptrndata->docs[m]->language == "eng") {
	    	l[m][n] = 0; 
		temp_lang = 0;
	    } 
	    else if (ptrndata->docs[m]->language == "spa") {
	        l[m][n] = 1;
		temp_lang = 1;
	    }
	    else {
		//count_cs++;
	        temp_lang = (int)(((double)random() / RAND_MAX) * 2);
	        l[m][n] = temp_lang;
	    }
	    int temp_bg = (int)(((double)random() / RAND_MAX) * 2);
	    b[m][n] = temp_bg;
	    if (b[m][n] == 0) {
    	        // number of instances of word i assigned to topic j
    	        if (temp_lang == 0) {
	            nw0[ptrndata->docs[m]->words[n]][topic] += 1;
	            nw0sum[topic] += 1;
	        }
	        else {
		    nw1[ptrndata->docs[m]->words[n]][topic] += 1;
	            nw1sum[topic] += 1;
	        }
	    }
	    else {
		if (temp_lang == 0) {
		    nwb0[ptrndata->docs[m]->words[n]] += 1;
		    nwb0sum += 1;
		}
		else {
	            nwb1[ptrndata->docs[m]->words[n]] += 1;
		    nwb1sum += 1;
		}
	    }
	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    nl[m][temp_lang] += 1;
	    nb[temp_bg] += 1;
	    // total number of words assigned to topic j
        } 
        // total number of words in document i
        ndsum[m] = N;      
    }
    
    alpha = new double[K];
    for (k = 0; k < K; k ++) {
    	alpha[k] = alpha_init;
    }
    sum_alpha = K * alpha_init;
    beta0 = new double[V];
    beta1 = new double[V];
    for (n = 0; n < V; n++) {
    	beta0[n] = beta_init;
	beta1[n] = beta_init;
    }
    sum_beta0 = V * beta_init;
    sum_beta1 = V * beta_init;
    tau = new double[2];
    for (n=0; n < 2; n++){
        tau[n] = tau_init;
    }
    sum_tau = 2*tau_init;

    epsilon = new double[2];
    for (n=0; n < 2; n++) {
	epsilon[n] = epsilon_init;
    }
    sum_epsilon = 2*epsilon_init;

    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    }
    phi0 = new double*[K];
    phi1 = new double*[K];
    for (k = 0; k < K; k++) {
        phi0[k] = new double[V];
	phi1[k] = new double[V];
    }
    phi_b0 = new double[V];
    phi_b1 = new double[V];

    pi = new double*[M];
    for (m = 0; m < M; m ++) {
        pi[m] = new double[2];
    }
    psi =  new double[2]; 

    gamma_theta = new double*[M];
    for (m = 0; m < M; m++) {
        gamma_theta[m] = new double[K];
    }
    gamma_phi0 = new double*[K+1];
    gamma_phi1 = new double*[K+1];
    for (k = 0; k < K+1; k++) {
        gamma_phi0[k] = new double[V];
	gamma_phi1[k] = new double[V];
    }
    gamma_pi = new double*[M];
    for (m = 0; m < M; m ++) {
        gamma_pi[m] = new double[2];
    }
    gamma_pi_s = new double*[count_cs];
    for (m = 0; m < count_cs; m ++) {
        gamma_pi_s[m] = new double[2];
    }
    gamma_psi = new double*[1];
    gamma_psi[0] = new double[2];

    return 0;
}

int model::init_estc() {
    // estimating the model from a previously estimated one
    int m, n, w, k;

/*    p = new double[K];

    // load model, i.e., read z and ptrndata
    if (load_model(model_name)) {
	printf("Fail to load word-topic assignmetn file of the model!\n");
	return 1;
    }

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }
	
    nwsum = new int[K];
    for (k = 0; k < K; k++) {
	nwsum[k] = 0;
    }
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	ndsum[m] = 0;
    }

    for (m = 0; m < ptrndata->M; m++) {
	int N = ptrndata->docs[m]->length;

	// assign values for nw, nd, nwsum, and ndsum	
        for (n = 0; n < N; n++) {
    	    int w = ptrndata->docs[m]->words[n];
    	    int topic = z[m][n];
    	    
    	    // number of instances of word i assigned to topic j
    	    nw[w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        } 
        // total number of words in document i
        ndsum[m] = N;      
    }
	
    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    }
	
    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }    
*/
    return 0;        
}


void model::estimate() {
    if (twords > 0) {
	// print out top words per topic
	dataset::read_wordmap(dir + wordmapfile, &id2word);
    }

    printf("Sampling %d iterations!\n", niters);
    double log_lhood;
    int last_iter = liter;
    for (liter = last_iter + 1; liter <= niters + last_iter; liter++) {
	printf("Iteration %d ...\n", liter);
	// for all z_i
	for (int m = 0; m < M; m++) {
	    for (int n = 0; n < ptrndata->docs[m]->length; n++) {
		// (z_i = z[m][n])
		// sample from p(z_i|z_-i, w)
		    sampling(m, n);
	    }
	}
	if (savestep > 0) {
	    if (liter % savestep == 0) {
		// saving the model
		printf("Saving the model at iteration %d ...\n", liter);
		compute_theta();
		compute_phi();
		compute_pi();
		compute_psi();
		int idx = 0;
		for(int m = 0; m < M; m++) {
			if (ptrndata->docs[m]->language=="codeS") {
				for (int ll = 0; ll < 2; ll++) {
					gamma_pi_s[idx][ll] = gamma_pi[m][ll];
				}
				idx++;
			}
		}
		if (alpha_est)
			sum_alpha=hyperpara_est(alpha,gamma_theta,M,K,0.5,0.01,20);
		sum_tau=hyperpara_est(tau,gamma_pi_s,count_cs,2,0.5,0.0001,20);
//		fprintf(stderr, "tau = %f, %f\n", tau[0], tau[1]);
		sum_epsilon=hyperpara_est(epsilon,gamma_psi,1,2,0.5,0.01,20);
//		fprintf(stderr, "epsilon = %f, %f\n", epsilon[0], epsilon[1]);
		sum_beta0=hyperpara_est(beta0,gamma_phi0,K+1,V,0.5,0.01,20);
		sum_beta1=hyperpara_est(beta1,gamma_phi1,K+1,V,0.5,0.01,20);
/*		fprintf(stderr, "beta0 = ");
		for (int m = 0; m < V; m ++){
			fprintf(stderr, "%f ", beta0[m]);
		}
		fprintf(stderr, "\n");
		fprintf(stderr, "beta1 = ");
		for (int m = 0; m < V; m ++){
			fprintf(stderr, "%f ", beta1[m]);
		}
		fprintf(stderr, "\n");*/
		log_lhood = llhood();
		fprintf(stderr, "iter=%d, loglhood=%f\n", liter, log_lhood);
		save_model(utils::generate_model_name(liter));
	    }
	}
    }
    
    printf("Gibbs sampling completed!\n");
    printf("Saving the final model!\n");
    if (alpha_est)
	    sum_alpha=hyperpara_est(alpha,gamma_theta,M,K,0.5,0.01,20);
    sum_tau=hyperpara_est(tau,gamma_pi_s,count_cs,2,0.5,0.001,20);
    sum_epsilon=hyperpara_est(epsilon,gamma_psi,1,2,0.5,0.01,20);
    sum_beta0=hyperpara_est(beta0,gamma_phi0,K+1,V,0.5,0.01,20);
    sum_beta1=hyperpara_est(beta1,gamma_phi1,K+1,V,0.5,0.01,20);
    
    compute_theta();
    compute_phi();
    compute_pi();
    compute_psi();
    liter--;
    log_lhood = llhood();
    fprintf(stderr, "iter=%d, loglhood=%.9e\n", liter, log_lhood);
    save_model(utils::generate_model_name(-1));
}

void model::sampling(int m, int n) {
    // remove z_i from the count variables
    int index;
    int topic = z[m][n];
    int w = ptrndata->docs[m]->words[n];
    string g = ptrndata->docs[m]->language;
    int ll = l[m][n];
    int bg = b[m][n];
    nd[m][topic] -= 1;
    if (g == "codeS")
	nl[m][ll] -= 1;
    if (bg == 0) {
        if (ll == 0) {
	    nw0[w][topic] -= 1;
	    nw0sum[topic] -= 1;
        }
        else {
    	    nw1[w][topic] -= 1;
	    nw1sum[topic] -= 1;
	}
    }
    else {
        if (ll == 0) {
            nwb0[w] -= 1;
            nwb0sum -= 1;
	}
	else {
	    nwb1[w] -= 1;
	    nwb1sum -= 1;
	}	
    }
    // do multinomial sampling via cumulative method
    // sample z,l
    if (g == "codeS") {
	if (bg == 0) {
    	    for (int k = 0; k < K; k++) {
		p_z_l[k] = (nd[m][k] + alpha[k]) * (nl[m][0]+tau[0]) * (nw0[w][k] + beta0[w]) / (nw0sum[k] + sum_beta0);
    	    }
	    for (int k = 0; k < K; k++) {
		p_z_l[k+K] = (nd[m][k] + alpha[k]) * (nl[m][1]+tau[1]) * (nw1[w][k] + beta1[w]) / (nw1sum[k] + sum_beta1);
	    }
	}
        else {
	    for (int k = 0; k < K; k++) {
		p_z_l[k] = (nd[m][k] + alpha[k]) * (nl[m][0]+tau[0]) * (nwb0[w] + beta0[w]) / (nwb0sum + sum_beta0);
    	    }
	    for (int k = 0; k < K; k++) {
		p_z_l[k+K] = (nd[m][k] + alpha[k]) * (nl[m][1]+tau[1]) * (nwb1[w] + beta1[w]) / (nwb1sum + sum_beta1);
	    }    
	}
    	//sample a multinomial distribution
	index = multinomial(p_z_l,K*2);
    	z[m][n] = index % K;
    	l[m][n] = index / K;
    }
    else { 
	if (bg == 0) {
	    if (ll == 0) {
		for (int k = 0; k < K; k++) {
	    	    p[k] = (nd[m][k] + alpha[k]) * (nw0[w][k] + beta0[w]) / (nw0sum[k] + sum_beta0);
	    	}
	    }
	    else {
		for (int k = 0; k < K; k++) {
		    p[k] = (nd[m][k] + alpha[k]) * (nw1[w][k] + beta1[w]) / (nw1sum[k] + sum_beta1);
		}
	    }
	}
	else {
	    for (int k = 0; k < K; k++) {
	        p[k] = (nd[m][k] + alpha[k]);
	    }
	}
	index = multinomial(p,K);
	z[m][n] = index;
    }
    // add and exclude neccessary counts
    topic = z[m][n];
    ll = l[m][n];
    nd[m][topic] += 1;
    if (g == "codeS") {
	    nl[m][ll] += 1;
    }
    nb[bg] -= 1;

    // sample b
    if (ll == 0) {
        p[0] = (nb[0] + epsilon[0])*(nw0[w][topic]+beta0[w])/(nw0sum[topic]+sum_beta0);
        p[1] = (nb[1] + epsilon[1])*(nwb0[w]+beta0[w])/(nwb0sum+sum_beta0);
    }
    else {
        p[0] = (nb[0] + epsilon[0])*(nw1[w][topic]+beta1[w])/(nw1sum[topic]+sum_beta1);
        p[1] = (nb[1] + epsilon[1])*(nwb1[w]+beta1[w])/(nwb1sum+sum_beta1);
    }
    index = multinomial(p,2);
    b[m][n] = index;
    bg = b[m][n];
    // update counts
    nb[bg] += 1;
    if (bg == 0) {
        if (ll == 0) {
    	    nw0[w][topic] += 1;
    	    nw0sum[topic] += 1;
        }
        else {
	    nw1[w][topic] += 1;
	    nw1sum[topic] += 1;
        }
    }
    else {
	if (ll == 0) {
	    nwb0[w] += 1;
	    nwb0sum += 1;
	}
	else {
	    nwb1[w] += 1;
	    nwb1sum += 1;
	}
    }
}

int model::multinomial(double *p, int K) {
    int t;
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
	    p[k] += p[k-1];
    }
    // scaled sample because of unnormalized p[]
    double u = ((double)random() / RAND_MAX) * p[K-1];

    for (t = 0; t < K; t++) {
	    if (p[t] > u) {
		return t;
	    }
    }
    return K-1;
}
void model::compute_theta() {
    for (int m = 0; m < M; m++) {
	int ndsum = ptrndata->docs[m]->length;
	for (int k = 0; k < K; k++) {
	    gamma_theta[m][k] = nd[m][k]+alpha[k];
	    theta[m][k] =gamma_theta[m][k] / (ndsum + sum_alpha);
	}
    }
}

void model::compute_phi() {
    for (int k = 0; k < K; k++) {
	for (int w = 0; w < V; w++) {
	    gamma_phi0[k][w] = nw0[w][k]+beta0[w];
	    phi0[k][w] = gamma_phi0[k][w] / (nw0sum[k] + sum_beta0);
	    gamma_phi1[k][w] = nw1[w][k]+beta1[w];
	    phi1[k][w] = gamma_phi1[k][w] / (nw1sum[k] + sum_beta1);
	}
    }
    
    for (int w = 0; w < V; w++) {
        gamma_phi0[K][w] = nwb0[w]+beta0[w];
	phi_b0[w] = gamma_phi0[K][w] / (nwb0sum + sum_beta0);
	gamma_phi1[K][w] = nwb1[w]+beta1[w];
	phi_b1[w] = gamma_phi1[K][w] / (nwb1sum + sum_beta1);
    }
}

void model::compute_pi() {
    for (int m = 0; m < M; m++) {
	string g = ptrndata->docs[m]->language;
	if (g=="codeS") {
		int nlsum = ptrndata->docs[m]->length;
        	for (int l = 0; l < 2; l++) {
			gamma_pi[m][l] = nl[m][l]+tau[l];
			pi[m][l] = gamma_pi[m][l] / (nlsum + sum_tau);
		}
	}
    }
}

void model::compute_psi() {
    int nbsum = nb[0]+nb[1];
    gamma_psi[0][0] = nb[0]+epsilon[0];
    psi[0] = gamma_psi[0][0] / (nbsum + sum_epsilon);
    gamma_psi[0][1] = nb[1]+epsilon[1];
    psi[1] = gamma_psi[0][1] / (nbsum + sum_epsilon);
}

double model::llhood() {
    double llhood = 0.0;
    int w;
    double inner_sum,inner_inner_sum;
    for (int m = 0; m < M; m++) {
	string g = ptrndata->docs[m]->language;
	if (g == "codeS") {
        	for (int n = 0; n < ptrndata->docs[m]->length; n++) {
			w = ptrndata->docs[m]->words[n];
			inner_sum = 0.0;
			for (int k = 0; k < K; k++) {
				inner_inner_sum = 0.0;
				inner_inner_sum += pi[m][0]*phi0[k][w];
				inner_inner_sum += pi[m][1]*phi1[k][w];
				inner_sum += theta[m][k]*inner_inner_sum;
			}
			inner_sum *= psi[0];
			inner_inner_sum = 0.0;
			inner_inner_sum += pi[m][0]*phi_b0[w];
			inner_inner_sum += pi[m][1]*phi_b1[w];
			inner_sum += psi[1]*inner_inner_sum;
			llhood += log(inner_sum);
		}
	}
	else if(g == "eng") {
        	for (int n = 0; n < ptrndata->docs[m]->length; n++) {
			w = ptrndata->docs[m]->words[n];
			inner_sum = 0.0;
			for (int k = 0; k < K; k++) {
				inner_inner_sum = phi0[k][w];
				inner_sum += theta[m][k]*inner_inner_sum;

			}
			inner_sum *= psi[0];
			inner_inner_sum = phi_b0[w];
			inner_sum += psi[1]*inner_inner_sum;
			llhood += log(inner_sum);
		}
	}
	else if(g == "spa") {
        	for (int n = 0; n < ptrndata->docs[m]->length; n++) {
			w = ptrndata->docs[m]->words[n];
			inner_sum = 0.0;
			for (int k = 0; k < K; k++) {
				inner_inner_sum = phi1[k][w];
				inner_sum += theta[m][k]*inner_inner_sum;
			}
			inner_sum *= psi[0];
			inner_inner_sum = phi_b1[w];
			inner_sum += psi[1]*inner_inner_sum;
			llhood += log(inner_sum);
		}
	}
    }
    return llhood;
}

int model::init_inf() {
    int m, n, w, k;

    // + read testing data
    pnewdata = new dataset;
    if (pnewdata->read_newdata(dir + dfile, dir + wordmapfile)) {
        printf("Fail to read test data!\n");
        return 1;
    }
		
    // + allocate memory and assign values for variables
    M = pnewdata->M;
    newM = pnewdata->M;
    V = pnewdata->V;
    printf("k = %d, m = %d, v = %d\n", K, M, V);
    // K: from command line or default value
    // niters, savestep: from command line or default values

    alpha = new double[K];
    tau = new double[2];
    psi =  new double[2]; 
    phi0 = new double*[K];
    phi1 = new double*[K];
    for (k = 0; k < K; k++) {
    	phi0[k] = new double[V];
    	phi1[k] = new double[V];
    }
    phi_b0 = new double[V];
    phi_b1 = new double[V];
    if (load_for_inf(utils::generate_model_name(-1))) {
    	printf("fail to load the model!\n");
	return 1;
    }

    /*printf("printing out phi0 and phi1:\n");
    for (k = 0; k < K; k ++) {
    	for (w = 0; w < V; w++) {
		printf("%.9e, %.9e\t", phi0[k][w], phi1[k][w]);
	}
	printf("\n");
    }
    printf("printing out phi_bg0 and phi_bg1:\n");
    for (w = 0; w < V; w++) {
	    printf("%f, %f\t", phi_b0[w], phi_b1[w]);
    }
    printf("\n");
*/

    /*    beta0 = new double[V];
    beta1 = new double[V];
    for (n = 0; n < V; n++) {
    	beta0[n] = beta_init;
	beta1[n] = beta_init;
    }
    sum_beta0 = V * beta_init;
    sum_beta1 = V * beta_init;

    epsilon = new double[2];
    for (n=0; n < 2; n++) {
	epsilon[n] = epsilon_init;
    }
    sum_epsilon = 2*epsilon_init;
*/
    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    	for (k = 0; k < K; k ++)
		theta[m][k] = 0.0;
    }

    pi = new double*[M];
    for (m = 0; m < M; m ++) {
        pi[m] = new double[2];
	for (n = 0; n < 2; n++)
		pi[m][n] = 0.0;
    }

/*    gamma_theta = new double*[M];
    for (m = 0; m < M; m++) {
        gamma_theta[m] = new double[K];
    }
 
    gamma_psi = new double*[1];
    gamma_psi[0] = new double[2];
    gamma_phi0 = new double*[K+1];
    gamma_phi1 = new double*[K+1];
    for (k = 0; k < K+1; k++) {
        gamma_phi0[k] = new double[V];
	gamma_phi1[k] = new double[V];
    }
    gamma_pi = new double*[M];
    for (m = 0; m < M; m ++) {
        gamma_pi[m] = new double[2];
    }
    gamma_pi_s = new double*[count_cs];
    for (m = 0; m < count_cs; m ++) {
        gamma_pi_s[m] = new double[2];
    }*/
    p = new double[K];
    p_z_l = new double[K*2];
    for (k = 0; k < K; k ++) {
        p[k] = 0.0;
    }
    for (k = 0; k< 2*K; k ++) {
        p_z_l[k] = 0.0;
    }
    
/*    nw0 = new int*[V];
    nw1 = new int*[V];
    for (w = 0; w < V; w++) {
        nw0[w] = new int[K];
        nw1[w] = new int[K];
	for (k = 0; k < K; k++) {
    	    nw0[w][k] = 0;
            nw1[w][k] = 0;
	}
    }
    nwb0 = new int[V];
    nwb1 = new int[V];
    for (w = 0; w < V; w++) {
	nwb0[w] = 0;
	nwb1[w] = 0;
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }

    nl = new int*[M];
    for (m = 0; m < M; m++) {
        nl[m] = new int[2];
        for (k = 0; k < 2; k++) {
    	    nl[m][k] = 0;
        }
    }
    nb = new int[2];
    nb[0] = nb[1] = 0;

    nw0sum = new int[K];
    nw1sum = new int[K];
    for (k = 0; k < K; k++) {
	nw0sum[k] = 0;
    	nw1sum[k] = 0;
    }
    nwb0sum = 0;
    nwb1sum = 0;
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	ndsum[m] = 0;
    }
    
    srandom(time(0)); // initialize for random number generation
    z = new int*[M];
    l = new int*[M];
    b = new int*[M];
    for (m = 0; m < pnewdata->M; m++) {
	int N = pnewdata->docs[m]->length;
	z[m] = new int[N];
	l[m] = new int[N];
	b[m] = new int[N];
        // initialize for z
        if (pnewdata->docs[m]->language == "codeS") 
		count_cs++;
	for (n = 0; n < N; n++) {
    	    int topic = (int)(((double)random() / RAND_MAX) * K);
	    z[m][n] = topic;
	    int temp_lang;
    	    if (pnewdata->docs[m]->language == "eng") {
	    	l[m][n] = 0; 
		temp_lang = 0;
	    } 
	    else if (pnewdata->docs[m]->language == "spa") {
	        l[m][n] = 1;
		temp_lang = 1;
	    }
	    else {
		//count_cs++;
	        temp_lang = (int)(((double)random() / RAND_MAX) * 2);
	        l[m][n] = temp_lang;
	    }
	    int temp_bg = (int)(((double)random() / RAND_MAX) * 2);
	    b[m][n] = temp_bg;
	    if (b[m][n] == 0) {
    	        // number of instances of word i assigned to topic j
    	        if (temp_lang == 0) {
	            nw0[pnewdata->docs[m]->words[n]][topic] += 1;
	            nw0sum[topic] += 1;
	        }
	        else {
		    nw1[pnewdata->docs[m]->words[n]][topic] += 1;
	            nw1sum[topic] += 1;
	        }
	    }
	    else {
		if (temp_lang == 0) {
		    nwb0[pnewdata->docs[m]->words[n]] += 1;
		    nwb0sum += 1;
		}
		else {
	            nwb1[pnewdata->docs[m]->words[n]] += 1;
		    nwb1sum += 1;
		}
	    }
	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    nl[m][temp_lang] += 1;
	    nb[temp_bg] += 1;
	    // total number of words assigned to topic j
        } 
        // total number of words in document i
        ndsum[m] = N;      
    }
  */  
    return 0;
}

void model::inference() {
    if (twords > 0) {
	// print out top words per topic
	dataset::read_wordmap(dir + wordmapfile, &id2word);
    }

    printf("Sampling %d iterations for inference!\n", niters);
   
    double aggr_perplexity = 0.0; 
    // for all newz_i
    int validCount = 0;
    for (int m = 0; m < newM; m++) {
	    // (newz_i = newz[m][n])
	    // sample from p(z_i|z_-i, w)
	    printf("Processing the %dth document!\n", m);
	    if (pnewdata->docs[m]->length <= 1)
		    continue;
	    validCount += 1;
	    aggr_perplexity += perplexity_per_doc(m, pnewdata->docs[m]);
    	    //printf("perplexity = %f\n", aggr_perplexity);
    }
    
    printf("Gibbs sampling for inference completed!\n");
    printf("Saving the inference outputs!\n");
    inf_liter--;
    fprintf(stderr, "%f\n", aggr_perplexity/validCount);
}

int model::perplexity_per_doc(int m, document * doc) {
	int n = doc->length;
	int half = n/2;
	int index;
	int i, j, k, w;
	int sumK, sumL;
	int * dz = new int[half]; 
	int * dl = new int[half]; 
	int * db = new int[half]; 
	int * dnk = new int[K]; 
	int * dnl = new int[2];
	for (i = 0; i < K; i ++)
		dnk[i] = 0;
	for (i = 0; i < 2; i ++)
		dnl[i] = 0;

	for (i = 0; i < half; i++) {
    	    dz[i] = (int)(((double)random() / RAND_MAX) * K);
	    dnk[dz[i]] += 1;
    	    dl[i] = (int)(((double)random() / RAND_MAX) * 2);
    	    dnl[dl[i]] += 1;
	    db[i] = (int)(((double)random() / RAND_MAX) * 2);
	}
	//printf("doc length = %d, half = %d, K = %d, M = %d, m = %d, burn-in iterations = %d\n", n, half, K, M, m, burn_iters);

    	for (inf_liter = 1; inf_liter <= niters; inf_liter++) {
	//    printf("Iteration %d ...\n", inf_liter);
	    for (i = 0; i < half; i ++) {
	        dnk[dz[i]] -= 1;
		dnl[dl[i]] -= 1;
		//sample z,l
    		w = doc->words[i];
		//printf("w_id = %d,\t", w);
		if (db[i] == 0) {
		    for (k = 0; k < K; k ++)
			p_z_l[k] = (dnk[k] + alpha[k])*(dnl[0] + tau[0])*phi0[k][w];
		    for (k = 0; k < K; k ++)
			p_z_l[k+K] = (dnk[k] + alpha[k])*(dnl[1] + tau[1])*phi1[k][w];
		}	
		else {
		    for (k = 0; k < K; k ++)
			p_z_l[k] = (dnk[k] + alpha[k])*(dnl[0] + tau[0])*phi_b0[w];
		    for (k = 0; k < K; k ++)
   			p_z_l[k+K] = (dnk[k] + alpha[k])*(dnl[1] + tau[1])*phi_b1[w];
		}
		//printf("complete assignment!!\n");
		//sample a multinomial distribution
		index = multinomial(p_z_l,K*2);
		dz[i] = index % K;
		dl[i] = index / K;
		// update the counts to include the newly sampled assignments of the current token
		dnk[dz[i]] += 1;
		dnl[dl[i]] += 1;
		//printf("complete sampling l and z!!\n");
		// sample b
		if (dl[i] == 0) {
			p[0] = psi[0]*phi0[dz[i]][w];
			p[1] = psi[1]*phi_b0[w];
		}
		else {
			p[0] = psi[0]*phi1[dz[i]][w];
			p[1] = psi[1]*phi_b1[w];
		}
		db[i] = multinomial(p, 2); 
		//printf("complete sampling b!!\n");
	    }
	    if (inf_liter > burn_iters) {
	   	sumK = half + sum_alpha;
		for (k = 0; k < K; k ++)
		    theta[m][k] += (dnk[k] + alpha[k])/sumK;
	        sumL = half + sum_tau;	
		for (j = 0; j < 2; j ++)
		    pi[m][j] += (dnl[j] + tau[j])/sumL;
	    }
	}
	for (k = 0; k < K; k ++) {
		theta[m][k] /= niters - burn_iters;
	}
	for (j = 0; j < 2; j ++){
		pi[m][j] /= niters - burn_iters;
	}
	//given theta, pi and phi, compute llhood for the rest half
	//according to sum_i log(sum_k(phi[k][doc[i]] * theta_burnIn[k]))
	double llhood = 0.0;
	for (i = half; i < n; i ++) {
    		w = doc->words[i];
		//printf("w_id = %d,\t", w);
		double inner_sum = 0.0;
		for (k = 0; k < K; k ++) {
		        double inner_inner_sum = 0.0;
			inner_inner_sum += pi[m][0]*phi0[k][w];
			inner_inner_sum += pi[m][1]*phi1[k][w];
			inner_sum += inner_inner_sum * theta[m][k];
		}
		inner_sum *= psi[0];
		double inner_inner_sum = 0.0;
		inner_inner_sum += pi[m][0]*phi_b0[w];
		inner_inner_sum += pi[m][1]*phi_b1[w];
		inner_sum += psi[1]*inner_inner_sum;
		llhood += log(inner_sum);
		//printf("inner_sum = %f, temp log-likelihood = %f\n", inner_sum, llhood);
	}
	double perplexity = exp(-llhood / (n - half));
	
	if (dz) {
		delete dz;
	}	
	if (dl) {
		delete dl;
	}	
	if (db) {
		delete db;
	}	
	if (dnk) {
		delete dnk;
	}	
	if (dnl) {
		delete dnl;
	}	
	
	return perplexity;

}

void model::compute_newtheta() {
/*    for (int m = 0; m < newM; m++) {
	for (int k = 0; k < K; k++) {
	    newtheta[m][k] = (newnd[m][k] + alpha) / (newndsum[m] + K * alpha);
	}
    }*/
}

void model::compute_newphi() {
/*    map<int, int>::iterator it;
    for (int k = 0; k < K; k++) {
	for (int w = 0; w < newV; w++) {
	    it = pnewdata->_id2id.find(w);
	    if (it != pnewdata->_id2id.end()) {
		newphi[k][w] = (nw[it->second][k] + newnw[w][k] + beta) / (nwsum[k] + newnwsum[k] + V * beta);
	    }
	}
    }*/
}

/**hyper parameter estimation using Newton-Raphson methods
 * Input: 
 * alpha0[i], (i = 1~K. K=#of topics, or languages)
 * gamma[d][i] = #of topic (or language) i occurs in document d (average over several gibbs samples), d=1~M
 * Output:
 * alpha_sum: suitable for alpha, beta, etc.
 */
double model::hyperpara_est(double * alpha, double ** gamma, int MM, int KK, double step, double converge, int MaxIter){
	int i, d, iter;
	double alpha_sum, llhood, max_gradient, z, c_up, c_down, c;
	double * gamma_sum, * psi_gamma, * gradient, * h, * s;
        utils myutil;
	gamma_sum = new double[MM];
	for (i = 0; i < MM; i++) {
		gamma_sum[i] = 0.0;	
	}
	gradient = new double[KK];
	h = new double[KK];
	s = new double[KK];
	psi_gamma = new double[KK];
	for (i = 0; i < KK; i ++) {
		gradient[i] = 0.0;
		h[i] = 0.0;
		s[i] = 0.0;
		psi_gamma[i] = 0.0;
	}
	for (d = 0; d < MM; d++) {
		for (i = 0; i < KK; i ++) {
			gamma_sum[d] += gamma[d][i];
		}
	}
	for (i = 0; i < KK; i ++ ) {
		for (d = 0; d < MM; d ++) {
			psi_gamma[i] += myutil.digamma(gamma[d][i]) - myutil.digamma(gamma_sum[d]);
		}
	}
	iter = 0;
	while (iter < MaxIter){
		iter += 1;
		alpha_sum = 0.0;
		for (i = 0; i < KK; i ++){
			alpha_sum += alpha[i];
		}
                //compute log-likelihood llhood = sum_d(log(gammaFunc(alpha_sum)) - sum_i(log(gammaFunc(alpha[i])) + sum_i((alpha[i]-1)*(digamma(gamma[d][i]) - digamma(sum_gamma[d]))))
		llhood = MM * myutil.log_gamma(alpha_sum);
		for (i = 0; i < KK; i ++) {
			llhood -= MM * myutil.log_gamma(alpha[i]);
			llhood += (alpha[i]-1) * psi_gamma[i];
		}
//compute the gradient compute gradient g (KK by 1), g[i] = MM*(digamma(sum_a)-digamma(alpha[i])) + sum_d{ digamma(gamma[d][i])-digamma(sum_ga[d]) }
		max_gradient = 0;
		for (i = 0; i < KK; i++) {
			gradient[i] = MM*(myutil.digamma(alpha_sum)-myutil.digamma(alpha[i])) + psi_gamma[i];
			if (fabs(gradient[i]) > max_gradient)
				max_gradient = fabs(gradient[i]);
		}
//		print 'iter=', iter, 'likelihood=', llhood, 'max-gradient=', max_gradient
		if (max_gradient <= converge)
			return alpha_sum;
//compute Newton direction(s=H^-1*g) using Newton-Raphson methods
//s[i] = (g[i]-c)/h[i]
//h[i] = -M*trigamma(alpha[i])
//c = sum_i{g[i]/h[i]} / (1/z + sum_i{1/h[i]})
//z = M*trigamma(sum_i{alpha[i]})
		z = MM * myutil.trigamma(alpha_sum);
		for (i = 0; i < KK; i ++) {
			h[i]=-MM * myutil.trigamma(alpha[i]);
		}
		c_up = 0;
		c_down = 1/z;
		for (i = 0; i < KK; i ++) {
			c_up += gradient[i]/h[i];
			c_down += 1/h[i];
		}
		c = c_up/c_down;
		for (i = 0; i < KK; i ++) {
			s[i]=(gradient[i]-c)/h[i];//update alpha according to Newton's method
		}
		for (i = 0; i < KK; i ++) {
		    alpha[i] -= step * s[i];	
	            if (alpha[i] < 0) {
		        alpha[i] = 1e-2;
		    }	
		}
	}
	alpha_sum = 0.0;
	for (i = 0; i < KK; i ++){
		alpha_sum += alpha[i];
	}
	return alpha_sum;
}
