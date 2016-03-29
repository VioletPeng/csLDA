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

#include <stdio.h>
#include <cstdlib>
#include <string>
#include <map>
#include "strtokenizer.h"
#include "utils.h"
#include "model.h"

using namespace std;

int utils::parse_args(int argc, char ** argv, model * pmodel) {
    int model_status = MODEL_STATUS_UNKNOWN;
    string dir = "";
    string model_name = "";
    string dfile = "";
    double alpha = -1.0;
    double beta = -1.0;
    double tau = -1.0;
    double epsilon = -1.0;
    int K = 0;
    int niters = 0;
    int burn_iters = 0;
    int savestep = 0;
    int twords = 0;
    int withrawdata = 0;

    int i = 0; 
    while (i < argc) {
	string arg = argv[i];
	
	if (arg == "-est") {
	    model_status = MODEL_STATUS_EST;
	    
	} else if (arg == "-estc") {
	    model_status = MODEL_STATUS_ESTC;
	    
	} else if (arg == "-inf") {
	    model_status = MODEL_STATUS_INF;
	    
	} else if (arg == "-dir") {
	    dir = argv[++i];	    
	    
	} else if (arg == "-dfile") {
	    dfile = argv[++i];	    
	    
	} else if (arg == "-model") {
	    model_name = argv[++i];	    	    
	    
	} else if (arg == "-alpha") {
	    alpha = atof(argv[++i]);	    
	    
	} else if (arg == "-beta") {
	    beta = atof(argv[++i]);	    
	    
	} else if (arg == "-tau") {
	    tau = atof(argv[++i]);

	} else if (arg == "-epsilon") {
	    epsilon = atof(argv[++i]);

	} else if (arg == "-ntopics") {
	    K = atoi(argv[++i]);	    
	    
	} else if (arg == "-niters") {
	    niters = atoi(argv[++i]);	    
	    
	} else if (arg == "-burniters") {
	    burn_iters = atoi(argv[++i]);	    
	    
	} else if (arg == "-savestep") {
	    savestep = atoi(argv[++i]);
	    
	} else if (arg == "-twords") {
	    twords = atoi(argv[++i]);
	    
	} else if (arg == "-withrawdata") {
	    withrawdata = 1;
	
	} else {
	    // any more?
	}	
		
	i++;
    }
    
    if (model_status == MODEL_STATUS_EST) {
	if (dfile == "") {
	    printf("Please specify the input data file for model estimation!\n");
	    return 1;
	}
	
	pmodel->model_status = model_status;
	
	if (K > 0) {
	    pmodel->K = K;
	}
	
	if (alpha >= 0.0) {
	    pmodel->alpha_init = alpha;
	} else {
	    // default value for alpha
	    pmodel->alpha_init = 50.0 / pmodel->K;
	}
	
	if (beta >= 0.0) {
	    pmodel->beta_init = beta;
	}
	
	if (tau >= 0.0) {
	    pmodel->tau_init = tau;
	}

	if (epsilon >= 0.0) {
	    pmodel->epsilon_init = epsilon;
	}

	if (niters > 0) {
	    pmodel->niters = niters;
	}
	
	if (savestep > 0) {
	    pmodel->savestep = savestep;
	}
	
	if (twords > 0) {
	    pmodel->twords = twords;
	}
	
	pmodel->dfile = dfile;
	
	string::size_type idx = dfile.find_last_of("/");			
	if (idx == string::npos) {
	    pmodel->dir = "./";
	} else {
	    pmodel->dir = dfile.substr(0, idx + 1);
	    pmodel->dfile = dfile.substr(idx + 1, dfile.size() - pmodel->dir.size());
	    printf("dir = %s\n", pmodel->dir.c_str());
	    printf("dfile = %s\n", pmodel->dfile.c_str());
	}
    } 
    
    if (model_status == MODEL_STATUS_ESTC) {
	if (dir == "") {
	    printf("Please specify model directory!\n");
	    return 1;
	}
	
	if (model_name == "") {
	    printf("Please specify model name upon that you want to continue estimating!\n");
	    return 1;
	}	

	pmodel->model_status = model_status;

	if (dir[dir.size() - 1] != '/') {
	    dir += "/";
	}
	pmodel->dir = dir;

	pmodel->model_name = model_name;

	if (niters > 0) {
	    pmodel->niters = niters;
	}
	
	if (savestep > 0) {
	    pmodel->savestep = savestep;
	}
	
	if (twords > 0) {
	    pmodel->twords = twords;
	}
	
	// read <model>.others file to assign values for ntopics, alpha, beta, etc.
	if (read_and_parse(pmodel->dir + pmodel->model_name + pmodel->others_suffix, pmodel)) {
	    return 1;
	}	
    } 
    
    if (model_status == MODEL_STATUS_INF) {
	if (dir == "") {
	    printf("Please specify model directory please!\n");
	    return 1;
	}
	
	if (model_name == "") {
	    printf("Please specify model name for inference!\n");
	    return 1;
	}	

	if (dfile == "") {
	    printf("Please specify the new data file for inference!\n");
	    return 1;
	}
	
	pmodel->model_status = model_status;

	if (dir[dir.size() - 1] != '/') {
	    dir += "/";
	}
	pmodel->dir = dir;
	
	pmodel->model_name = model_name;

	pmodel->dfile = dfile;

	if (niters > 0) {
	    pmodel->niters = niters;
	} else {
	    // default number of Gibbs sampling iterations for doing inference
	    pmodel->niters = 20;
	}
	
	if (burn_iters > 0) {
	    pmodel->burn_iters = burn_iters;
	} else {
	    // default number of Gibbs sampling iterations for doing inference
	    pmodel->niters = 100;
	}
	
	if (twords > 0) {
	    pmodel->twords = twords;
	}
	
	if (withrawdata > 0) {
	    pmodel->withrawstrs = withrawdata;
	}
		
	// read <model>.others file to assign values for ntopics, alpha, beta, etc.
	if (read_and_parse(pmodel->dir + pmodel->model_name + pmodel->others_suffix, pmodel)) {
	    return 1;
	}
    }
    
    if (model_status == MODEL_STATUS_UNKNOWN) {
	printf("Please specify the task you would like to perform (-est/-estc/-inf)!\n");
	return 1;
    }
    
    return 0;
}

int utils::read_and_parse(string filename, model * pmodel) {
    // open file <model>.others to read:
    // alpha=?
    // beta=?
    // tau=?
    // ntopics=?
    // ndocs=?
    // nwords=?
    // citer=? // current iteration (when the model was saved)
    
    FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
	printf("Cannot open file: %s\n", filename.c_str());
	return 1;
    }
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    while (fgets(buff, BUFF_SIZE_SHORT - 1, fin)) {
	line = buff;
	strtokenizer strtok(line, "= \t\r\n");
	int count = strtok.count_tokens();
	
	if (count != 2) {
	    // invalid, ignore this line
	    continue;
	}

	string optstr = strtok.token(0);
	string optval = strtok.token(1);
	
	if (optstr == "alpha") {
	    pmodel->alpha_init = atof(optval.c_str());
	    
	} else if (optstr == "beta") {	    
	    pmodel->beta_init = atof(optval.c_str());
	
	} else if (optstr == "tau") {	    
	    pmodel->tau_init= atof(optval.c_str());
	
	} else if (optstr == "ntopics") {
	    pmodel->K = atoi(optval.c_str());
	
	} else if (optstr == "ndocs") {	   
	    pmodel->M = atoi(optval.c_str());
	 
	} else if (optstr == "nwords") {
	    pmodel->V = atoi(optval.c_str());
	
	} else if (optstr == "liter") {
	    pmodel->liter = atoi(optval.c_str());
	
	} else {
	    // any more?
	}
    }
    
    fclose(fin);
    
    return 0;
}

string utils::generate_model_name(int iter) {
    string model_name = "model-";

    char buff[BUFF_SIZE_SHORT];
    
    if (0 <= iter && iter < 10) {
	sprintf(buff, "0000%d", iter);
    } else if (10 <= iter && iter < 100) {
	sprintf(buff, "000%d", iter);
    } else if (100 <= iter && iter < 1000) {
	sprintf(buff, "00%d", iter);
    } else if (1000 <= iter && iter < 10000) {
	sprintf(buff, "0%d", iter);
    } else {
	sprintf(buff, "%d", iter);
    }
    
    if (iter >= 0) {
	model_name += buff;
    } else {
	model_name += "final";
    }
    
    return model_name;
}

void utils::sort(vector<double> & probs, vector<int> & words) {
    for (int i = 0; i < probs.size() - 1; i++) {
	for (int j = i + 1; j < probs.size(); j++) {
	    if (probs[i] < probs[j]) {
		double tempprob = probs[i];
		int tempword = words[i];
		probs[i] = probs[j];
		words[i] = words[j];
		probs[j] = tempprob;
		words[j] = tempword;
	    }
	}
    }
}

void utils::quicksort(vector<pair<int, double> > & vect, int left, int right) {
    int l_hold, r_hold;
    pair<int, double> pivot;
    
    l_hold = left;
    r_hold = right;    
    int pivotidx = left;
    pivot = vect[pivotidx];

    while (left < right) {
	while (vect[right].second <= pivot.second && left < right) {
	    right--;
	}
	if (left != right) {
	    vect[left] = vect[right];
	    left++;
	}
	while (vect[left].second >= pivot.second && left < right) {
	    left++;
	}
	if (left != right) {
	    vect[right] = vect[left];
	    right--;
	}
    }

    vect[left] = pivot;
    pivotidx = left;
    left = l_hold;
    right = r_hold;
    
    if (left < pivotidx) {
	quicksort(vect, left, pivotidx - 1);
    }
    if (right > pivotidx) {
	quicksort(vect, pivotidx + 1, right);
    }    
}


/*
 * given log(a) and log(b), return log(a + b)
 *
 */

double utils::log_sum(double log_a, double log_b)
{
  double v;

  if (log_a < log_b)
  {
      v = log_b+log(1 + exp(log_a-log_b));
  }
  else
  {
      v = log_a+log(1 + exp(log_b-log_a));
  }
  return(v);
}

 /**
   * Proc to calculate the value of the trigamma, the second
   * derivative of the loggamma function. Accepts positive matrices.
   * From Abromowitz and Stegun.  Uses formulas 6.4.11 and 6.4.12 with
   * recurrence formula 6.4.6.  Each requires workspace at least 5
   * times the size of X.
   *
   **/

double utils::trigamma(double x)
{
    double p;
    int i;

    x=x+6;
    p=1/(x*x);
    p=(((((0.075757575757576*p-0.033333333333333)*p+0.0238095238095238)
         *p-0.033333333333333)*p+0.166666666666667)*p+1)/x+0.5*p;
    for (i=0; i<6 ;i++)
    {
        x=x-1;
        p=1/(x*x)+p;
    }
    return(p);
}


/*
 * taylor approximation of first derivative of the log gamma function
 *
 */

double utils::digamma(double x)
{
    double p;
    x=x+6;
    p=1/(x*x);
    p=(((0.004166666666667*p-0.003968253986254)*p+
	0.008333333333333)*p-0.083333333333333)*p;
    p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
    return p;
}


double utils::log_gamma(double x)
{
     double z=1/(x*x);

    x=x+6;
    z=(((-0.000595238095238*z+0.000793650793651)
	*z-0.002777777777778)*z+0.083333333333333)/x;
    z=(x-0.5)*log(x)-x+0.918938533204673+z-log(x-1)-
	log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
    return z;
}



/*
 * make directory
 *
 */

void utils::make_directory(char* name)
{
    mkdir(name, S_IRUSR|S_IWUSR|S_IXUSR);
}


/*
 * argmax
 *
 */

int utils::argmax(double* x, int n)
{
    int i;
    double max = x[0];
    int argmax = 0;
    for (i = 1; i < n; i++)
    {
        if (x[i] > max)
        {
            max = x[i];
            argmax = i;
        }
    }
    return(argmax);
}
