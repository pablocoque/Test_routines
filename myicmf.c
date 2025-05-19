#include <math.h>
#include <gsl/gsl_randist.h>
#include "proto.h"

static double envelope(double m);
static double inverse_envelope(double u, double norm, double min);

double ICMF(double mass, double *param){
    double mtrunc = param[0];
    double norm_ctn = param[1];
    return norm_ctn *exp(-mass/mtrunc)/pow(mass, 2);
}
  
double ICMF_mPDF(double mass, double *param){
    double mtrunc = param[0];
    double norm_ctn = param[1];
    return norm_ctn *exp(-mass/mtrunc)/mass;
}
  
double trapz_log(double min, double max, double *params, double (*func)(double, double*), int n){
    double integral = 0.;
    double h, dx, low, previous, current;
    h = (log10(max) - log10(min))/n;
    low = log10(min);
    for(int i = 1; i<n; i++){
        current = pow(10,low + i*h);
        previous = pow(10,low + (i-1)*h);
        dx = fabs(current - previous);
        integral += 0.5*dx*((*func)(current, params) + (*func)(previous, params));
    }
    return integral;
}
  
double draw_mass(double min, double max, double *params, double (*func)(double, double*), gsl_rng *rng){
    double u, guess1, guess2, criterion, env_norm, scale;
    
    u = gsl_rng_uniform(rng);
    env_norm = 1/(1/min - 1/max);
    guess1 = inverse_envelope(u, env_norm, min);
  
    guess2 = gsl_rng_uniform(rng);
    if ((*func)(min, params) > env_norm * envelope(min)){
        scale = (*func)(min, params)/(env_norm * envelope(min));
    } else {
        scale = 1.;
    }
    
    criterion = (*func)(guess1, params)/(scale * env_norm * envelope(guess1));
    
    if (guess2 <= criterion){
        return guess1;
    } else {
        return draw_mass(min, max, params, (*func), rng);
    }
}

static double envelope(double m){
    return 1/pow(m, 2);
}

static double inverse_envelope(double u, double norm, double min){
    return 1/(1/min - u/norm);
}