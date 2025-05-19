#include <gsl/gsl_randist.h>

#define GRAVITY 6.6738e-8
#define SOLAR_MASS 1.989e33
#define CLIGHT_REAL 2.99792458e10
#define STEFAN_BOLTZMANN 5.67374e-5
#define PARSEC 3.085678e18
#define SEC_PER_MEGAYEAR 3.15576e13

double ICMF(double mass, double *param);
double ICMF_mPDF(double mass, double *param);
double trapz_log(double min, double max, double *params, double (*func)(double, double*), int n);
double draw_mass(double min, double max, double *params, double (*func)(double, double*), gsl_rng *rng);
double cfe_local(double rho, double sigma, double cs, double tview);