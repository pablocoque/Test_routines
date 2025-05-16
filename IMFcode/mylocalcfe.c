#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>

#define GRAVITY 6.6738e-8
#define SOLAR_MASS 1.989e33
#define CLIGHT_REAL 2.99792458e10
#define STEFAN_BOLTZMANN 5.67374e-5
#define PARSEC 3.085678e18
#define SEC_PER_MEGAYEAR 3.15576e13

double f_tff(double rho){
  return sqrt(3.*M_PI/32./(GRAVITY)/rho); /* free fall time*/
}

double f_mach(double sigmaloc, double csloc){
  return sigmaloc/csloc; /* Mach number*/
}

double f_sigrho(double mach, double beta0){
  double b = 0.5;
  return sqrt(log(1. + 3. * pow(b, 2.) * pow(mach, 2.) * beta0/(beta0 + 1.))); /* dispersion */
}

double f_xcrit(double qvir, double mach){
  double phix = 1.12;
  return pow(M_PI, 2.) * pow(phix, 2.)/15. * qvir * pow(mach, 2.);
}

double f_sfrff(double qvir, double mach, double beta0, double ecore, double sflaw){
  double phit = 1.91;
  double xcrit, sigrho, f;
  
  xcrit = f_xcrit(qvir, mach);
  sigrho = f_sigrho(mach, beta0);
  if (sflaw == 1) {
    f = .5 * ecore/phit * (1. + erf((-2.* log(xcrit) + pow(sigrho, 2.))/(pow(2., 1.5) * sigrho))); /* specific star formation rate per free-fall time*/
  } else {
    f = 0.012;
  }
  return f;
}

double f_dpdx(double x, double mulnx, double sig){
  return 1./(sqrt(2.*M_PI*pow(sig, 2.))*x)*exp(-.5*pow(log(x)-mulnx, 2.)/pow(sig, 2.)); /* overdensity PDF */
}

double f_surfg(double rholoc, double sigmaloc){
  double phiP = 3.;
  return sqrt(2.*rholoc*pow(sigmaloc, 2.)/(M_PI*GRAVITY*phiP)); /* gas surface density */
}

double f_fstar(double rholoc, double sigmaloc, double csloc, double x, double ecore, double beta0, double qvir, double tsn, double tview, double surfGMC, int sflaw, int radfb){
  double phifb=1.6e-5; /* feedback efficiency */
  double kappa0=2.4e-5; /* opacity constant */
  double psi=.3; /* light-to-mass ratio */
  double phitrap=.2; /* trapping ratio */
  double surfg,surffb,mach,sfrff,rhog,tff,efb,einc,efbrad,fstar;
  double epsarr[4];

  surfg = f_surfg(rholoc,sigmaloc); /*!estimate of gas surface density */
  surffb = fmax(surfGMC,surfg); /*!surface density on which radiative feedback acts */
  mach = f_mach(sigmaloc,csloc); /*!Mach number, see above */
  sfrff = f_sfrff(qvir,mach,beta0,ecore,sflaw); /*!specific star formation rate per free-fall time, see above */
  rhog = x*rholoc; /*!gas volume density */
  tff = f_tff(rhog); /*!free-fall time, see above */

  if (radfb == 0 || radfb == 2){
    efb = 0.5*sfrff*tsn/tff*(1.+sqrt(1.+4.*tff*pow(sigmaloc, 2.)/(phifb*sfrff*pow(tsn, 2.)*x))); /* SN feedback*/
  } else {
    efb = 1.;
  }

  einc = sfrff*tview/tff;

  if (radfb > 0){
    efbrad = 2. * STEFAN_BOLTZMANN/(phitrap*pow(kappa0, 2.)*psi*pow(surffb, 3.))
    *(sqrt(1. + 2.*M_PI*CLIGHT_REAL*GRAVITY*phitrap*pow(kappa0, 2.)*pow(surffb, 4.)/(1.*STEFAN_BOLTZMANN))-1.); /* radiative feedback*/
  } else {
    efbrad = 1. ;
  }

  epsarr[0] = ecore; /* SFEs for maximum, SNfeedback, radiative feedback, incomplete */
  epsarr[1] = efb; /* SFEs for maximum, SNfeedback, radiative feedback, incomplete */
  epsarr[2] = efbrad; /* SFEs for maximum, SNfeedback, radiative feedback, incomplete */
  epsarr[3] = einc; /* SFEs for maximum, SNfeedback, radiative feedback, incomplete */
  fstar = fmin(epsarr[0], fmin(epsarr[1], fmin(epsarr[2], epsarr[3])));
  return fstar;
}

double f_integrate(double xsurv, double mulnx, double sig, double rholoc, double sigmaloc, double csloc, double ecore, double beta0, double qvir, double tsn, double tview, double surfGMC, int cce, int sflaw, int radfb){
  double xmin1, xmin2, xmax, f1, dx, xg, fstar, bound, integral, dpdx, f2, frac;
  double xarr[1000];

  xmin1 = exp(mulnx - 5.*sig);
  xmax  = exp(mulnx + 10.*sig);
  if (cce > 0 && xsurv < xmin1){
    xmin1 = xsurv;
  }
  if (cce > 0 && xsurv > xmax){
    xmax = xsurv;
  }

  for(int ix = 0; ix < 1000; ix++){
    xarr[ix] = xmin1*pow(xmax/xmin1, (ix + 0.5)/1000.);
  }
  
  f1 = 0.;
  for(int ix = 0; ix < 1000; ix++){
    dx = xarr[ix]*(pow(xmax/xmin1,1./(2.*1000.))-pow(xmax/xmin1,-1./(2.*1000.))); /* step size */
    xg = xarr[ix]; /*overdensity*/
    fstar = f_fstar(rholoc,sigmaloc,csloc,xg,ecore,beta0,qvir,tsn,tview,surfGMC,sflaw,radfb); /* local SFE */
    bound = fstar/ecore; /*local bound fraction*/
    if(cce == 0){
      bound = 1.; /*if not calculating the cruel cradle effect but the naturally bound fraction of SF, the denominator should contain all SF*/
    }
    if(cce == 2){
      bound = 1.; /*if calculating the cruel cradle effect with respect to all SF, the denominator should contain all SF*/
    }
    integral = bound*fstar*xarr[ix]; /*!integral part 1*/
    dpdx = f_dpdx(xg,mulnx,sig); /*!overdensity PDF, i.e. integral part 2 */
    f1 = f1+integral*dpdx*dx; /* !integral */
  }

  if (cce > 0){
    xmin2 = xsurv; /*if calculating the cruel cradle effect set minimum overdensity to critical overdensity*/
  } else {
    xmin2 = xmin1;
  }

  for(int ix = 0; ix < 1000; ix++){
    xarr[ix] = xmin2*pow(xmax/xmin2,(ix + 0.5)/1000.);
  }

  f2 = 0.;
  for(int ix = 0; ix < 1000; ix++){
    dx = xarr[ix]*(pow(xmax/xmin2,1./(2.*1000.))-pow(xmax/xmin2,-1./(2.*1000.))); /* step size */
    xg = xarr[ix]; /*overdensity*/
    fstar = f_fstar(rholoc,sigmaloc,csloc,xg,ecore,beta0,qvir,tsn,tview,surfGMC,sflaw,radfb); /* local SFE */
    bound = fstar/ecore; /*local bound fraction*/
    if(cce == 2){
      bound = 1.; /*if calculating the cruel cradle effect with respect to all SF, the denominator should contain all SF*/
    }
    integral=bound*fstar*xarr[ix]; /*!integral part 2*/
    dpdx = f_dpdx(xg,mulnx,sig); /*!overdensity PDF, i.e. integral part 2 */
    f2 = f2+integral*dpdx*dx; /* !integral */
  }

  if (f1 == 0.){
    frac = 0.;
  } else {
    frac = f2/f1;
  }
  return frac;
}

double cfe_local_mode(double rholoc, double sigmaloc, double csloc)
{
  double surfg,mach,sigrho,mulnx,fbound;
  
  /* Model parameters */
  int sflaw = 0; /*!star formation law - NOTE: set to 0 for Elmegreen(2002) and to 1 for Krumholz & McKee (2005)*/
  double qvir = 1.3; /*!giant molecular cloud virial parameter*/
  double tsn = 3.; /*!time of the first supernova*/
  double tview = 10.; /*!time at which CFE is determined*/
  double surfGMC = 100;/*!giant molecular cloud surface density*/
  double ecore = 0.5; /*!maximum (protostellar core) star formation efficiency*/
  double beta0 = 1e10;/*!if turbulent-to-magnetic pressure ratio is not specified, set to turbulent-only*/
  int radfb = 0; /*!SN radiative feedback mode - NOTE: set to 0 for supernovae only, to 1 for radiative only, and to 2 for both*/

  tsn = tsn*SEC_PER_MEGAYEAR;
  tview = tview*SEC_PER_MEGAYEAR;
  surfGMC = surfGMC*SOLAR_MASS/pow(PARSEC, 2.);
  surfg = f_surfg(rholoc,sigmaloc); /*!estimate of gas surface density*/
  if(surfg > surfGMC){
    surfGMC=surfg;
  }

  /*!CALCULATE DERIVED PARAMETERS*/
  mach = f_mach(sigmaloc,csloc); /*!Mach number*/
  sigrho = f_sigrho(mach,beta0); /*!dispersion of overdensity PDF*/
  mulnx = -.5*pow(sigrho,2.); /*!logarithmic mean of overdensity PDF*/
  /*!CALCULATE F_BOUND*/
  fbound = f_integrate(0.,mulnx,sigrho,rholoc,sigmaloc,csloc,ecore,beta0,qvir,tsn,tview,surfGMC,0,sflaw,radfb); /*!naturally bound part of star formation*/
  
  return fbound; /* CFE local without Cruddle effect*/
}

double ICMF(double mass, double *param){
  // double *arr = (double *)param;
  double mtrunc = param[0];
  double norm_ctn = param[1];
  return norm_ctn *exp(-mass/mtrunc)/pow(mass, 2);
}

double ICMF_mPDF(double mass, double *param){
  // double *arr = (double *)params;
  double mtrunc = param[0];
  double norm_ctn = param[1];
  return norm_ctn *exp(-mass/mtrunc)/mass;
}

// double trapezoidal(double min, double max, void *params, double (*func)(double, void*), int n){
//   double sum = 0., integral;
//   double h, x;
//   h = fabs(max - min)/n;
//   for(int i = 1; i<n; i++){
//     x = min + i*h;
//     sum = sum + (*func)(x, params);
//   }
//   integral = 0.5*h*((*func)(min, params) + (*func)(max, params) + 2*sum);
//   return integral;
  // h=fabs(max-min)/n;
  // for(int i=1;i<n;i++){
  //   x=min+i*h;
  //   if(i%2==0){
  //     sum=sum+2*(*func)(x, params);
  //   }
  //   else{
  //     sum=sum+4*(*func)(x, params);
  //   }
  // }
  // integral=(h/3)*((*func)(min, params)+(*func)(max, params)+sum);
  // return integral;
// }

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

// double fraction_trapz_log(double min, double max, double lowl, double highl,void *params, double (*func)(double, void*), int n){
//   double upintegral = 0., downintegral = 0.;
//   double h, dx, low, previous, current;
//   h = (log10(max/SOLAR_MASS) - log10(min/SOLAR_MASS))/n;
//   low = log10(min/SOLAR_MASS);
//   for(int i = 1; i<n; i++){
//     current = pow(10,low + i*h)*SOLAR_MASS;
//     previous = pow(10,low + (i-1)*h)*SOLAR_MASS;
//     dx = fabs(current - previous);
//     if ((low + (i-1)*h)>log10(lowl/SOLAR_MASS) && (low + i*h)<log10(highl/SOLAR_MASS)){
//       upintegral += 0.5*dx*((*func)(current, params) + (*func)(previous, params));
//     }
//     downintegral += 0.5*dx*((*func)(current, params) + (*func)(previous, params));
//   }
//   return upintegral/downintegral;
// }

double envelope(double m){
  return 1/pow(m, 2);
}

double inverse_envelope(double u, double norm, double min){
  return 1/(1/min - u/norm);
}

double draw_mass(double min, double max, double *params, double (*func)(double, double*), gsl_rng *rng, int *n){
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
  // guess1 = (log10(max)- log10(min)) * gsl_rng_uniform(rng) + log10(min);
  // if (log10((*func)(max, params)) > __DBL_MIN_10_EXP__ ){
  //   guess2 = (log10((*func)(min, params))- log10((*func)(max, params))) * gsl_rng_uniform(rng) + log10((*func)(max, params));
  // } else {
  //   guess2 = (log10((*func)(min, params))- __DBL_MIN_10_EXP__ )* gsl_rng_uniform(rng) + __DBL_MIN_10_EXP__;
  // }
  // criterion = log10((*func)(pow(10,guess1), params));
  
  if (guess2 <= criterion){
    return guess1;
  } else {
    // printf("regression called\n");
    *n += 1;
    return draw_mass(min, max, params, (*func), rng, n);
    // return 0.;
  }
}