#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include "proto.h"

static double f_xcrit(double mach);
static double f_sigrho(double mach);
static double f_sfrff(double mach);
static double f_dpdx(double x, double mulnx, double sig);
static double f_surfg(double rholoc, double sigmaloc);
static double f_fstar(double rholoc, double sigmaloc, double csloc, double surfdens, double x, double tview);
static double f_phiad(double x);
static double f_xcce(double sigmaloc, double surfGMC, double tview);
static double f_integrate(double xsurv, double mulnx, double sig, double rholoc, double sigmaloc, double csloc, double surfdens, double tview, int cce);

/* Model parameters */
static double qvir = 1.3;                                // giant molecular cloud virial parameter
static double tsn = 3.*SEC_PER_MEGAYEAR;                 // time of the first supernova
static double surfGMC = 100.*SOLAR_MASS/(PARSEC*PARSEC); // giant molecular cloud surface density
static double ecore = 0.5;                               // maximum (protostellar core) star formation efficiency
static double beta0 = 1e10;                              // if turbulent-to-magnetic pressure ratio is not specified,
                                                         // set to turbulent-only
static int radfb = 1;                                    // SN radiative feedback mode - NOTE: set to 0 for supernovae only,
                                                         // to 1 for radiative only, and to 2 for both
static int sflaw = 1;                                    // star formation law - NOTE: set to 0 for Elmegreen(2002) and 
                                                         // to 1 for Krumholz & McKee (2005)

/*  Functions to obtain CFE based on Kruijsen 2012 model */

/** \brief Critical overdensity for SF
 * 
 *  critical overdensity for SF in the KM05 sSFR_ff as a function of Mach
 *  number and GMC virial ratio, from Krumholz&McKee2005
 * 
 *  \param mach Mach number
*/
static double f_xcrit(double mach){
  double phix = 1.12;
  return M_PI * M_PI * phix * phix /15. * qvir * mach * mach;
}

/** \brief Dispersion of the overdensity PDF
 *  
 *  dispersion of overdensity PDF as a function of Mach number and magnetic
 *  pressure ratio, based on e.g. Padoan&Nordlund2011
 * 
 *  \param mach Mach number
*/
static double f_sigrho(double mach){
  double b = 0.5;
  return sqrt(log(1. + 3. * b * b * mach * mach * beta0/(beta0 + 1.)));
}

/** \brief Specific star formation rate per free-fall time
 * 
 *  It implements two models following Elmegreen2002 or Krumholz&McKee2005
 * 
 *  \param mach Mach number
*/
static double f_sfrff(double mach){
  double phit = 1.91;
  double xcrit, sigrho, f;

  xcrit = f_xcrit(mach);
  sigrho = f_sigrho(mach);
  if (sflaw == 1) {
    f = 0.5 * ecore/phit * (1. + erf((-2.* log(xcrit) + sigrho * sigrho)/(pow(2., 1.5) * sigrho)));
  } else {
    f = 0.012;
  }
  return f;
}

/** \brief Gives the overdensity PDF of the ISM
 *
 *  This is obtained as a function of the overdensity or density contrast x, its
 *  logarithmic mean and its dispersion.
 * 
 *  \param x         density contrast
 *  \param mulnx     logarithmic mean of overdensity PDF
 *  \param sig       dispersion of overdensity PDF
*/
static double f_dpdx(double x, double mulnx, double sig){
  return 1. / (sqrt(2. * M_PI * sig * sig) * x) * exp(-0.5 * pow(log(x) - mulnx, 2.) / (sig * sig));
}

/** \brief Estimate of the local gas surface density
 * 
 *  Function to estimate the local gas surface density from the local density
 *  and velocity dispersion. The estimate is based on rough relation between
 *  gas surface density, volume density and velocity dispersion assuming an
 *  equilibrium disk. Follows the model of Krumholz&McKee2005 (eq 32).
 * 
 *  \param rholoc    local pressure
 *  \param sigmaloc  local velocity dispersion
*/
static double f_surfg(double rholoc, double sigmaloc){
  double phiP = 3.; // Contribution of stars gravity to the pressure of the disk.
  double surfg;

  surfg = sqrt(2. * rholoc * sigmaloc * sigmaloc / (M_PI * GRAVITY * phiP));
  return surfg;
}

/** \brief Obtain local naturally bound fraction of star formation
 * 
 *  Function to calculate the local star formation efficiency for a density contrast
 *  of the ISM. The SFE is obtained following 4 scenarios: maximum sfe, including SN
 *  feedback, including radiative feedback and a incomplete sfe. Then, the lowest SFE
 *  is returned.
 * 
 *  \param rholoc    local pressure
 *  \param sigmaloc  local velocity dispersion
 *  \param csloc     local sound speed
 *  \param surfdens  local estimated surface density
 *  \param x         density contrast
*/
static double f_fstar(double rholoc, double sigmaloc, double csloc, double surfdens, double x, double tview){
  double phifb=0.16; // feedback efficiency in cm2 s-3
  double kappa0=2.4e-4; // opacity constant in cm2 g-1 K-2
  double psi=3e3; // light-to-mass ratio in cm2 s-1
  double phitrap=0.2; // trapping ratio
  double surffb,mach,sfrff,rhog,tff,efb,einc,efbrad,fstar;

  surffb = fmax(surfGMC,surfdens); // surface density on which radiative feedback acts
  mach = sigmaloc/csloc; // Mach number, from Krumholz&McKee2005
  sfrff = f_sfrff(mach); // specific star formation rate per free-fall time
  rhog = x*rholoc; // gas volume density
  tff = sqrt(3.*M_PI/32./GRAVITY/rhog); // free-fall time

  // SFE for SNfeedback
  if (radfb == 0 || radfb == 2){
    efb = 0.5*sfrff*tsn/tff*(1+sqrt(1.+4*tff*sigmaloc*sigmaloc/(phifb*sfrff*tsn*tsn*x)));
  } else {
    efb = 1;
  }

  // SFE for incomplete feedback
  einc = sfrff*tview/tff;

  // SFE for radiative feedback
  if (radfb > 0){
    efbrad = 2. * STEFAN_BOLTZMANN/(phitrap*kappa0*kappa0*psi*pow(surffb, 3.)) *
             (sqrt(1. + 2.*M_PI*CLIGHT_REAL*GRAVITY*phitrap*kappa0*kappa0*pow(surffb, 4.)/STEFAN_BOLTZMANN)-1.);
  } else {
    efbrad = 1. ;
  }

  fstar = fmin(ecore, fmin(efb, fmin(efbrad, einc)));
  return fstar;
}

static double f_phiad(double x){
  double phit = 3.1 * sqrt((qvir/1.3) * (x/1e4));
  return exp(-2. * phit);
}

static double f_xcce(double sigmaloc, double surfGMC, double tview){
  double phiad, xfit, xfit0, diff, diffx;
  int niter = 0;
  int maxiter = 100;
  double xmin = 1e-4;
  double xmax = 1e8;
  double accuracy = 1e-6;
  double xarr[100], xarr2[100];
  
  xfit = xmin*xmin/xmax;
  xfit0 = xfit*accuracy;
  
  while(fabs(log10(xfit/xfit0)) > accuracy){
    xfit0 = xfit;
    for(int ix = 0; ix < 100; ix++){
      xarr[ix] = pow(10. ,(ix/100.)*log10(xmax/xmin) + log10(xmin));
      phiad = f_phiad(xarr[ix]);
      xarr2[ix] = 87.5*sqrt(M_PI)*0.7*1.5*GRAVITY*2.8*surfGMC*phiad*tview/sigmaloc;
    }
    diff = 1e30;
    for(int ix = 0; ix < 100; ix++){
      diffx = fabs(xarr2[ix] - xarr[ix]);
      if(diffx < diff && ix > 0){
        diff = diffx;
        xfit = xarr[ix];
        xmin = xarr[ix-1];
        xmax = xarr[ix+1];
      }
    }

    niter++;
    
    if(niter == maxiter){
      printf("CFE: x_cce not converging, niter: %d, diff: %e \n", niter, diff);
    }
  }

  return xfit;
}

/** \brief Integrate the overdensity PDF
 * 
 *  Function to integrate the overdensity PDF and obtain the fraction of stars born in
 *  bound systems.
 * 
 *  \param mulnx     logarithmic mean of overdensity PDF
 *  \param sig       dispersion of overdensity PDF
 *  \param rholoc    local density
 *  \param sigmaloc  local velocity dispersion
 *  \param csloc     local sound speed
 *  \param surfdens  local estimated surface density
*/
static double f_integrate(double xsurv, double mulnx, double sig, double rholoc, double sigmaloc, double csloc, double surfdens, double tview, int cce){
  double xmin, xmax, f1, dx, xg, bound, integral, dpdx, f2, frac, fstar;

  xmin = exp(mulnx - 5. * sig);
  xmax = exp(mulnx + 10. * sig);

  if (cce == 1 && xsurv < xmin) xmin = xsurv;
  if (cce == 1 && xsurv > xmax) xmax = xsurv;

  f1 = 0.;
  for(int ix = 0; ix < 1000; ix++){
    xg = xmin * pow(xmax / xmin, (ix + 0.5) / 1000.); // overdensity
    dx = xg * (pow(xmax / xmin, 1. / (2. * 1000.)) - pow(xmax / xmin, -1. / (2. * 1000.))); // step size
    fstar = f_fstar(rholoc,sigmaloc,csloc,surfdens,xg, tview); // local SFE
    bound = fstar / ecore; // local bound fraction
    if (cce == 0) bound = 1.; // if not calculating the cruel cradle effect but the naturally bound fraction of SF, the denominator should contain all SF
    integral = bound * fstar * xg; // integral part 1
    dpdx     = f_dpdx(xg, mulnx, sig); // overdensity PDF, i.e. integral part 2 
    f1      += integral * dpdx * dx; // integral 
  }

  if (cce == 1) xmin = xsurv;

  f2 = 0.;
  for(int ix = 0; ix < 1000; ix++){
    xg = xmin * pow(xmax / xmin, (ix + 0.5) / 1000.); // overdensity
    dx = xg * (pow(xmax / xmin, 1. / (2. * 1000.)) - pow(xmax / xmin, -1. / (2. * 1000.))); // step size
    fstar = f_fstar(rholoc,sigmaloc,csloc,surfdens,xg, tview); // local SFE 
    bound = fstar / ecore; // local bound fraction
    if(bound > 1){
      printf("CFE: SFE higher than maximum protostellar core SFE, ratio: %f \n", bound);
    }
    integral = bound * fstar * xg;    // integral part 2
    dpdx     = f_dpdx(xg, mulnx, sig);    // overdensity PDF, i.e. integral part 2 
    f2      += integral * dpdx * dx; // integral 
  }

  if(f1 == 0.){
    frac = 0.;
  } else {
    frac = f2 / f1;
  }
  
  return frac;
}

/** \brief Obtain the cluster formation efficiency for a gas cell
 * 
 *  Given the local properties (density and pressure)of a gas cell where
 *  a star was formed the cluster formation efficiency is calculated. This
 *  is a ported version from the model explained in Kruijssen2012 and
 *  publicly available. The implementation differs on the estimate of the
 *  surface density of the gas.
 * 
 *  \param igas index of the gas cell where star formation happened
*/
double cfe_local(double rholoc, double sigmaloc, double csloc, double tview){
  double mach, sigrho, mulnx, fbound, xcce, fcce, surfdens;

  /* CALCULATE DERIVED PARAMETERS*/
  mach     = sigmaloc / csloc; // Mach number
  sigrho   = f_sigrho(mach); // dispersion of overdensity PDF
  mulnx    = -.5 * pow(sigrho, 2.); // logarithmic mean of overdensity PDF
  surfdens = f_surfg(rholoc, sigmaloc); // Estimate of local gas surface density

  /* CALCULATE CFE*/
  fbound = f_integrate(0., mulnx, sigrho, rholoc, sigmaloc, csloc, surfdens, tview, 0); // naturally bound part of star formation

  xcce = f_xcce(sigmaloc, surfGMC, tview); // x_cce
  fcce = f_integrate(xcce, mulnx, sigrho, rholoc, sigmaloc, csloc, surfdens, tview, 1); // bound part after cruel cradle effect
  
  return fcce; /* CFE local after Cruddle Craddle effect*/
}
/* End of functions to obtain CFE based on Kruijsen 2012 model */