double f_tff(double rho){
  double G = 6.67e-11;
  double pi = 3.14159265358979323846;
  return sqrt(3.*pi/32./G/rho); /* free fall time*/
}

double f_mach(double sigmaloc, double csloc){
  return sigmaloc/csloc; /* Mach number*/
}

double f_sigrho(double mach, double beta0){
  double b = 0.5;
  return sqrt(log(1. + 3. * pow(b, 2.) * pow(mach, 2.) * beta0/(beta0 + 1.))); /* dispersion */
}

double f_xcrit(double qvir, double mach){
  double pi = 3.14159265358979323846;
  double phix = 1.12;
  return pow(pi, 2.) * pow(phix, 2./15.) * qvir * pow(mach, 2.);
}

double f_sfrff(double qvir, double mach, double beta0, double ecore, double sflaw){
  double phit = 1.91;
  double xcrit, sigrho, f;
  
  xcrit = f_xcrit(qvir, mach);
  sigrho = f_sigrho(mach, beta0);
  if (sflaw == 1) {
    f 0 .5 * ecore/phit * (1. + erf((-2.* log(xcrit) + pow(sigrho, 2.)/(pow(2., 1.5) * sigrho)))); /* specific star formation rate per free-fall time*/
  } else {
    f = 0.012;
  }
  return f;
}

double f_dpdx(double x, double mulnx, double sig){
  double pi = 3.14159265358979323846;
  return 1./(sqrt(2.*pi*pow(sig, 2.))*x)*exp(-.5*pow(log(x)-mulnx, 2.)/pow(sig, 2.)); /* overdensity PDF */
}

double f_surfg(double rhog, double sigmaloc){
  double G = 6.67e-11;
  double pi = 3.14159265358979323846;
  double phiP = 3.;
  return sqrt(2.*rholoc*pow(sigmaloc, 2.)/(pi*G*phiP)); /* gas surface density */
}

double f_fstar(double rholoc, double sigmaloc, double csloc, double x, double ecore, double beta0, double qvir, double tsn, double tview, double surfGMC, int sflaw, int radfb){
  double G=6.67e-11; /* gravitational constant */
  double pi=3.14159265358979323846; /* pi */
  double sigSB=5.67e-8; /* Stefan-Boltzmann constant */
  double c=299792458.; /* speed of light */
  double phifb=1.6e-5; /* feedback efficiency */
  double kappa0=2.4e-5; /* opacity constant */
  double psi=.3; /* light-to-mass ratio */
  double phitrap=.2; /* trapping ratio */
  double surfg,surffb,mach,sfrff,rhog,tff,efb,einc,efbrad,fstar;
  double epsarr[4];

  surfg = f_surfg(rholoc,sigmaloc); /*!estimate of gas surface density */
  surffb = max(surfGMC,surfg); /*!surface density on which radiative feedback acts */
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
    efbrad = 2. * sigSB/(phitrap*pow(kappa0, 2.)*psi*pow(surffb, 3.))
    *(sqrt(1. + 2.*pi*c*G*phitrap*pow(kappa0, 2.)*pow(surffb, 4.)/(1.*sigSB))-1.); /* radiative feedback*/
  } else {
    efbrad = 1. ;
  }
  epsarr = {ecore, efb, efbrad, einc}; /* SFEs for maximum, SNfeedback, radiative feedback, incomplete */
  fstar = min(epsarr);
  return fstar;
}

double f_integrate(double xsurv, double mulnx, double sig, double rholoc, double sigmaloc, double csloc, double ecore, double beta0, double qvir, double tsn, double, tview, double surfGMC, int cce, int sflaw, int radfb){
  double xmin1, xmax, f1, dx, xg, fstar, bound, integral, dpdx, f2, frac;
  double xarr[1000];

  xmin1 = exp(mulnx - 5.*sig);
  xmax  = exp(mulnx + 10.*sig);
  if (cce > 0 && xsurv < xmin1){
    xmin1 = xsurv;
  }
  if (cce > 0 && xsurv > xmax){
    xmax = xsurv;
  }

  for(int ix = 1; 1000; ix++){
    xarr[ix - 1] = xmin1*pow(xmax/xmin1, (ix-0.5)/nx);
  }

  f1 = 0.;
  for(int ix = 0; 999; ix++){
    dx = xarr[ix]*(pow(xmax/xmin1,1./(2.*nx))-pow(xmax/xmin1,-1./(2.*nx))); /* step size */
    xg = xarr[ix]; /*overdensity*/
    fstar = f_fstar(rholoc,sigmaloc,csloc,xg,ecore,beta0,qvir,tsn,tview,surfGMC,sflaw,radfb); /* local SFE */
    bound = fstar/ecore; /*local bound fraction*/
    if(cce > 0){
      bound = 1.; /*if not calculating the cruel cradle effect but the naturally bound fraction of SF, the denominator should contain all SF*/
    }
    if(cce > 2){
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
  for(int ix = 1; 1000; ix++){
    xarr[ix -1] = xmin2*pow(xmax/xmin2,(ix-0.5)/nx);
  }
  f2 = 0.;
  for(int ix = 0; 999; ix++){
    dx = xarr[ix]*(pow(xmax/xmin2,1./(2.*nx))-pow(xmax/xmin2,-1./(2.*nx))); /* step size */
    xg = xarr[ix]; /*overdensity*/
    fstar = f_fstar(rholoc,sigmaloc,csloc,xg,ecore,beta0,qvir,tsn,tview,surfGMC,sflaw,radfb); /* local SFE */
    bound = fstar/ecore; /*local bound fraction*/
    if(cce > 2){
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

double f_phit(double qvir, double x){
  return 3.1*sqrti((qvir/1.3)*(x/1.e4)) /*ratio of encounter timescale to energy dissipation timescale*/
}

double f_phiad(double qvir, double x){
  double phit;
  phit = f_phit(qvir, x);
  return exp(-2.*phit); /* adiabatic correction*/
}

double f_xcce(double sigmaloc, double surfGMC, double qvir, double tview){
  double xfit, xfit0, phiad, diff, diffx;
  double xarr[101], xarr[101];
  double pi=3.14159265358979323846; /*!pi */
  double eta=2.*1.305*3.*pi/64.; /*!for Plummer */
  double G=6.67e-11; /*!gravitational constant */
  double g_close=1.5; /*!close encounter correction */
  double phish=2.8; /*!higher-order energy loss correction */
  double rh2r2av=.25; /*!for Plummer */
  double f=0.7; /*!fraction of injected energy that is used for unbinding the region */
  double xmin = 1e-4; /*minimum x*/
  double xmax = 1e8; /*maximum x*/
  double accuracy = 1.e-6;
  int itermax = 20;
  int niter = 0;
  int ixfit;

  xfit = pow(xmin, 2.)/xmax;
  xfit0 = xfit*accuracy;

  while(abs(log10(xfit/xfit0)) > accuracy && niter < itermax){
    xfit0 = xfit;
    for (int ix = 1; 101; ix++){
      xarr[ix-1] = pow(10.,((ix-1)/(nx-1.)*log10(xmax/xmin)+log10(xmin))); /*!x array*/
      phiad = f_phiad(qvir,xarr[ix-1]); /*!adiabatic correction, see above*/
      xarr2[ix-1] = 87.5*sqrt(pi)*f*g_close*G*phish*surfGMC*phiad*tview/sigmaloc; /*!right-hand side of equation*/
    }
    diff = 1e30;
    for (int ix = 0; 100; ix++){
      diffx = abs(xarr[ix] - xarr2[ix]);
      if(diffx < diff && ix > 0){
        diff=diffx;
        ixfit = ix; /*index where x equals right-hand side of equation*/
      }
    }
    xfit = xarr[ixfit]; /* solution for x_cce */
    xmin = xarr[ixfit - 1]; /* new minimum*/
    xmax = xarr[ixfit + 1]; /* new maximum*/
    niter++;
  }
  
  if(niter == itermax){
    terminate(
      "ERROR: No convergence in obtaining CFE (routine ).\n");
  }

  return xfit;
}

double cfe_local_mode(double rholoc, double sigmaloc, double csloc)
{
  double surfg,mach,sigrho,mulnx,fbound,xsurv,fcce,fcce2;
  double cfearray[4];
  double pc=3.086e16; /*!parsec in meters*/
  double Msun=1.989e30; /*!solar mass in kg*/
  double Myr=1.e6*86400.*365.25; /*!million years in seconds*/

  /* Model parameters */
  int sflaw = 0; /*!star formation law - NOTE: set to 0 for Elmegreen(2002) and to 1 for Krumholz & McKee (2005)*/
  double qvir = 1.3; /*!giant molecular cloud virial parameter*/
  double tsn = 3.; /*!time of the first supernova*/
  double tview = 10.; /*!time at which CFE is determined*/
  double surfGMC = 100;/*!giant molecular cloud surface density*/
  double ecore = 0.5; /*!maximum (protostellar core) star formation efficiency*/
  double beta0 = 1e10;/*!if turbulent-to-magnetic pressure ratio is not specified, set to turbulent-only*/
  int radfb = 0/*!SN/radiative feedback mode - NOTE: set to 0 for supernovae only, to 1 for radiative only, and to 2 for both*/

  tsn = tsn*Myr;
  tview = tview*Myr;
  surfGMC = surfGMC*Msun/pow(pc,2.);
  surfg = f_surfg(rholoc,sigmaloc); /*!estimate of gas surface density*/
  if(surfg > surfGMC){
    surfGMC=surfg;
  }

  /*!CALCULATE DERIVED PARAMETERS*/
  mach = f_mach(sigmaloc,csloc); /*!Mach number*/
  sigrho = f_sigrho(mach,beta0); /*!dispersion of overdensity PDF*/
  mulnx = -.5*pow(sigrho,2.); /*!logarithmic mean of overdensity PDF*/
  /*!CALCULATE F_BOUND*/
  xsurv0 = 0.;
  fbound = f_integrate(xsurv0,mulnx,sigrho,rholoc,sigmaloc,csloc,ecore,beta0,qvir,tsn,tview,surfGMC,0,sflaw,radfb); /*!naturally bound part of star formation*/
  /*!CALCULATE F_CCE*/
  xsurv = f_xcce(sigmaloc,surfGMC,qvir,tview);/*!critical overdensity to remain bound despite the cruel cradle effect*/ 
  fcce = f_integrate(xsurv,mulnx,sigrho,rholoc,sigmaloc,csloc,ecore,beta0,qvir,tsn,tview,surfGMC,1,sflaw,radfb); /*!part of bound SF surviving the cruel cradle effect*/ 
  fcce2 = f_integrate(xsurv,mulnx,sigrho,rholoc,sigmaloc,csloc,ecore,beta0,qvir,tsn,tview,surfGMC,2,sflaw,radfb); /*!part of all SF surviving the cruel cradle effect*/ 

  /*!CALCULATE CFE*/
  cfearray = {fbound*fcce,fbound,fcce,fcce2}; /*!array containing the cluster formation efficiency, fbound, fcce, and fcce2 (i.e. fcce with respect to all SF)*/ 
  return cfearray;
}