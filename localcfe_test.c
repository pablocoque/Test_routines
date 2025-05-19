#include "proto.h"
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(){
    gsl_rng *random_generator;
    random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
    double rhoMW = 0.03 * SOLAR_MASS/pow(PARSEC, 3); // volume density in Msun/pc^3
    double rhohigh = 0.3 * SOLAR_MASS/pow(PARSEC, 3); // volume density in Msun/pc^3
    double sigmaMW = 7.e5;
    double csMW = 0.2e5;
    double mGMC = 2e9*SOLAR_MASS; // in Msun units
    double sfe = 0.1;
    double mstar = 5e5*SOLAR_MASS; // in Msun units
    
    double cfe1, cfe2, mtrunc, norm_ctn, meanmc, Nclus, mclus, f, totclusmass;
    int Nreal, N;
    double mmin = 1e2*SOLAR_MASS;
    double mmax = 1e8*SOLAR_MASS;
    double lowlimit = 5e3*SOLAR_MASS;
    double params[2];
    double tview;
    double min_tview = 1e-2*SEC_PER_MEGAYEAR;
    double max_tview = 100.*SEC_PER_MEGAYEAR;
    
    double tview_step = (log10(max_tview) - log10(min_tview))/999.;
    for (int i=0; i<1000; i++){
        tview = pow(10,log10(min_tview) + i*tview_step);
        cfe1 = cfe_local(rhoMW, sigmaMW, csMW, tview);
        cfe2 = cfe_local(rhohigh, sigmaMW, csMW, tview);
        printf("%.5e %.5e %.5e\n", tview, cfe1, cfe2);
    }
    
    // mtrunc = cfe*sfe*mGMC;
    // printf("Truncation mass: %e Msun. \n", mtrunc/SOLAR_MASS);
    // params[0] = mtrunc;
    // params[1] = 1.;

    // norm_ctn = 1/trapz_log(mmin, mmax, params, ICMF, 1000);
    // f = norm_ctn * trapz_log(lowlimit, fmin(mstar,mmax), params, ICMF, 1000);
    // params[1] = norm_ctn;
    // printf("Expected fraction %f \n", f);
    // meanmc = trapz_log(mmin, mmax, params, ICMF_mPDF, 1000);

    // printf("Mean cluster mass= %e MSun. \n", meanmc/SOLAR_MASS);

    // Nclus = cfe*mstar/meanmc;
    // printf("Expected number of clusters=%f \n", Nclus);
    // for (int j=0; j<10;j++){
    //     totclusmass = 0.;
    //     N = 0;
    //     Nreal = gsl_ran_poisson(random_generator, Nclus);
    //     printf("Actual number of cluster masses to be sampled=%d \n", Nreal);
        
    //     for(int i=1; i<=Nreal; i++){
    //         mclus = draw_mass(mmin, mmax, params, ICMF, random_generator);
    //         if (mclus > lowlimit && totclusmass + mclus < mstar){
    //             // printf("Cluster %d has mass %e Msun.\n", i, mclus);
    //             totclusmass += mclus;
    //             N += 1;
    //         } else {
    //             continue;
    //         }
    //     }

    //     printf("TEST: Nclus*mc = %e. \n", Nreal*meanmc/SOLAR_MASS);
    //     printf("TEST: CFE*Mstar = %e. \n", cfe*mstar/SOLAR_MASS);

    //     printf("Total mass in clusters %e Msun in %d clusters from star particle of mass %e Msun.\n", totclusmass/SOLAR_MASS, N, mstar/SOLAR_MASS);
    //     printf("Expected fraction of sampled masses %f \n", f);
    //     printf("Obtained fraction of sampled masses %f \n", 1.0*N/Nreal);
    //     printf("%d, %d, %f \n", Nreal, N, 1.0*N/Nreal);
    // }

    return 0;
}