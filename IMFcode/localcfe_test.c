#include "mylocalcfe.c"

int main(){
    gsl_rng *random_generator;
    random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
    double rhoMW = 6.804626e-25;
    double sigmaMW = 7.739466e5;
    double csMW = 9.991607e5;
    // double mGMC = 0.009*1.6e8; // in Msun units
    double factor = 1e-10;
    double mGMC = 7e-2; // in 1e10 Msun units
    double sfe = 0.1;
    double mstar = 6.6e-2;
    
    double cfe, mtrunc, norm_ctn, meanmc, Nclus, mclus, f, totclusmass, mdrawn;
    int Nreal;
    int Nreg = 0, Ntot = 0, N = 0;
    double mmin = 1e2*factor;
    double mmax = 1e8*factor;
    double lowlimit = 5e3*factor;
    double params[2];
    // double (*fun_ptr1)(double, void *) = &ICMF;
    // double (*fun_ptr2)(double, void *) = &ICMF_mPDF;

    cfe = cfe_local_mode(rhoMW, sigmaMW, csMW);
    printf("TEST: Obtained CFE local mode is %g. \n", cfe);
    
    printf("TESTING WITH mGMC %e and Mstar %e \n", mGMC, mstar);

    mtrunc = cfe*sfe*mGMC;
    params[0] = mtrunc;
    params[1] = 1.;

    norm_ctn = 1/trapz_log(mmin, mmax, params, &ICMF, 1000);
    // printf("Norm: %f \n", norm_ctn);
    f = norm_ctn * trapz_log(lowlimit, fmin(mstar,mmax), params, &ICMF, 1000);
    params[1] *= norm_ctn;
    printf("ICMF params: %e %e. \n", params[0], params[1]);
    printf("Expected fraction %f \n", f);

    meanmc = trapz_log(mmin, mmax, params, &ICMF_mPDF, 1000);
    printf("Mean cluster mass= %e MSun. \n", meanmc);

    Nclus = cfe*mstar/meanmc;
    printf("Expected number of clusters=%f \n", Nclus);
    for (int j=0; j<1000;j++){
        totclusmass = 0.;
        mdrawn = 0.;
        Nreal = gsl_ran_poisson(random_generator, Nclus);
        Ntot += Nreal;
        // printf("Actual number of cluster masses to be sampled=%d \n", Nreal);
        
        for(int i=1; i<=Nreal; i++){
            mclus = draw_mass(mmin, mmax, params, &ICMF, random_generator, &Nreg);
            mdrawn += mclus;

            if (mclus > lowlimit){
                // printf("Cluster %d has mass %e Msun.\n", i, mclus);
                totclusmass += mclus;
                N += 1;
            } else {
                continue;
            }
        }

        // printf("Total mass in clusters %e Msun in %d clusters from star particle of mass %e Msun.\n", totclusmass, N, mstar);
        // printf("%d, %d, %f, %e, %e \n", Nreal, N, 1.0*N/Nreal, mdrawn, totclusmass);
    }

    printf("Sampling Ntot %d, accepted %d \n", Ntot, N);
    printf("Regression called %d times \n", Nreg);

    gsl_rng_free(random_generator);
    return 0;
}