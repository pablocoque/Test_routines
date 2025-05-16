#include "mylocalcfe.c"

int main(){
    gsl_rng *random_generator;
    random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
    double rhoMW = 1.154e-24;
    double sigmaMW = 1.3079e6;
    double csMW = 1.6885e6;
    double mGMC = 2e9; // in Msun units
    double sfe = 0.1;
    double mstar = 2e8;
    
    double cfe, mtrunc, norm_ctn, meanmc, Nclus, mclus, f, totclusmass;
    int Nreal, N;
    double mmin = 1e2;
    double mmax = 1e8;
    double lowlimit = 3e5;
    double params[2];
    double (*fun_ptr1)(double, void *) = &ICMF;
    double (*fun_ptr2)(double, void *) = &ICMF_mPDF;

    cfe = cfe_local_mode(rhoMW, sigmaMW, csMW);
    printf("TEST: Obtained CFE local mode is %g. \n", cfe);
    
    mtrunc = cfe*sfe*mGMC;
    // printf("Truncation mass: %e Msun. \n", mtrunc);
    params[0] = mtrunc;

    norm_ctn = 1/trapz_log(mmin, mmax, params, (*fun_ptr1), 1000);
    f = norm_ctn * trapz_log(lowlimit, fmin(mstar,mmax), params, (*fun_ptr1), 1000);
    params[1] = norm_ctn;
    printf("Expected fraction %f \n", f);
    meanmc = trapz_log(mmin, mmax, params, (*fun_ptr2), 1000);

    // printf("Mean cluster mass= %e MSun. \n", meanmc);

    Nclus = cfe*mstar/meanmc;
    // printf("Expected number of clusters=%f \n", Nclus);
    for (int j=0; j<99;j++){
        totclusmass = 0.;
        N = 0;
        Nreal = gsl_ran_poisson(random_generator, Nclus);
        // printf("Actual number of cluster masses to be sampled=%d \n", Nreal);
        
        for(int i=1; i<=Nreal; i++){
            mclus = draw_mass(mmin, mmax, params, (*fun_ptr1), random_generator);
            if (mclus > lowlimit && totclusmass + mclus < mstar){
                // printf("Cluster %d has mass %e Msun.\n", i, mclus);
                totclusmass += mclus;
                N += 1;
            } else {
                continue;
                // if (mclus<mstar) printf("Cluster %d discarded, too low. \n", i);
                // if (mclus>mstar) printf("Cluster %d discarded, too high. \n", i);
            }
            // if (totclusmass > mstar){
            //     printf("Total star mass reached. No more masses to be sampled.\n");
            //     break;
            // }
        }

        // printf("TEST: Nclus*mc = %e. \n", Nreal*meanmc);
        // printf("TEST: CFE*Mstar = %e. \n", cfe*mstar);

        // totclusmass *= SOLAR_MASS;
        // printf("Total mass in clusters %e Msun in %d clusters from star particle of mass %e Msun.\n", totclusmass, N, mstar);
        // printf("Expected fraction of sampled masses %f \n", f);
        // printf("Obtained fraction of sampled masses %f \n", 1.0*N/Nreal);
        // printf("%d, %d, %f \n", Nreal, N, 1.0*N/Nreal);
    }

    return 0;
}