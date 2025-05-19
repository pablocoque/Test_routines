
! TEST PROGRAM - READ THE HEADER OF F_CFE.F90 FOR DETAILS

PROGRAM testCFElocal
    use CFElocalmod
    REAL pc,msun,myr,rhoMW,sigmaMW,csMW,cfe(1:4)
    
    pc=3.1e16 !parsec in meters
    Msun=2.e30 !solar mass in kg
    Myr=3.16e13 !million years in seconds
    rhoMW= 0.03*msun/pc**3. !volume density
    sigmaMW= 7.e3 !velocity dispersion
    csMW=0.2e3 !sound speed
    cfe=f_cfelocal(rhoMW,sigmaMW,csMW) !CFE
    PRINT*,cfe(1)*100.
END PROGRAM