! rotate_1.f90
SUBROUTINE flip1 ( exiold, eyiold, eziold, dphimx, dcosmx, exinew, eyinew, ezinew  )

c    *******************************************************************
c    ** makes a random change in the polar angles phi and theta.      **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** REAL      dcosmx  maximum change in COS(theta)                **
c    ** REAL      dphimx  maximum change in phi                       **
c    ** REAL      phiold  phi in the old state                        **
c    ** REAL      phinew  phi in the NEW trial state                  **
c    ** REAL      cosold  COS(theta) in the old state                 **
c    ** REAL      cosnew  COS(theta) in the NEW trial state           **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** flip1 makes a random change in phi and COS(theta).            **
c    ** the maximum allowed changes in these variables are controlled **
c    ** by the parameters dphimx and dcosmx respectively.             **
c    ** phi and theta are the euler angles describing the orientation **
c    ** of the axial vector. this method can be readily extended to   **
c    ** polyatomics by changing the third euler angle , psi, in the   **
c    ** range zero to twopi. it would be faster to PASS the variables **
c    ** phiold and cosold directly to flip1 IF they are available in  **
c    ** the main PROGRAM. similarly phinew and cosnew could be        **
c    ** passed directly back through the SUBROUTINE header            **
c    *******************************************************************

        REAL        exiold, eyiold, eziold, exinew, eyinew, ezinew
        REAL        dphimx, dcosmx

        REAL        cosnew, cosold, phinew, phiold, sinnew
        REAL        twopi, pi
        REAL        ranf, dummy

        PARAMETER ( twopi = 6.2831853 )

c    *******************************************************************

c    ** convert the axial vector to the euler angles **

        cosold = eziold
        phiold = ATAN2 ( eyiold, exiold )

c    ** perform the displacements **

        phinew = phiold + ( 2.0 * ranf ( dummy ) - 1.0 ) * dphimx
        phinew = phinew - ANINT ( phinew / twopi ) * twopi
        cosnew = cosold + ( 2.0 * ranf ( dummy ) - 1.0 ) * dcosmx
        cosnew = cosnew - ANINT ( cosnew / 2.0 ) * 2.0
        sinnew = SQRT ( 1.0 - cosnew * cosnew )

c    ** convert the euler angles to axial vectors **

        exinew = COS ( phinew ) * sinnew
        eyinew = SIN ( phinew ) * sinnew
        ezinew = cosnew

        RETURN
        END
