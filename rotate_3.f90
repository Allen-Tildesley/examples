! rotate_3.f90
SUBROUTINE flip3 ( exiold, eyiold, eziold, dotmin, exinew, eyinew, ezinew )

c    *******************************************************************
c    ** chooses a random displacement on the surface of a unit sphere.**
c    **                                                               **
c    ** reference:                                                    **
c    **                                                               **
c    ** marsaglia, ann maths stat 43, 645, 1972.                      **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** REAL      dot                   dot product of old and NEW    **
c    **                                 axial vectors                 **
c    ** REAL      dotmin                PARAMETER to adjust maximum   **
c    **                                 displacement. dotmin should   **
c    **                                 be less than one              **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** flip3 uses a rejection technique to create a trial            **
c    ** orientation of molecule i subject to the constraint that      **
c    ** the cosine of the angle between the old and NEW axial         **
c    ** vectors is greater than ( 1.0 - dotmin ).                     **
c    *******************************************************************

        REAL    exiold, eyiold, eziold, exinew, eyinew, ezinew, dotmin

        REAL    dot, xi1, xi2, xi, xisq
        REAL    ranf, dummy

c    *******************************************************************

c    ** initialise dot **

        dot  = 0.0

c    ** iterative loop **

1000    IF ( ( 1.0 - dot ) .GE. dotmin ) THEN

c       ** initialise xisq **

           xisq = 1.0

c       ** inner iterative loop **

2000       IF ( xisq .GE. 1.0 ) THEN

              xi1  = ranf ( dummy ) * 2.0 - 1.0
              xi2  = ranf ( dummy ) * 2.0 - 1.0
              xisq = xi1 * xi1 + xi2 * xi2

              GOTO 2000

           ENDIF

           xi = SQRT ( 1.0 - xisq )
           exinew = 2.0 * xi1 * xi
           eyinew = 2.0 * xi2 * xi
           ezinew = 1.0 - 2.0 * xisq
           dot    = exinew * exiold + eyinew * eyiold + ezinew * eziold

           GOTO 1000

        ENDIF

        RETURN
        END



