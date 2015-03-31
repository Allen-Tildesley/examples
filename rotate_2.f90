! rotate_2.f90
SUBROUTINE flip2 ( exiold, eyiold, eziold, dgamax, exinew, eyinew, ezinew         )

c    *******************************************************************
c    ** performs random rotation about space-fixed axes.              **
c    **                                                               **
c    ** reference:                                                    **
c    **                                                               **
c    ** barker and watts, chem phys letts 3, 144, 1969.               **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** REAL      dgamax  maximum angular displacement in radians     **
c    ** INTEGER   iaxis   space fixed axis for rotation               **
c    **                   (1 = x, 2 = y, 3 = z)                       **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** flip2 chooses one of the three space-fixed axes at random     **
c    ** and rotates the molecule around this axis by dgamma radians.  **
c    ** the maximum angular displacement is dgamax. this method can   **
c    ** readily extended to polyatomic molecules.                     **
c    *******************************************************************

        REAL    exiold, eyiold, eziold, exinew, eyinew, ezinew, dgamax

        REAL    cosdg, sindg, dgamma
        REAL    ranf, dummy
        INTEGER iaxis

c    *******************************************************************

c    ** choose a space fixed axis at random **

        iaxis = INT ( 3.0 * ranf ( dummy ) ) + 1

c    ** choose a random rotation **

        dgamma = ( 2.0 * ranf ( dummy ) - 1.0 ) * dgamax

c    ** set up the rotation matrix **

        cosdg = COS ( dgamma )
        sindg = SIN ( dgamma )

c    ** perform rotations **

        IF ( iaxis .EQ. 1 ) THEN

           exinew = exiold
           eyinew = cosdg * eyiold + sindg * eziold
           ezinew = cosdg * eziold - sindg * eyiold

        ELSE IF ( iaxis .EQ. 2 ) THEN

           exinew = cosdg * exiold - sindg * eziold
           eyinew = eyiold
           ezinew = cosdg * eziold + sindg * exiold

        ELSE

           exinew = cosdg * exiold + sindg * eyiold
           eyinew = cosdg * eyiold - sindg * exiold
           ezinew = eziold

        ENDIF

        RETURN
        END
