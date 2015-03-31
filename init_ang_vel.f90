! init_ang_vel.f90
SUBROUTINE init_ang_vel ( temp, inert )

c    *******************************************************************
c    ** centre of mass and angular velocities for linear molecules    **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER n                   the number of molecules           **
c    ** REAL    rx(n),ry(n),rz(n)   positions                         **
c    ** REAL    vx(n),vy(n),vz(n)   velocities                        **
c    ** REAL    ex(n),ey(n),ez(n)   orientations                      **
c    ** REAL    ox(n),oy(n),oz(n)   space-fixed angular velocities    **
c    ** REAL    temp                reduced temperature               **
c    ** REAL    inert               reduced moment of inertia         **
c    **                                                               **
c    ** supplied routines:                                            **
c    **                                                               **
c    ** SUBROUTINE comvel ( temp )                                    **
c    **    sets the centre of mass velocities for a configuration of  **
c    **    linear molecules at a given temperature.                   **
c    ** SUBROUTINE angvel ( temp, inert )                             **
c    **    sets the angular velocities for a configuration of linear  **
c    **    molecules at a given temperature.                          **
c    ** REAL FUNCTION ranf ( dummy )                                  **
c    **    returns a uniform random variate on the range zero to one  **
c    ** REAL FUNCTION gauss ( dummy )                                 **
c    **    returns a uniform random normal variate from a             **
c    **    distribution WITH zero mean and unit variance.             **
c    **                                                               **
c    ** units:                                                        **
c    **                                                               **
c    ** we assume unit molecular mass and employ lennard-jones units  **
c    **       property                      units                     **
c    **       rx, ry, rz           (epsilon/m)**(1.0/2.0)             **
c    **       ox, oy, oz           (epsilon/m*sigma**2)**(1.0/2.0)    **
c    **       inert                 m*sigma**2                        **
c    *******************************************************************

        SUBROUTINE comvel ( temp )

        COMMON / block1 / rx, ry, rz, ex, ey, ez,
     :                    vx, vy, vz, ox, oy, oz

c    *******************************************************************
c    ** translational velocities from maxwell-boltzmann distribution  **
c    **                                                               **
c    ** the distribution is determined by temperature and (unit) mass.**
c    ** this routine is general, and can be used for atoms, linear    **
c    ** molecules, and non-linear molecules.                          **
c    **                                                               **
c    ** routine referenced:                                           **
c    **                                                               **
c    ** REAL FUNCTION gauss ( dummy )                                 **
c    **    returns a uniform random normal variate from a             **
c    **    distribution WITH zero mean and unit variance.             **
c    *******************************************************************

        INTEGER     n
        PARAMETER ( n = 108 )

        REAL        rx(n), ry(n), rz(n), ex(n), ey(n), ez(n)
        REAL        vx(n), vy(n), vz(n), ox(n), oy(n), oz(n)
        REAL        temp

        REAL        rtemp, sumx, sumy, sumz
        REAL        gauss, dummy
        INTEGER     i

c    *******************************************************************

        rtemp = SQRT ( temp )

        DO 100 i = 1, n

           vx(i) = rtemp * gauss ( dummy )
           vy(i) = rtemp * gauss ( dummy )
           vz(i) = rtemp * gauss ( dummy )

100     CONTINUE

c    ** remove net momentum **

        sumx = 0.0
        sumy = 0.0
        sumz = 0.0

        DO 200 i = 1, n

           sumx = sumx + vx(i)
           sumy = sumy + vy(i)
           sumz = sumz + vz(i)

200     CONTINUE

        sumx = sumx / REAL ( n )
        sumy = sumy / REAL ( n )
        sumz = sumz / REAL ( n )

        DO 300 i = 1, n

           vx(i) = vx(i) - sumx
           vy(i) = vy(i) - sumy
           vz(i) = vz(i) - sumz

300     CONTINUE

        RETURN
        END



        SUBROUTINE angvel ( temp, inert )

        COMMON / block1 / rx, ry, rz, ex, ey, ez,
     :                    vx, vy, vz, ox, oy, oz

c    *******************************************************************
c    ** angular velocities from the maxwell-boltzmann distribution.   **
c    **                                                               **
c    ** the distribution is determined by temperature and inertia.    **
c    ** this routine is specific to linear molecules.                 **
c    ** it chooses the direction of the angular velocity randomly but **
c    ** perpendicular to the molecular axis. the square of the        **
c    ** magnitude of the angular velocity is chosen from an           **
c    ** exponential distribution. there is no attempt to set the      **
c    ** total angular momentum to zero.                               **
c    **                                                               **
c    ** routine referenced:                                           **
c    **                                                               **
c    ** REAL FUNCTION ranf ( dummy )                                  **
c    **    returns a uniform random variate on the range zero to one  **
c    *******************************************************************

        INTEGER     n
        PARAMETER ( n = 108 )

        REAL        rx(n), ry(n), rz(n), ex(n), ey(n), ez(n)
        REAL        vx(n), vy(n), vz(n), ox(n), oy(n), oz(n)
        REAL        temp, inert

        REAL        norm, dot, osq, o, mean
        REAL        xisq, xi1, xi2, xi
        REAL        ranf, dummy
        INTEGER     i

c       ****************************************************************

        mean = 2.0 * temp / inert

c    ** set direction of the angular velocity **

        DO 100 i = 1, n

c       ** choose a random vector in space **

           xisq = 1.0

1000       IF ( xisq .GE. 1.0 ) THEN

              xi1  = ranf ( dummy ) * 2.0 - 1.0
              xi2  = ranf ( dummy ) * 2.0 - 1.0
              xisq = xi1 * xi1 + xi2 * xi2

              go to 1000

           ENDIF

           xi    = SQRT ( 1.0 - xisq )
           ox(i) = 2.0 * xi1 * xi
           oy(i) = 2.0 * xi2 * xi
           oz(i) = 1.0 - 2.0 * xisq

c       ** constrain the vector to be perpendicular to the molecule **

           dot   = ox(i) * ex(i) + oy(i) * ey(i) + oz(i) * ez(i)
           ox(i) = ox(i) - dot * ex(i)
           oy(i) = oy(i) - dot * ey(i)
           oz(i) = oz(i) - dot * ez(i)

c       ** renormalize **

           osq   = ox(i) * ox(i) + oy(i) * oy(i) + oz(i) * oz(i)
           norm  = SQRT ( osq )
           ox(i) = ox(i) / norm
           oy(i) = oy(i) / norm
           oz(i) = oz(i) / norm

c       ** choose the magnitude of the angular velocity **

           osq   = - mean * LOG ( ranf ( dummy ) )
           o     = SQRT ( osq )
           ox(i) = o * ox(i)
           oy(i) = o * oy(i)
           oz(i) = o * oz(i)

100     CONTINUE

        RETURN
        END



        REAL FUNCTION gauss ( dummy )

c    *******************************************************************
c    ** random variate from the standard normal distribution.         **
c    **                                                               **
c    ** the distribution is gaussian WITH zero mean and unit variance.**
c    **                                                               **
c    ** reference:                                                    **
c    **                                                               **
c    ** knuth d, the art of computer programming, (2nd edition        **
c    **    addison-wesley), 1978                                      **
c    **                                                               **
c    ** routine referenced:                                           **
c    **                                                               **
c    ** REAL FUNCTION ranf ( dummy )                                  **
c    **    returns a uniform random variate on the range zero to one  **
c    *******************************************************************

        REAL        a1, a3, a5, a7, a9
        PARAMETER ( a1 = 3.949846138, a3 = 0.252408784 )
        PARAMETER ( a5 = 0.076542912, a7 = 0.008355968 )
        PARAMETER ( a9 = 0.029899776                   )

        REAL        sum, r, r2
        REAL        ranf, dummy
        INTEGER     i

c    *******************************************************************

        sum = 0.0

        DO 10 i = 1, 12

           sum = sum + ranf ( dummy )

10      CONTINUE

        r  = ( sum - 6.0 ) / 4.0
        r2 = r * r

        gauss = (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 +a1 )
     :          * r

        RETURN
        END



        REAL FUNCTION ranf ( dummy )

c    *******************************************************************
c    ** returns a uniform random variate in the range 0 to 1.         **
c    **                                                               **
c    **                 ***************                               **
c    **                 **  warning  **                               **
c    **                 ***************                               **
c    **                                                               **
c    ** good random number generators are machine specific.           **
c    ** please USE the one recommended for your machine.              **
c    *******************************************************************

        INTEGER     l, c, m
        PARAMETER ( l = 1029, c = 221591, m = 1048576 )

        INTEGER     seed
        REAL        dummy
        SAVE        seed
        DATA        seed / 0 /

c    *******************************************************************

        seed = MOD ( seed * l + c, m )
        ranf = REAL ( seed ) / m

        RETURN
        END



