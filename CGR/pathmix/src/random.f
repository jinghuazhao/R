      DOUBLE PRECISION FUNCTION RANU()
C---
C--- Two entry points are defined:
C---
C--- X = RANU() -- Return a uniformly distributed random number in the
C---               closed range 0.0 to 1.0
C---
C--- CALL RNDMIZ(SEED) -- set starting seed for RANU (details below)
C---
C---
C--- The following algorithm derives from the "Lehmer Congruence" method as
C--- outlined in K. D. Tocher's The Art of Simulation (1963, pp 75-81).  My
C--- implementation is summarized as follows:
C---
C--- A modulus "m" is chosen from the Mersenne number series compatible with
C--- the word length of the machine, such that m = 2^n - 1 and such that m
C--- is prime.  A "seed" is maintained as an n bit integer in the range 1 to
C--- m - 1.  The update consists of multiplying the seed by a prime raised
C--- to a power which is prime to m - 1, and then extracting modulus m of
C--- the result.  The updated seed is then normalized to 1.0 and returned
C--- as the pseudo-random result.
C---
C--- This method has the properties that a given seed in the sequence is
C--- non-repeatable for m - 1 cycles, and gives low serial correlations.
C---
C--- Note: I have made the following assumptions regarding number formats:
C---    unsigned long -- precise to at least 31 bits
C---    double        -- mantissa precise to at least 46 bits.
C---
C--- written February, 1986 by Skip Russell
C--- converted from C to Fortran June 1988 by SR
C---

      DOUBLE PRECISION  SEED
      real              tarray(2), result
      equivalence       ( tarray(1), iequiv )

C     define two constants:  2^31 - 1 and 7^5
      DOUBLE PRECISION  FMOD, FMUL
      PARAMETER   ( FMOD=2.D0**31-1, FMUL=7.D0**5 ) 

C     define default starting seed for random number generator
      DATA RSEED  / 12345.6789 /

C     get new seed
      RSEED = MOD( RSEED*FMUL, FMOD )

C     normalize to 0,1
      RANU =  RSEED / FMOD
      RETURN

      ENTRY RNDMIZ( SEED )
C---
C--- Set the starting seed for the ranu function (below).  If the seed
C--- is specified as 0, we choose a pseudo-random seed based on the
C--- machine time.
C---
      IF ( SEED .NE. 0 ) THEN
         RSEED = SEED
      ELSE

C if unspecified, use the system time as the seed

         call dtime( tarray, result )
         SEED = iequiv + 1.0

      END IF
      END
