c============================================================================
c
!        SUBROUTINE GAMMAFunc (X,GA)
!C
!C       ==================================================
!C       Purpose: Compute the gamma function Gamma(x)
!C       Input :  x  --- Argument of Gamma(x)
!C                       ( x is not equal to 0,-1,-2, ... )
!C       Output:  GA --- Gamma(x)
!C       ==================================================
!C
!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!        DIMENSION G(26)
!        PI=3.141592653589793D0
!        IF (X.EQ.INT(X)) THEN
!           IF (X.GT.0.0D0) THEN
!              GA=1.0D0
!              M1=X-1
!              DO 10 K=2,M1
!10               GA=GA*K
!           ELSE
!              GA=1.0D+300
!           ENDIF
!        ELSE
!           IF (DABS(X).GT.1.0D0) THEN
!              Z=DABS(X)
!              M=INT(Z)
!              R=1.0D0
!              DO 15 K=1,M
!15               R=R*(Z-K)
!              Z=Z-M
!           ELSE
!              Z=X
!           ENDIF
!           DATA G/1.0D0,0.5772156649015329D0,
!     &          -0.6558780715202538D0, -0.420026350340952D-1,
!     &          0.1665386113822915D0,-.421977345555443D-1,
!     &          -.96219715278770D-2, .72189432466630D-2,
!     &          -.11651675918591D-2, -.2152416741149D-3,
!     &          .1280502823882D-3, -.201348547807D-4,
!     &          -.12504934821D-5, .11330272320D-5,
!     &          -.2056338417D-6, .61160950D-8,
!     &          .50020075D-8, -.11812746D-8,
!     &          .1043427D-9, .77823D-11,
!     &          -.36968D-11, .51D-12,
!     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
!           GR=G(26)
!           DO 20 K=25,1,-1
!20            GR=GR*Z+G(K)
!           GA=1.0D0/(GR*Z)
!           IF (DABS(X).GT.1.0D0) THEN
!              GA=GA*R
!              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
!           ENDIF
!        ENDIF
!
!       return
!
!      end
c
c============================================================================
c
      SUBROUTINE GAMMAFunc (X,GA)
      use parameters
      use GAMMA_MOD
      implicit none
      real(pr) :: X,GA
      CALL TGAMMA(X, GA)
      return
      end
c
c============================================================================
c

