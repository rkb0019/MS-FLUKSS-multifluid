! pure subroutine dfSplineEval( NEND, X, Y, B, C, D, M, U, S, DS, IERR )
subroutine dfSplineEval( NEND, X, Y, B, C, D, M, U, S, DS, IERR )

!                                                         == dfSplineEval.f90 ==
!//////////////////////////////////////////////////////////////////////////////!
!==                                                                          ==!
!==  PURPOSE:                                                                ==!
!==    dfSplineEval evaluates of a cubic spline                              ==!
!==                                                                          ==!
!==  REFERENCES:                                                             ==!
!==                                                                          ==!
!==  ARGUMENTS:                                                              ==!
!==    NEND       [INPUT]                                                    ==!
!==               number of the data points;                                 ==!
!==    X          [INPUT]                                                    ==!
!==               vector containing the abscissae of the data points.        ==!
!==               X must be ordered so that X(I) < X(I+1);                   ==!
!==    Y          [INPUT]                                                    ==!
!==               vector containing the ordinates of the data points;        ==!
!==    B,C,D      [INPUT]                                                    ==!
!==               spline coefficients. The value of the spline               ==!
!==               approximation at T is                                      ==!
!==                 S(T) = ((D(I)*DX + C(I))*DX + B(I))*DX + Y(I)            ==!
!==               where X(I) <= T < X(I+1) and DX = T - X(I);                ==!
!==    M          [INPUT]                                                    ==!
!==               number of points at which the cubic spline is to be        ==!
!==               evaluated;                                                 ==!
!==    U          [INPUT]                                                    ==!
!==               vector containing the abscissae of the  points             ==!
!==               at which the cubic spline is to be evaluated.              ==!
!==                                                                          ==!
!==    S          [OUTPUT]                                                   ==!
!==               vector of the values of the spline approximation at U;     ==!
!==    DS         [OUTPUT]                                                   ==!
!==               vector of the values of approximation of first             ==!
!==               derivative of the spline at U;                             ==!
!==    IERR       [OUTPUT]                                                   ==!
!==               return code:                                               ==!
!==               =   0  normal run;                                         ==!
!==               = -1 if U(I) < X(1)                    (Warning)           ==!
!==               = -2 if U(I) > X(NEND)                 (Warning)           ==!
!==               = -3 if U(I) < X(1) and U(J) > X(NEND) (Warning)           ==!
!==  EXTERNALS:                                                              ==!
!==                                                                          ==!
!==  DIAGNOSTICS:                                                            ==!
!==                                                                          ==!
!==  NOTES:                                                                  ==!
!==                                                                          ==!
!//////////////////////////////////////////////////////////////////////////////!
!==                                                                          ==!
!==  CREATION DATE:  09.02.1999                                              ==!
!==                                                                          ==!
!==  REVISIONS HISTORY:                                                      ==!
!==                                                                          ==!
!//////////////////////////////////////////////////////////////////////////////!
  implicit  none

  integer,        intent(IN )                   :: NEND, M
  integer,        intent(OUT)                   :: IERR
  real(kind = 8), intent(IN ), dimension (NEND) :: X, Y, B, C, D
  real(kind = 8), intent(IN ), dimension (M   ) :: U
  real(kind = 8), intent(OUT), dimension (M   ) :: S, DS

!  integer,automatic        :: I, JERR, K, KERR
!  real(kind = 8),automatic :: DD, DX, HB, HD, HDC

  integer        :: I, JERR, K, KERR
  real(kind = 8) :: DD, DX, HB, HD, HDC



!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

  JERR   = 0
  KERR   = 0
  I      = 1
  do   K = 1,M
10  continue
    DX     = U(K) - X(I)
    if( DX == 0.0D0 ) then
      S(K)   = Y(I)
      DS(K)  = B(I)
      cycle
    endif

    if( DX < 0.0D0 ) then
      if( I == 1 ) THEN
        JERR   = 1
      else
        I      = I + 1
        goto   10
      endif
    else
20    continue
      if( I >= NEND ) then
        if( U(K) - X(NEND) > epsilon( DX ) ) KERR   = 2
        DX     = U(K) - X(NEND-1)
        I      = NEND - 1
      else
        DD     = U(K) - X(I+1)
        if( DD >= 0.0D0 ) then
          DX     = DD
          I      = I + 1
          goto   20
        endif
        if( DX == 0.0D0 ) then
          S(K)   = Y(I)
          DS(K)  = B(I)
          cycle
        endif
      endif
    endif

    HD     = D(I)*DX
    HDC    = HD + C(I)
    HB     = B(I)

    S(K)   = (HDC*DX + HB)*DX + Y(I)
    DS(K)  = (HD + HDC + HDC)*DX + HB
  enddo
  IERR   =-(JERR + KERR)

end subroutine dfSplineEval
