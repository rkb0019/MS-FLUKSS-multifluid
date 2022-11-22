!pure subroutine dfLinearApprox( NEND, X, Y, M, U, S, IERR )
subroutine dfLinearApprox( NEND, X, Y, M, U, S, IERR )

!                                                       == dfLinearApprox.f90 ==
!//////////////////////////////////////////////////////////////////////////////!
!==                                                                          ==!
!==  PURPOSE:                                                                ==!
!==    dfLinearApprox computes values of discrete function in intermedian    ==!
!==  points with linear approximation.                                       ==!
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
!==    M          [INPUT]                                                    ==!
!==               number of points at which the values are to be evaluated;  ==!
!==    U          [INPUT]                                                    ==!
!==               vector containing the abscissae of the  points             ==!
!==               at which the values are to be evaluated.                   ==!
!==                                                                          ==!
!==    S          [OUTPUT]                                                   ==!
!==               vector of the values of the linear approximation at U;     ==!
!==    IERR       [OUTPUT]                                                   ==!
!==               return code:                                               ==!
!==               =   0  normal run;                                         ==!
!==               = -1 if U(I) < X(1)                    (Warning)           ==!
!==               = -2 if U(I) > X(NEND)                 (Warning)           ==!
!==               = -3 if U(I) < X(1) and U(J) > X(NEND) (Warning)           ==!
!==               = -4 if the abscissae are't ordered.                       ==!
!==  EXTERNALS:                                                              ==!
!==                                                                          ==!
!==  DIAGNOSTICS:                                                            ==!
!==                                                                          ==!
!==  NOTES:                                                                  ==!
!==                                                                          ==!
!//////////////////////////////////////////////////////////////////////////////!
!==                                                                          ==!
!==  CREATION DATE:  09.06.2002                                              ==!
!==                                                                          ==!
!==  REVISIONS HISTORY:                                                      ==!
!==                                                                          ==!
!//////////////////////////////////////////////////////////////////////////////!
  implicit  none

  integer,        intent(IN )                   :: NEND, M
  integer,        intent(OUT)                   :: IERR
  real(kind = 8), intent(IN ), dimension (NEND) :: X, Y
  real(kind = 8), intent(IN ), dimension (M   ) :: U
  real(kind = 8), intent(OUT), dimension (M   ) :: S

  interface dfGetInterval
    pure integer function dfGetInterval( NEND, X, XX )
      implicit  none
      integer,        intent(IN)                   :: NEND
      real(kind = 8), intent(IN), dimension (NEND) :: X
      real(kind = 8), intent(IN)                   :: XX
    end function dfGetInterval
  end interface dfGetInterval

!  integer,automatic :: I, JERR, K, KERR
  integer :: I, JERR, K, KERR


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

  if( X(1) <= X(NEND) ) then
    do   I = 2, NEND
      if( X(I-1) > X(I) ) then
        IERR   =-4
        return
      endif
    enddo
  else
    do   I = 2, NEND
      if( X(I-1) < X(I) ) then
        IERR   =-4
        return
      endif
    enddo
  endif

  JERR   = 0
  KERR   = 0
  do   K = 1, M
    I      = dfGetInterval( NEND, X, U(K) )
    if( I == 0 ) then
      S(K)   = Y(1)
      JERR   = 1
    else if( I == NEND ) then
      S(K)   = Y(NEND)
      KERR   = 2
    else
      S(K)   = Y(I) + (Y(I+1) - Y(I))/(X(I+1) - X(I))*(U(K) - X(I))
    endif
  enddo
  IERR   =-(JERR + KERR)

end subroutine dfLinearApprox
