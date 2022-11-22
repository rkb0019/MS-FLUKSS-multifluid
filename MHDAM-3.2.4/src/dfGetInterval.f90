! pure integer function dfGetInterval( NEND, X, XX )
integer function dfGetInterval( NEND, X, XX )

!                                                        == dfGetInterval.f90 ==
!//////////////////////////////////////////////////////////////////////////////!
!==                                                                          ==!
!==  PURPOSE:                                                                ==!
!==    dfGetInterval finds the largest integer I in 1 <= I <= NEND such      ==!
!==  that  X(I) <= XX (or X(I) >= XX). Precisely,                            ==!
!==                                                                          ==!
!==                        XX <= X(1)              0                         ==!
!==         if  X(K)    <= XX <= X(K+1)  then  I = K                         ==!
!==             X(NEND) <= XX                      NEND.                     ==!
!==                                                                          ==!
!==    In case of wrong usage the function returns error code:               ==!
!==  dfGetInterval = -1  if NEND < 1;                                        ==!
!==                = -2  if X is't ordered.                                  ==!
!==                                                                          ==!
!==  REFERENCES:                                                             ==!
!==                                                                          ==!
!==  ARGUMENTS:                                                              ==!
!==    NEND       [INPUT]                                                    ==!
!==               number of the data points;                                 ==!
!==    X          [INPUT]                                                    ==!
!==               vector containing the data points.                         ==!
!==               X must be monotone ordered so that X(I) <= X(I+1) for all  ==!
!==               I or X(I) >= X(I+1) for all I from [1,NEND];               ==!
!==    XX         [INPUT]                                                    ==!
!==               argument.                                                  ==!
!==  EXTERNALS:                                                              ==!
!==                                                                          ==!
!==  DIAGNOSTICS:                                                            ==!
!==                                                                          ==!
!==  NOTES:                                                                  ==!
!==                                                                          ==!
!//////////////////////////////////////////////////////////////////////////////!
!==                                                                          ==!
!==  CREATION DATE:  08.06.2002                                              ==!
!==                                                                          ==!
!==  REVISIONS HISTORY:                                                      ==!
!==                                                                          ==!
!//////////////////////////////////////////////////////////////////////////////!
  implicit  none

  integer,        intent(IN)                   :: NEND
  real(kind = 8), intent(IN), dimension (NEND) :: X
  real(kind = 8), intent(IN)                   :: XX

!  integer,automatic :: iLow, iHigh, iMiddle
  integer :: iLow, iHigh, iMiddle


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

  if( NEND < 1 ) then
    dfGetInterval  = -1
    return
  endif

  if( X(1) < X(NEND) ) then
!                                                     Monotonly increased points
    if( XX < X(1) ) then
      dfGetInterval  = 0
    else if( XX > X(NEND) ) then
      dfGetInterval  = NEND
    endif

    iLow   = 1
    iHigh  = NEND
10  continue
      iMiddle  = (iLow + iHigh)/2
      if( iMiddle == iLow ) then
        dfGetInterval  = iLow
        return
      endif
      if( XX < X(iMiddle) ) then
        iHigh  = iMiddle
      else
        iLow   = iMiddle
      endif
    goto   10
  else if( X(1) > X(NEND) ) then
!                                                     Monotonly decreased points
    if( XX > X(1) ) then
      dfGetInterval  = 0
    else if( XX < X(NEND) ) then
      dfGetInterval  = NEND
    endif

    iLow   = 1
    iHigh  = NEND
20  continue
      iMiddle  = (iLow + iHigh)/2
      if( iMiddle == iLow ) then
        dfGetInterval  = iLow
        return
      endif
      if( XX > X(iMiddle) ) then
        iHigh  = iMiddle
      else
        iLow   = iMiddle
      endif
    goto   20
  else
    if( XX < X(1) ) then
      dfGetInterval  = 0
    else
      dfGetInterval  = NEND
    endif
  endif

end function dfGetInterval
