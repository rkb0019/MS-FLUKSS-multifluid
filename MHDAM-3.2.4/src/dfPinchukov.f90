! pure subroutine dfPinchukov( NEND, X, Y, B, C, D, iBCBGN, iBCEND, IERR )
subroutine dfPinchukov( NEND, X, Y, B, C, D, iBCBGN, iBCEND, IERR )

!                                                          == dfPinchukov.f90 ==
!//////////////////////////////////////////////////////////////////////////////!
!==                                                                          ==!
!==  PURPOSE:                                                                ==!
!==    dfPinchukov - interpolatory approximation by monotonised cubic        ==!
!==  splines with arbitrary end conditions                                   ==!
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
!==    iBCBGN     [INPUT]                                                    ==!
!==               flag of given left end condition:                          ==!
!==               = 0 any end condition is not given. In this case condition ==!
!==                 of continuity of third derivative at X(2) will be used;  ==!
!==               =+/-i value of i derivative must be used (0<i<4)           ==!
!==               If iBCBGN > 0 than value of derivative must be given in    ==!
!==               B(1). If iBCBGN < 0 than value of derivative will be       ==!
!==               calculated from devided difference;                        ==!
!==    iBCEND     [INPUT]                                                    ==!
!==               flag of given right end condition;                         ==!
!==                                                                          ==!
!==    B,C,D      [OUTPUT]                                                   ==!
!==               spline coefficients. The value of the spline approximation ==!
!==               at T is                                                    ==!
!==                        S(T) = ((D(I)*DX + C(I))*DX + B(I))*DX + Y(I)     ==!
!==               where X(I) <= T < X(I+1) and DX = T - X(I);                ==!
!==    IERR       [OUTPUT]                                                   ==!
!==               return code:                                               ==!
!==               =   0  normal run;                                         ==!
!==               =  -1  if NEND < 2;                                        ==!
!==               =  -2  if input abscissae are not ordered.                 ==!
!==  EXTERNALS:                                                              ==!
!==                                                                          ==!
!==  DIAGNOSTICS:                                                            ==!
!==                                                                          ==!
!==  NOTES:                                                                  ==!
!==                                                                          ==!
!//////////////////////////////////////////////////////////////////////////////!
!==                                                                          ==!
!==  CREATION DATE:  05.01.2002                                              ==!
!==                                                                          ==!
!==  REVISIONS HISTORY:                                                      ==!
!==                                                                          ==!
!//////////////////////////////////////////////////////////////////////////////!
  implicit  none

  integer,        intent(IN )                   :: NEND, iBCBGN, iBCEND
  integer,        intent(OUT)                   :: IERR
  real(kind = 8), intent(IN ), dimension (NEND) :: X, Y
  real(kind = 8), intent(OUT), dimension (NEND) :: B, C, D

  real(kind = 8), parameter :: Gamma  = 1.41421356237309504880168872420969808D0

!  integer,automatic        :: i, im, ip, NM1, NM2, NM3, bgnBC, endBC
!  real(kind = 8),automatic :: an, bn, rn, h, h1, h2, DN, DNM1, DNM2, DD, DI, DP, M, P

  integer        :: i, im, ip, NM1, NM2, NM3, bgnBC, endBC
  real(kind = 8) :: an, bn, rn, h, h1, h2, DN, DNM1, DNM2, DD, DI, DP, M, P



!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

  IERR   = 0

  if( NEND < 2 ) then
!                                                 Terminal - too few data points
    IERR   = -1
    return
  endif

  if( NEND == 2 ) then
    if( X(1) >= X(2) ) then
      IERR   = -2
      return
    endif
!                             Special case NEND == 2 => use linear interpolation
    B(1)   = (Y(2) - Y(1))/(X(2) - X(1))
    C(1)   = 0.0D0
    D(1)   = 0.0D0
    return
  endif

  do   I = 2,NEND
    IM     = I - 1
    H      = X(I) - X(IM)
    if( H <= 0.0D0 ) then
!                                           Terminal - x not monotone increasing
      IERR   = -2
      return
    endif

    C(I)   = H
    D(I)   = (Y(I) - Y(IM))/H
  enddo

  NM1    = NEND - 1
  NM2    = NEND - 2
  NM3    = NEND - 3

  DN     = D(NEND)
  DNM1   = D(NM1)
  DNM2   = D(NM2)

  bgnBC  = iBCBGN
  if( (NEND == 3) .and. (iBCBGN == -3) ) then
    bgnBC  = 3
    B(1)   = 0.0D0
  endif

  select case( abs( bgnBC ) )
    case ( 0 )
      D(1)   = C(3)
      C(1)   = C(2) + C(3)
      B(1)   = (D(2)*C(3)*(2.0D0*C(3) + 3.0D0*C(2)) + C(2)*C(2)*D(3))/C(1)
    case ( 1 )
      if( bgnBC < 0 ) then
        B(1)   = D(2)
      endif
      D(1)   = 1.0D0
      C(1)   = 0.0D0
    case ( 2 )
      if( bgnBC < 0 ) then
        B(1)   = 2.0D0*(D(3) - D(2))/(X(3) - X(1))
      endif
      D(1)   = 2.0D0
      C(1)   = 1.0D0
      B(1)   = 3.0D0*D(2) - 0.5D0*C(2)*B(1)
    case ( 3 )
      if( bgnBC < 0 ) then
        h1     = (D(3) - D(2))/(X(3) - X(1))
        h2     = (D(4) - D(3))/(X(4) - X(2))
        B(1)   = 6.0D0*(h2 - h1)/(X(4) - X(1))
      endif
      D(1)   = 1.0D0
      C(1)   = 1.0D0
      B(1)   = 2.0D0*D(2) + C(2)*C(2)*B(1)/6.0D0
  endselect

  do   i = 2, NM1 
    im     = i - 1
    ip     = i + 1
    h1     = C(ip)/D(im)

    DI     = D(i)
    if( DI < 0.0D0 ) DI    =-DI
    DP     = D(ip) 
    if( DP < 0.0D0 ) DP    =-DP

    DD     = Gamma*min( DI, DP )*(C(i) + C(ip))
    h      = C(i)*DI + C(ip)*DP

    p      = 1.0D0
    if( (h > 0.0D0) .and. (DD > 0.0D0) ) then
      p      = min( 1.0D0, DD/h )
    endif

    M      = max( -DD, min( C(ip)*D(i) + C(i)*D(ip), DD ) )
    B(i)   =  3.0D0           * M             - h1*B(im)
    D(i)   = (3.0D0/p - 1.0D0)*(C(ip) + C(i)) - h1*C(im)
  enddo

  endBC  = iBCEND
  if( NEND == 3 ) then
    if( (iBCEND == -3) .or. ((iBCEND == 0) .and. (iBCBGN == 0)) ) then
      endBC    = 3
      B(NEND)  = 0.0D0
    endif
  endif

  select case( abs( endBC ) )
    case ( 0 )
      an     = C(NEND) + C(NM1)
      bn     = C(NM1)
      rn     = (DN*C(NM1)*(2.0D0*C(NM1) + 3.0D0*C(NEND)) + C(NEND)*C(NEND)*DNM1)/an
    case ( 1 )
      if( endBC < 0 ) then
        B(NEND)  = DN
      endif
      an     = 0.0D0
      bn     = 1.0D0
      rn     = B(NEND)
    case ( 2 )
      if( endBC < 0 ) then
        B(NEND)  = 2.0D0*(DN - DNM1)/(X(NEND) - X(NM2))
      endif
      an     = 1.0D0
      bn     = 2.0D0
      rn     = 3.0D0*DN + 0.5D0*C(NEND)*B(NEND)
    case ( 3 )
      if( endBC < 0 ) then
        h1       = (DN   - DNM1)/(X(NEND) - X(NM2))
        h2       = (DNM1 - DNM2)/(X(NM1)  - X(NM3))
        B(NEND)  = 6.0D0*(h1 - h2)/(X(NEND) - X(NM3))
      endif
      an     = 1.0D0
      bn     = 1.0D0
      rn     = 2.0D0*DN + C(NEND)*C(NEND)*B(NEND)/6.0D0
  endselect

  B(NEND)  = (rn*D(NM1) - an*B(NM1))/(bn*D(NM1) - an*C(NM1))

  do  i = NM1, 1, -1
    B(i)   = (B(i) - C(i)*B(i+1))/D(i)
  enddo

  do   I = 2, NEND
    IM     = I - 1
    H      = C(I)
    H1     = (Y(I) - Y(IM))/H
    H2     = B(I) + B(IM) - H1 - H1
    D(IM)  = H2/(H*H)
    C(IM)  = (H1 - B(IM) - H2)/H
  enddo

end subroutine dfPinchukov
