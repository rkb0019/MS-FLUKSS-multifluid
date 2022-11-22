
! This subroutine generates random numbers that have standard 
! Gaussian distribution N(0,1). Prior using this subroutine 
! one must 
! initialize standard random number generator which produces
! uniformly distribited 
! numbers on [0,1]. 
! In Fortran it can be done by calling "random_seed"

      subroutine  GaussDistribution(rn)

	REAL*8 :: rn, d_PI
      REAL*8 :: harv1,harv2
  
      parameter (d_PI     = 3.14159265358979323846264338327950288D0 )

	call random_number(harv1)

	rn = 2D0 * harv1 - 1D0 
!	call random_number(harv2)
!      
!	do
!	  if (harv1 == 0.D0) then
!	    call random_number(harv1) 
!	    cycle
!	  else 
!	    rn = sqrt(-2.d0*log(harv1))*cos(2.0D0*d_PI*harv2)
!	    exit
!     endif
!	enddo

	return
      end