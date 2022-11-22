
 subroutine fft2dmain(nx,ny,dx,phi)
  
! 2d fft subroutine
! isign = 1 for forward transform
! isign =-1 for inverse transform

       implicit             none
	  
       integer nx,ny
	   real(kind=8) :: dx

       integer           :: ix, iy, idum1=-1298, idum2=298, idum3=836       
       real(kind=8)      :: x,y,pi=3.1459, lx=10., ly=10.
       real(kind = 8)    :: phi(nx,ny)
       real(kind = 8)    :: ran1, ran2, ran3
       real(kind = 8)    :: fas1, fas2, fas3
       real(kind = 8)    :: kx_max, ky_max, k_max, k, phi_aver
       real(kind = 8)    :: amp=1.e3

       
      
       kx_max = 2.*pi/lx * dfloat(nx-1)/2.
       ky_max = 2.*pi/ly * dfloat(ny-1)/2.
       k_max  = 2./3. * sqrt(kx_max**2 + ky_max**2)



        do iy = 1, ny
           do ix = 1, nx

              kx_max       = 2.*pi/lx * dfloat(ix-1)
              ky_max       = 2.*pi/ly * dfloat(iy-1)
              k            = sqrt(kx_max**2 + ky_max**2)

              
!                 if (( k < k_max/1.5).and.(k>k_max/4)) then
                 if (( k < k_max/1.5).and.(k>k_max/2)) then
                 fas2         = 2.*pi*dble(ran2(idum2))
                 fas3         = 2.*pi*dble(ran3(idum3))

				 call random_number(fas2)
				 fas2         = 2.*pi*fas2

                 
!                 phi(ix,iy)    = amp*cos(fas2)/(1 + k**2)
				 phi(ix,iy)    = amp*sin(fas2)/(1 + k**2)
                 
                 else
                    phi(ix,iy)    = 0.

                 end if


              enddo
        enddo
    

!!Inverse Fourier Transform ####
        CALL    fft2d(-1, phi, nx, ny)


end subroutine fft2dmain


subroutine fft2dmain2()
  
! 2d fft subroutine
! isign = 1 for forward transform
! isign =-1 for inverse transform

       implicit             none
	         

       integer           :: ix, iy, nx, ny, idum1=-1298, idum2=298, idum3=836
       parameter            (nx=180, ny=180)
       real(kind=8)      :: x,y,pi=3.14, lx=10., ly=10.
       real(kind = 8)    :: phi(nx,ny)
       real(kind = 8)    :: ran1, ran2, ran3
       real(kind = 8)    :: fas1, fas2, fas3
       real(kind = 8)    :: kx_max, ky_max, k_max, k, phi_aver
       real(kind = 8)    :: amp=1.e3


       call random_seed

       open(88, file='output.dat', status='unknown')


       kx_max = 2.*pi/lx * dfloat(nx-1)/2.
       ky_max = 2.*pi/ly * dfloat(ny-1)/2.
       k_max  = 2./3. * sqrt(kx_max**2 + ky_max**2)



        do iy = 1, ny
           do ix = 1, nx

              kx_max       = 2.*pi/lx * dfloat(ix-1)
              ky_max       = 2.*pi/ly * dfloat(iy-1)
              k            = sqrt(kx_max**2 + ky_max**2)

              
                 if( k < k_max/1.5) then
                 fas2         = 2.*pi*dble(ran2(idum2))
                 fas3         = 2.*pi*dble(ran3(idum3))

				 call random_number(fas2)
				 fas2         = 2.*pi*fas2

                 
!                 phi(ix,iy)    = amp*cos(fas2)/(1 + k**2)
				 phi(ix,iy)    = amp*sin(fas2)/(1 + k**2)
                 
                 else
                    phi(ix,iy)    = 0.

                 end if


              enddo
        enddo
    
	    write(97,*) 'VARIABLES="X" "Y" "phi"'
        write(97,*) 'ZONE I= ',nx, 'J= ',ny, ' DATAPACKING=POINT'
        do iy = 1, ny
           do ix = 1, nx
                    write(97,*) ix,' ',iy,' ',phi(ix,iy) 
              enddo              
        enddo
        CALL spectra(nx, ny, k_max, phi)

!!Inverse Fourier Transform ####
        CALL    fft2d(-1, phi, nx, ny)


!! Phi is a real variable #####
        do iy = 1, ny
           do ix = 1, nx
                    write(99,*)phi(ix,iy) 
              enddo
              write(99,*)
        enddo

        phi_aver = 0.D0

		!! Phi is a real variable #####
		write(98,*) 'VARIABLES="X" "Y" "phi"'
        write(98,*) 'ZONE I= ',nx, 'J= ',ny, ' DATAPACKING=POINT'
        do iy = 1, ny
           do ix = 1, nx
                    write(98,*) ix,' ',iy,' ',phi(ix,iy) 
					phi_aver = phi_aver + phi(ix,iy) 
              enddo              
        enddo



end subroutine fft2dmain2






subroutine               spectra(nx, ny, kmax, ph) 
  implicit               none
  integer           ::   ix, iy, nx, ny
  real(kind = 8)    ::   kmax,dk, toten, km, kp, NN, t1, kk, k1, cut, pi=3.14, nrm
  real(kind = 8)    ::   kx_max, ky_max, k, lx=10., ly=10.
  REAL(kind = 8)    ::   ph(nx,ny)

 
  
! this file creats energy spectrum for density fluctuations (spc**)
 
      
      open (88,file='spc.dat',status='unknown')

     
      cut = 2./3. *2./3. * kmax
      kmax      = cut 
      kmax      = dsqrt(kmax)
      dk        = 0.2d0
      k1        = 0.0
      nrm = 1.d0 /(nx*ny)
 

      do  k1    = 2*dk, kmax, dk
     
      toten     = 0.0
     

      NN        = 0.d0
      km        = k1 - dk
      kp        = k1 + dk

      do ix = 1, nx
         do iy = 1, ny
            kx_max       = 2.*pi/lx * dfloat(ix-1)
            ky_max       = 2.*pi/ly * dfloat(iy-1)
            k            = sqrt(kx_max**2 + ky_max**2)
    
            if((k .ge. km) .and. (k .le. kp)) then
            toten   = toten   +  k*abs(ph(ix,iy))**2

            NN   = NN + 1.d0
         END IF
         
      enddo
   enddo

            if(NN .eq. 0.d0) goto 13
   
            toten    = .5   * nrm * toten
            toten    = 2.d0 * pi  * k1 * toten/NN

         
            write (88,'(2(1pe15.7))') k1, toten
            

13            CONTINUE
         end do

            close(88)

          
end subroutine spectra








subroutine     fft2d(isign, psi, nx, ny)
  integer      isign,nx,ny
  REAL*8       psi(nx,ny)
  REAL*8       psit(ny,nx)
  REAL(KIND=8), DIMENSION (:), ALLOCATABLE :: wsave



  ALLOCATE(wsave(2*nx+15))


   CALL drffti(nx,wsave)

   if (isign .eq. 1) then

      DO j = 1, ny
         CALL drfftf(nx,psi(:,j),wsave)
      END DO

      CALL   mattrans(nx, psi)

      DO j = 1, ny
         CALL drfftf(nx,psi(:,j),wsave)
      END DO

      CALL   mattrans(nx, psi)

   end if


!    if (isign .eq. -1) then
!
!      DO j = 1, ny
!         CALL drfftb(nx,psi(:,j),wsave)
!      END DO
!
!      CALL   mattrans(nx, psi)
!
!      DO j = 1, ny
!         CALL drfftb(nx,psi(:,j),wsave)
!      END DO
!
!      CALL   mattrans(nx, psi)
!
!       psi = psi / dfloat(nx*ny)
!
!    end if

    if (isign .eq. -1) then

      DO j = 1, ny
         CALL drfftb(nx,psi(:,j),wsave)
      END DO

      CALL   mattrans_g(nx, ny, psi, psit)

      DO j = 1, nx
         CALL drfftb(ny,psit(:,j),wsave)
      END DO

      CALL   mattrans_g(ny, nx, psit, psi)

      psi = psi / dfloat(nx*ny)

    end if

    DEALLOCATE(wsave)


end subroutine fft2d


subroutine             mattrans(ndim, matrix)
  implicit             none
  integer           :: ix, iy, ndim
  REAL(kind = 8) :: matrix(ndim, ndim)        ! original matrix
  REAL(kind = 8) :: trans_matrix(ndim, ndim)  ! transposed matrix

!Note : A double transpose of  a square matrix is the original matrix



  do iy= 1, ndim
     do ix= 1, ndim

        trans_matrix(iy,ix) = matrix(ix,iy)

     end do
  end do

  matrix = trans_matrix;

end subroutine mattrans

subroutine             mattrans_g(nx, ny, morig, mtran)
  implicit             none
  integer           :: ix, iy, nx, ny
  REAL(kind = 8) :: morig(nx, ny)            ! original matrix
  REAL(kind = 8) :: mtran(ny, nx)            ! transposed matrix

  do iy= 1, ny
     do ix= 1, nx

        mtran(iy,ix) = morig(ix,iy)

     end do
  end do
end subroutine mattrans_g


      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
      DOUBLE PRECISION ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.d0/MBIG)
!     PARAMETER (MBIG=4000000.d0,MSEED=1618033.d0,MZ=0.d0,FAC=1.d0/MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
    END FUNCTION ran3


        
FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      DOUBLE PRECISION ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,      &
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)



END  FUNCTION ran2



  
FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      DOUBLE PRECISION ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1.d0/IM,IQ=127773,IR=2836, &
      NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
     
END  FUNCTION ran1  

