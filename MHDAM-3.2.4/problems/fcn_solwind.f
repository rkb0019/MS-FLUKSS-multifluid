      subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
      implicit none
      
      real*8 VELE,RHOE,PRESE,rfcn
      common /fcnsw/ VELE,RHOE,PRESE,rfcn

      integer n,ldfjac,iflag
      double precision x(n),fvec(n),fjac(ldfjac,n)
      real*8 re,ue,pe,gam,r
      re = RHOE
      pe = PRESE
      ue = VELE
      gam = 5D0/3D0
      r   = rfcn
      
      if(iflag.EQ.1) then
         fvec(1)=x(1)*x(2)*r*r-re*ue
         fvec(2)=(x(2)*x(2)-ue*ue)/2.0+gam/(gam-1.0)*(x(3)/x(1)-pe/re)
         fvec(3)=x(3)-pe*(x(1)/re)**gam
       elseif(iflag.EQ.2) then
         fjac(1,1)=x(2)*r*r
         fjac(1,2)=x(1)*r*r
         fjac(1,3)=0.0
         fjac(2,1)=-gam/(gam-1.0)*x(3)/x(1)/x(1)
         fjac(2,2)=x(2)
         fjac(2,3)=gam/(gam-1.0)/x(1)
         fjac(3,1)=-gam*pe/re*(x(1)/re)**(gam-1.0)
         fjac(3,2)=0.0
         fjac(3,3)=1.0
      endif
      return
      end
