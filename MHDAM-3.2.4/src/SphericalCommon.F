#include "CONSTANTS.H"
#include "SphericalIntegrator.H"


c Module that stores data for polar (2D) and spherical (3D) grids
      MODULE SphericalData
      
      IMPLICIT NONE
      
      ! number of ghost cells
      INTEGER, PARAMETER, PUBLIC :: m_nGS = 4
      
      REAL*8  :: m_rlo(0:CH_SPACEDIM-1), m_rhi(0:CH_SPACEDIM-1)
      INTEGER :: m_ilo(0:CH_SPACEDIM-1), m_ihi(0:CH_SPACEDIM-1)
      INTEGER    m_numlevels      
      INTEGER :: m_size(0:CH_SPACEDIM-1)
      INTEGER, DIMENSION (:), ALLOCATABLE :: m_ref_ratios
      INTEGER, DIMENSION (:), ALLOCATABLE :: m_ref_factor
      REAL*8 , DIMENSION (:), ALLOCATABLE :: m_inv_factor ! = one/m_ref_factor
                  
      ! array of dr
      REAL*8, DIMENSION (:,:), ALLOCATABLE :: m_dr 
      ! distances to cells, cell center for 'i'th cell is located at m_rn(i)+0.5*m_dr(i)
      REAL*8, DIMENSION (:,:), ALLOCATABLE :: m_rn      
      ! distances to cell cells
      REAL*8, DIMENSION (:,:), ALLOCATABLE :: m_rc      
      
      
      ! array of dphi
      REAL*8, DIMENSION (:), ALLOCATABLE :: m_dphi ! phi will always be uniform (or at least now)
      ! defined for cell centers trigonometric functions (sin and cos of (j+half)*dphi)
      REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: m_phic 
      ! defined for nodes trigonometric functions (sin and cos of j*dphi)
      REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: m_phin
      
      ! array of dtheta
      REAL*8, DIMENSION (:,:), ALLOCATABLE :: m_dtheta
      ! defined for cell centers trigonometric functions (sin and cos of (j+half)*dtheta)
      REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: m_thetac 
      ! defined for nodes trigonometric functions (sin and cos of j*dtheta)
      REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: m_thetan
      
      contains 
      
! Prepare grid where dr is constant
      subroutine PrepareUniformR(a_R0, a_Rmax, rsize)       
        implicit none      
        REAL*8 :: a_R0, a_Rmax,dr
        integer :: rsize
                                      
        dr = (a_Rmax-a_R0)/rsize      
        m_dr(:,0) = dr        
                    
      end subroutine PrepareUniformR
      
      
! Prepare polar grid in which all cells have three edges of the same length
      subroutine PrepareUniformCells(a_R0, a_Rmax, rsize, a_dphi)       
        implicit none        
        REAL*8  :: a_R0, a_Rmax, a_dphi
        integer :: rsize        
        
        integer :: j
       
        m_dr(0,0) = a_R0*a_dphi
        m_rn(0,0) = a_R0        
        do j=1,rsize-1
          m_rn(j,0) = m_rn(j-1,0)+m_dr(j-1,0)
          m_dr(j,0) = m_rn(j,0)*a_dphi          
        enddo          
      
      end subroutine PrepareUniformCells
      
! Prepare grid with exponential distribution
      subroutine PrepareExponentialR(a_R0, a_Rmax, imin, imax, a_coeff)       
        implicit none        
        REAL*8  :: a_R0, a_Rmax, a_coeff
        integer :: imin, imax
        
        REAL*8  :: rsizef
        integer :: i
                
        rsizef = real(imax-imin+1)
               
        m_rn(imin,0) = a_R0        
        do i=imin+1,imax+1
          m_rn(i,0) = a_R0+(a_Rmax-a_R0)*
     &    (exp(a_coeff*real(i-imin)/rsizef)-1D0)/(exp(a_coeff)-1D0)
          m_dr(i-1,0) = m_rn(i,0)-m_rn(i-1,0)
        enddo          
                    
      end subroutine PrepareExponentialR

      
! Fill m_dr, m_rc, m_rn arrays assuming that m_dr(0:r0_max) is filled.      
      subroutine FillRarrays()
        implicit none      
        integer :: rsize     
        integer sizemax,l,j,ref_r
        
                         
        
        rsize = m_size(0)
                
        m_rn(0,0) = m_rlo(0)
        m_rc(0,0) = m_rlo(0)+half*m_dr(0,0)
        do j=1,rsize-1
          m_rn(j,0) = m_rn(j-1,0)+m_dr(j-1,0)          
          m_rc(j,0) = m_rn(j,0)+half*m_dr(j,0)
        enddo  
        m_rn(rsize,0) = m_rn(rsize-1,0)+m_dr(rsize-1,0)
        
        sizemax = rsize
        do l=1,m_numlevels
          ref_r   = m_ref_ratios(l-1)
          sizemax = sizemax*ref_r
          m_rn(0,l) = m_rlo(0)
          do j=0,sizemax-1
            m_dr(j,l)   = m_dr(j/ref_r,l-1)/ref_r
            m_rc(j,l)   = m_rn(j,l)+half*m_dr(j,l)
            m_rn(j+1,l) = m_rn(j,l)+m_dr(j,l)        
          enddo      
        enddo        
                
! Filling ghost cells              
        sizemax = rsize
        do l=0,m_numlevels        
          m_dr(-4:-1,l) = m_dr(0,l)
          m_rn(-4,l) = m_rlo(0)-4.d0 *m_dr(0,l)
          m_rc(-4,l) = m_rlo(0)-3.5d0*m_dr(0,l)        
          do j=-3,-1
            m_rn(j,l) = m_rn(j-1,l)+ m_dr(j-1,l)
            m_rc(j,l) = m_rn(j,l)  + half*m_dr(j,l)          
          enddo
          m_dr(sizemax:sizemax+3,l) = m_dr(sizemax-1,l)
          do j=sizemax,sizemax+3
            m_rc(j,l)   = m_rn(j,l)+half*m_dr(j,l)
            m_rn(j+1,l) = m_rn(j,l)+m_dr(j,l)                
          enddo        
          sizemax = sizemax*m_ref_ratios(l)        
        enddo        
                      
      end subroutine FillRArrays
      
            
            
! Fill m_phic, m_phin arrays assuming that m_phi is filled.      
      subroutine FillPhiArrays()
        implicit none              
        integer i,j,ref_factor,jmin,jmax
        
        REAL*8 phi,dphi,dphiLo,dphiHi,angle
                        
        m_dphi(m_ilo(1)-4:m_ilo(1)-1) = m_dphi(m_ilo(1))
        m_dphi(m_ihi(1)+1:m_ihi(1)+4) = m_dphi(m_ihi(1))
       
        ref_factor = 1
        do i=0,m_numlevels  
          phi  = m_rlo(1)
          jmin = m_ilo(1)*ref_factor
          jmax = m_size(1)*ref_factor-1
          do j = jmin,  jmax            
            dphi = m_dphi(j/ref_factor)/ref_factor
            m_phin(j,ISIN,i) = sin(phi)
            m_phin(j,ICOS,i) = cos(phi)            
            phi = phi + half*dphi
            m_phic(j,ISIN,i) = sin(phi)
            m_phic(j,ICOS,i) = cos(phi)
            phi = phi + half*dphi            
          end do  
          m_phin(jmax+1,ISIN,i) = sin(phi)
          m_phin(jmax+1,ICOS,i) = cos(phi)            
           
          ref_factor = ref_factor*m_ref_ratios(i)                  
        enddo

        
! Filling ghost cells                                      
        dphiLo   = m_dphi(m_ilo(1))        
        dphiHi   = m_dphi(m_ihi(1))    
        ref_factor = 1
        do i=0,m_numlevels          
          jmin = m_ilo(1)*ref_factor
          jmax = m_size(1)*ref_factor-1
          do j = jmin-4,  jmin-1
            angle = (j+0.5D0)*dphiLo + m_rlo(1)
            m_phic(j,ISIN,i) = sin(angle)
            m_phic(j,ICOS,i) = cos(angle)
            
            angle = j        *dphiLo + m_rlo(1)
            m_phin(j,ISIN,i) = sin(angle)
            m_phin(j,ICOS,i) = cos(angle)
          end do  
          do j=jmax+1,jmax+4
            angle = (j-jmax-0.5D0)*dphiHi + m_rhi(1)
            m_phic(j,ISIN,i) = sin(angle)
            m_phic(j,ICOS,i) = cos(angle)
            angle = (j-jmax      )*dphiHi + m_rhi(1)
            m_phin(j+1,ISIN,i) = sin(angle)
            m_phin(j+1,ICOS,i) = cos(angle)
          end do                      
          dphiLo     = dphiLo/m_ref_ratios(i)        
          dphiHi     = dphiHi/m_ref_ratios(i)                  
          ref_factor = ref_factor*m_ref_ratios(i)                  
        enddo                
      
      end subroutine FillPhiArrays
      
      
      subroutine FillThetaArrays()
        implicit none              
#if CH_SPACEDIM == 3              
        integer level,j,ref_factor,jmin,jmax,ref_r
                
        REAL*8 theta,dtheta
                                
        level = 0
        
c        m_dtheta(m_ilo(2)-1,level) = m_dtheta(m_ilo(2)  ,level)
c        m_dtheta(m_ilo(2)-2,level) = m_dtheta(m_ilo(2)+1,level)
c        m_dtheta(m_ilo(2)-3,level) = m_dtheta(m_ilo(2)+2,level)
c        m_dtheta(m_ilo(2)-4,level) = m_dtheta(m_ilo(2)+3,level)
                
c        m_dtheta(m_ihi(2)+1,level) = m_dtheta(m_ihi(2)  ,level)
c        m_dtheta(m_ihi(2)+2,level) = m_dtheta(m_ihi(2)-1,level)
c        m_dtheta(m_ihi(2)+3,level) = m_dtheta(m_ihi(2)-2,level)
c        m_dtheta(m_ihi(2)+4,level) = m_dtheta(m_ihi(2)-3,level)
                
        theta= m_rlo(2)     
        jmin = m_ilo(2)
        jmax = m_size(2)
        do j = jmin,  jmax            
          dtheta = m_dtheta(j,level)
          m_thetan(j,ISIN,level) = sin(theta)
          m_thetan(j,ICOS,level) = cos(theta)            
          theta = theta + half*dtheta
          m_thetac(j,ISIN,level) = sin(theta)
          m_thetac(j,ICOS,level) = cos(theta)
          theta = theta + half*dtheta            
        end do  
        m_thetan(jmax+1,ISIN,level) = sin(theta)
        m_thetan(jmax+1,ICOS,level) = cos(theta)                        
        
        ! Higher levels                
        ref_factor = m_ref_ratios(level)
        do level=1,m_numlevels  
          theta = m_rlo(2)     
          ref_r = m_ref_ratios(level-1)
          jmin  = m_ilo(2)*ref_factor
          jmax  = m_size(2)*ref_factor-1
          do j = jmin,  jmax            
            dtheta = m_dtheta(j/ref_r,level-1)/ref_r
            m_dtheta(j,level) = dtheta
            m_thetan(j,ISIN,level) = sin(theta)
            m_thetan(j,ICOS,level) = cos(theta)            
            theta = theta + half*dtheta
            m_thetac(j,ISIN,level) = sin(theta)
            m_thetac(j,ICOS,level) = cos(theta)
            theta = theta + half*dtheta            
          end do  
          m_thetan(jmax+1,ISIN,level) = sin(theta)
          m_thetan(jmax+1,ICOS,level) = cos(theta)                        
          ref_factor = ref_factor*m_ref_ratios(level)                  
        enddo
        
! Filling ghost cells              
        
        ref_factor = 1
        do level=0,m_numlevels                    
          jmin = m_ilo(2)*ref_factor
          jmax = m_size(2)*ref_factor-1
          do j = 1, 4                                
            m_thetac(jmin-j,ISIN,level) = -m_thetac(jmin+j-1,ISIN,level)
            m_thetac(jmin-j,ICOS,level) =  m_thetac(jmin+j-1,ICOS,level)                                   
            m_thetan(jmin-j,ISIN,level) = -m_thetan(jmin+j,ISIN,level)
            m_thetan(jmin-j,ICOS,level) =  m_thetan(jmin+j,ICOS,level)
            
            m_dtheta(jmin-j,level) = m_dtheta(jmin+j-1,level)

            m_thetac(jmax+j,ISIN,level) = -m_thetac(jmax-j+1,ISIN,level)
            m_thetac(jmax+j,ICOS,level) =  m_thetac(jmax-j+1,ICOS,level)
            m_thetan(jmax+1+j,ISIN,level) = -m_thetan(jmax-j+1,ISIN,level)
            m_thetan(jmax+1+j,ICOS,level) =  m_thetan(jmax-j+1,ICOS,level)            
            
            m_dtheta(jmax+j,level) = m_dtheta(jmax-j+1,level)
          end do  
          
          ref_factor = ref_factor*m_ref_ratios(level)                  
        enddo
#endif                      
      end subroutine FillThetaArrays

c delete this subroutine          
            
      
      END MODULE SphericalData

!     Init all necessary data
      subroutine InitSpherical(nlevels, lo, hi, sizes, ref_ratios, 
     &  rspacing, expCoeff)  
     
      USE SphericalData
      implicit none
          
      REAL*8 :: lo(0:CH_SPACEDIM-1), hi(0:CH_SPACEDIM-1), expCoeff 
      integer :: sizes(0:CH_SPACEDIM-1)      
      integer :: nlevels, ref_ratios(0:nlevels), rspacing
      
      integer sizemax,i
      
      REAL*8 dangle
      
      m_size(:) = sizes(:)
      m_rlo(:)  = lo(:)
      m_rhi(:)  = hi(:)
      do i=0,CH_SPACEDIM-1
        m_ilo(i) = 0
        m_ihi(i) = sizes(i)-1
      enddo
      
      m_numlevels = nlevels


                  
      ALLOCATE(
     &   m_ref_ratios(0:nlevels),
     &   m_ref_factor(0:nlevels),
     &   m_inv_factor(0:nlevels))
        
      m_ref_ratios(:) = ref_ratios(:)
      m_ref_factor(0) = 1
      m_inv_factor(0) = one
            
      do i=1,nlevels
        m_ref_factor(i) = m_ref_factor(i-1)*m_ref_ratios(i-1)
        m_inv_factor(i) = m_inv_factor(i-1)/m_ref_ratios(i-1)
      enddo
            
c Allocate arrays in 'r' direction      
      sizemax = sizes(0)
      do i=0,nlevels-1
        sizemax = sizemax*ref_ratios(i)
      enddo
      
      ALLOCATE(m_dr(-m_nGS:sizemax+m_nGS-1,0:nlevels),
     &         m_rc(-m_nGS:sizemax+m_nGS-1,0:nlevels),
     &         m_rn(-m_nGS:sizemax+m_nGS,0:nlevels))
          
      if (rspacing == 0) then
        call PrepareUniformR(lo(0), hi(0), sizes(0))
      else if (rspacing == 1) then
        dangle = (hi(1)-lo(1))/sizes(1)      
        call PrepareUniformCells(lo(0), hi(0), sizes(0), dangle)
      else if (rspacing == 2) then   
        call PrepareExponentialR(lo(0), hi(0), 0, sizes(0)-1, expCoeff)
      endif      
      
      call FillRarrays()
            
c Allocate arrays in 'phi' direction                           
      sizemax = sizes(1)
      do i=0,nlevels-1
        sizemax = sizemax*ref_ratios(i)
      enddo
      
      ALLOCATE(m_dphi(m_ilo(1)-m_nGS:m_ihi(1)+m_nGS),
     &         m_phic(-m_nGS:sizemax+m_nGS-1,0:1,0:nlevels), 
     &         m_phin(-m_nGS:sizemax+m_nGS  ,0:1,0:nlevels))
      
      dangle = (hi(1)-lo(1))/sizes(1)
      m_dphi(:) = dangle
      call FillPhiArrays()      
         
#if CH_SPACEDIM == 3
      sizemax = sizes(2)
      do i=0,nlevels-1
        sizemax = sizemax*ref_ratios(i)
      enddo

c Allocate arrays in 'theta' direction                                 
      ALLOCATE(m_dtheta(-m_nGS:sizemax+m_nGS-1,0:nlevels),
     &         m_thetac(-m_nGS:sizemax+m_nGS-1,0:1,0:nlevels), 
     &         m_thetan(-m_nGS:sizemax+m_nGS  ,0:1,0:nlevels))
      
      dangle = (hi(2)-lo(2))/sizes(2)
      m_dtheta(:,0) = dangle      
      call FillThetaArrays()      
#endif
      
      
       
      return 
      end


      subroutine FinalizeSpherical()  
      USE SphericalData
      implicit none
      
      DEALLOCATE(m_ref_ratios)
      DEALLOCATE(m_dr,m_rn,m_rc,m_dphi,m_phic,m_phin)
#if CH_SPACEDIM == 3      
      DEALLOCATE(m_dtheta,m_thetac,m_thetan)
#endif      
      return 
      end
      
      subroutine setgridspacing(a_du,a_dir)  
      USE SphericalData
      implicit none
            
      integer :: a_dir
      REAL*8  :: a_du(0:m_size(a_dir)-1)       
      
      if (a_dir == 0) then
        m_dr(0:m_size(a_dir)-1,0) = a_du(0:m_size(a_dir)-1)
        call FillRarrays()      
      endif
      
      if (a_dir == 1) then
        m_dphi(0:m_size(a_dir)-1) = a_du(0:m_size(a_dir)-1)
        call FillPhiArrays()      
      endif
      
      if (a_dir == 2) then
        m_dtheta(0:m_size(a_dir)-1,0) = a_du(0:m_size(a_dir)-1)
        call FillThetaArrays()      
      endif      
            
      return 
      end
      