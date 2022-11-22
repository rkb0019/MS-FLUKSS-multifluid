! This code follows the trajectories of a number of mega neutral particles
! and allows the motion to be altered through charge exchange events.
!
! Work in polar coordinates: 
!  theta=0 is positive z-axis, theta=pi/2 is positive x-axis, 
!  y is the ignorable coordinate. Assume symmetry across z-axis for theta > pi.
!
!  theta = acos(z/r)
!
! 3D velocity can be projected onto xz-plane by
!  vxr = (vx*xi + vy*yi)/sqrt(xi*xi + yi*yi)
!
! particle array: massi, xi, yi, zi, vxi, vyi, vzi, region
!  
! SI units are employed throughout this code.
!
! Here we do not use the Lipatov et al 1998 approach to counting down the
! charge-exchange time. Instead we use Prob(chex) = n_p sigma_ex ds.
!
!***********************************************************************
!***********************************************************************
!
MODULE mpi_defs_2D
  IMPLICIT NONE
  INCLUDE "mpif.h"
  INTEGER :: status(MPI_STATUS_SIZE), nprocMPI, my_rank, nthreads, my_thread
!$OMP threadprivate (nprocMPI, my_rank, nthreads, my_thread)
END MODULE mpi_defs_2D
!
!***********************************************************************
!

MODULE global_2D
!
  USE mpi_defs_2D
!
! Global variable declarations
!
  IMPLICIT NONE
  INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND(14,307)
  
  type charge_exchange_event
    INTEGER :: level
    INTEGER :: i
    INTEGER :: j
    REAL(KIND=dp) :: source_px, source_pz, source_mvsq
    REAL(KIND=dp) :: source_m
  end type charge_exchange_event
  
  type particle_storage
    INTEGER :: my_npart ! number of particles actually stored
    REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: particles
  end type particle_storage
  
  type node_workload
    INTEGER :: procID
    INTEGER :: part_disbalance
  end type node_workload
    
  
!
  
  INTEGER, PARAMETER, PUBLIC :: lenstr = 60
  INTEGER, PARAMETER, PUBLIC :: nsectors = 8
  REAL(KIND=dp), PARAMETER, PUBLIC :: au = 1.5e11_dp            ! m
  REAL(KIND=dp), PARAMETER, PUBLIC :: r_E_sq = au*au            ! m^2
  REAL(KIND=dp), PARAMETER, PUBLIC :: beta_ph_E = 9e-8_dp     ! s^{-1}
  REAL(KIND=dp), PUBLIC :: phi_rmax_sq
  REAL(KIND=dp), PARAMETER, PUBLIC :: pi = &
       3.1415926535897932384626433832795_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: twopi = 6.28318530717959_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: secs_per_year = 31536000.0_dp  
  REAL(KIND=dp), PUBLIC :: LISM_nH, LISM_vH, LISM_TH, k_B
  REAL(KIND=dp), PUBLIC :: LISM_v_thermal, LISM_density, boundary_layer
  REAL(KIND=dp), PUBLIC :: total_volume
  REAL(KIND=dp), PUBLIC :: neutral_mass, mega_number, macro_mass
  REAL(KIND=dp), PUBLIC :: time, starttime, end_time, sourcetime
  REAL(KIND=dp), PUBLIC :: timescale, timestep, my_tot_mass_loss
  REAL(KIND=dp), PUBLIC :: test_dens_max, test_temp_max, test_ux_max, test_uz_max
  REAL(KIND=dp), PUBLIC :: test_dens_min, test_temp_min, test_ux_min, test_uz_min
  REAL(KIND=dp), PUBLIC :: my_min_dt, my_max_dt, my_ave_dt
  INTEGER,       PUBLIC :: my_ntotal, my_n_chex
  INTEGER,       PUBLIC :: my_npart_boundary
  INTEGER(8),    PUBLIC :: npart
  INTEGER,       PUBLIC :: verbosity
  LOGICAL,       PUBLIC :: restart, steady_state_output, photoionize
  CHARACTER (LEN = lenstr), PUBLIC :: rawfile, gridfile, load_file, sourcefile
  
  INTEGER,       PUBLIC :: pout_unit
  CHARACTER (LEN = lenstr), PUBLIC ::pout_name
  
  INTEGER,       PUBLIC :: my_nchex_events, my_nchex_current
  type(charge_exchange_event), DIMENSION (:), ALLOCATABLE, PUBLIC :: my_chex_events
  type(particle_storage),      DIMENSION (:), ALLOCATABLE, PUBLIC :: my_particles
  type(node_workload),         DIMENSION (:), ALLOCATABLE, PUBLIC :: mpi_balancing
  

!$OMP threadprivate (my_nchex_events, my_nchex_current, my_chex_events)
!$OMP threadprivate (my_min_dt, my_max_dt, my_ave_dt)
!$OMP threadprivate (my_n_chex,my_tot_mass_loss)
!$OMP threadprivate (pout_name,pout_unit)


  REAL(KIND=dp), PUBLIC :: int_r(101,101), int_phi(101,101,101), drH, drp, dphi
  REAL(KIND=dp), PUBLIC :: int_Maxwellian(0:400)




  
!
! source grid variables
  REAL(KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: source_px, source_pz
  REAL(KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: source_m, source_mvsq
  INTEGER, DIMENSION (:,:,:), ALLOCATABLE :: source_nchex
  REAL(KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: g_source_m, g_source_mvsq
  REAL(KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: g_source_px, g_source_pz
  INTEGER, DIMENSION (:,:,:), ALLOCATABLE :: g_source_nchex
  REAL(KIND=dp), DIMENSION (:), ALLOCATABLE :: grids_dx, grids_boundaries_x, grids_boundaries_z
  INTEGER,       DIMENSION (:), ALLOCATABLE :: grids_nx, grids_nz
  INTEGER,       PUBLIC :: ngrids, grids_nx_tot, grids_nxmax, grids_nzmax
  INTEGER,       PUBLIC :: total_ndata
!
! plasma grid parameters  (same grid as sources)
  REAL(KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: p_dens, p_temp
  REAL(KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: p_ux, p_uz
  INTEGER, DIMENSION (:,:,:), ALLOCATABLE :: p_region
  REAL(KIND=dp), PUBLIC :: p_xmax, p_xmax_sq, p_zmin, p_zmax
!
! neutral grid parameters  (same grid as sources, except we also have REGIONS)
  REAL(KIND=dp), DIMENSION (:,:,:,:), ALLOCATABLE :: grid_dens, grid_temp
  REAL(KIND=dp), DIMENSION (:,:,:,:), ALLOCATABLE :: grid_ux, grid_uz
  INTEGER, DIMENSION (:,:,:,:), ALLOCATABLE :: grid_npart

END MODULE global_2D



!
!***********************************************************************
!***********************************************************************
!
MODULE neutrals_routines1_2D
!
! Lowest level subroutines
!
  USE global_2D
  USE random_gen_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_vr_phi, init_Maxwellian_vr
  PUBLIC :: get_vr_phi, get_Maxwellian_vr, inject_maxwellian
  PUBLIC :: plasma, exchange, collect_sources, source_chex, source_phi, collect_sources_node
  PUBLIC :: pout,init_pout

CONTAINS

  SUBROUTINE init_pout
    pout_unit = 100 + my_rank*nprocMPI + my_thread
    write  ( pout_name, '(a6,i2,a2,i2)' ) 'fort.n', my_rank, '.t', my_thread
    !open(unit = pout_unit,file=pout_name,form='formatted',status='unknown',action='write')
    !close(unit = pout_unit,status='delete')
    !open(unit = pout_unit,file=pout_name,status='unknown')
      !rewind(unit = pout_unit)
  end SUBROUTINE init_pout

  SUBROUTINE pout(str)
    character (len=*) str
    
 !   open(unit = pout_unit,file=pout_name,form='formatted',status='unknown', ACCESS = 'APPEND')
 !   write (pout_unit,*) str
 !   close(pout_unit)          
  end SUBROUTINE pout
  
  SUBROUTINE init_vr_phi  
    REAL(KIND=dp) :: rH, rp, int_p,phi
    INTEGER :: i,j,k
       
    drH = 0.1_dp
    drp = 0.03_dp
    !       dphi = 2.0_dp*pi/100.0_dp
    dphi = pi/100.0_dp

    DO i = 1, 101
      rH = (i-0.5_dp)*drH

      int_r(i,1) = 0.0_dp
      DO j = 2, 101
         rp = (j-0.5_dp)*drp

         int_p = 0.0_dp
         DO k = 2, 101
            phi = (k-0.5)*dphi
            int_p = int_p + &
    !            SQRT(rH*rH + rp*rp - 2.0_dp*rH*rp*COS(phi)) * rp * EXP(-rp*rp) * dphi
        SQRT(rH*rH + rp*rp - 2.0_dp*rH*rp*COS(phi)) * rp*rp * EXP(-rp*rp) * dphi
         END DO

         int_r(i,j) = int_r(i,j-1) + int_p*drp

      END DO
      int_r(i,:) = int_r(i,:)/int_r(i,101)
    END DO

    DO i = 1, 101
      rH = (i-0.5_dp)*drH
      
      DO j = 1, 101
         rp = (j-0.5_dp)*drp

         int_phi(i,j,1) = 0.0_dp
         DO k = 2, 101
            phi = (k-0.5_dp)*dphi
            int_phi(i,j,k) = int_phi(i,j,k-1) + &
    !            SQRT(rH*rH + rp*rp - 2.0_dp*rH*rp*COS(phi)) * rp * EXP(-rp*rp) * dphi
        SQRT(rH*rH + rp*rp - 2.0_dp*rH*rp*COS(phi)) * rp*rp * EXP(-rp*rp) * dphi
         END DO

         int_phi(i,j,:) = int_phi(i,j,:)/int_phi(i,j,101)

      END DO
    END DO
    
  end SUBROUTINE init_vr_phi
  
  SUBROUTINE get_vr_phi(vH,vT,vr,phi)
! Here we return the speed and direction (in the plasma velocity frame)
! of the new neutral after ch-ex. Follows Section 7.4 of Lipatov's book (2002)
    REAL(KIND=dp), INTENT(in)  :: vH, vT
    REAL(KIND=dp), INTENT(out) :: vr, phi
    !REAL(KIND=dp), SAVE :: int_r(101,101), int_phi(101,101,101), drH, drp, dphi
    !REAL(KIND=dp) :: rH, rp, int_p
    REAL(KIND=dp) :: j1, k1, dprand
    INTEGER :: i, j, k, dj, dk, jlow, jhigh, klow, khigh
    
    
    i = FLOOR(vH/vT * (100.0_dp*drH)) + 1
    i = MIN(i,101)
    i = MAX(i,1)

    CALL rand_number(dprand)

! use bisection to find speed
    jlow = 1
    jhigh = 101
    dj = jhigh - jlow
    j = 25
    DO WHILE (dj .GT. 1)
       j1 = int_r(i,j)
       IF (j1 .GT. dprand) THEN
          jhigh = j
       ELSE
          jlow = j
       END IF
       j = FLOOR((jhigh+jlow)/2.0_dp)
       dj = jhigh - jlow
    END DO
    j = jhigh

    vr = (j-1.3_dp)*vT*drp

    CALL rand_number(dprand)
! use bisection to find angle
       klow = 1
       khigh = 101
       dk = khigh - klow
       k = 50
       DO WHILE (dk .GT. 1)
          k1 = int_phi(i,j,k)
          IF (k1 .GT. dprand) THEN
             khigh = k
          ELSE
             klow = k
          END IF
          k = FLOOR((khigh+klow)/2.0_dp)
          dk = khigh - klow
       END DO
       k = khigh
       phi = (k-1.5_dp)*dphi

  END SUBROUTINE get_vr_phi
  
  
  SUBROUTINE init_Maxwellian_vr()
    INTEGER :: i
    REAL(KIND=dp) :: vr, dvr,tot_int
    
    dvr = 0.01_dp     
       
    int_Maxwellian(0) = 0.0_dp
    DO i = 1, 400
      vr = (i-0.5_dp)*dvr

      int_Maxwellian(i) = int_Maxwellian(i-1) + vr*vr * EXP(-vr*vr) * dvr

    END DO

    tot_int = int_Maxwellian(400)
    int_Maxwellian = int_Maxwellian / tot_int

  end SUBROUTINE init_Maxwellian_vr
  
!
!***********************************************************************
!
  SUBROUTINE get_Maxwellian_vr(vr)
! Here we randomly select a speed from a Maxwellian distribution.
! Follows Section 7.3 of Lipatov's book (2002)
    REAL(KIND=dp), INTENT(out) :: vr
    !REAL(KIND=dp), SAVE :: int_Maxwellian(0:400)
    REAL(KIND=dp) :: dvr, tot_int, i1, dprand
    INTEGER :: i, di, ilow, ihigh    

    dvr = 0.01_dp
    
    CALL rand_number(dprand)

! use bisection to find speed
    ilow = 0
    ihigh = 400
    di = ihigh - ilow
    i = 100
    DO WHILE (di .GT. 1)
       i1 = int_Maxwellian(i)
       IF (i1 .GT. dprand) THEN
          ihigh = i
       ELSE
          ilow = i
       END IF
       i = FLOOR((ihigh+ilow)/2.0_dp)
       di = ihigh - ilow
    END DO
    i = ihigh

    vr = (i-0.5_dp)*dvr

  END SUBROUTINE get_Maxwellian_vr
!
!***********************************************************************
!
  SUBROUTINE inject_maxwellian(x1,x2,z1,z2,v_thermal,my_number)
! inject a uniform Maxwellian distributions into a prescibed rectangular region
    REAL(KIND=dp), INTENT (in) :: x1, x2, z1, z2, v_thermal
    REAL(KIND=dp) :: mass0, x0, z0, vxy, vr, delta_x, delta_z, theta, phi
    REAL(KIND=dp), DIMENSION(4) :: dprand
    INTEGER :: i, my_number, ind
    REAL(KIND=dp) :: vx, vy, vz
 
    delta_x = x2 - x1
    delta_z = z2 - z1

    DO i = 1, my_number

       CALL rand_number(dprand(1))
       CALL rand_number(dprand(2))
       CALL rand_number(dprand(3))
       CALL rand_number(dprand(4))
       
       x0 = x1 + (1.0_dp-dprand(1))*delta_x
       z0 = z1 + dprand(2)*delta_z
       
       ind = my_particles(my_thread)%my_npart+i

       my_particles(my_thread)%particles(ind,1) = x0
       my_particles(my_thread)%particles(ind,2) = 0.0_dp
       my_particles(my_thread)%particles(ind,3) = z0

! Inject neutrals to have a Maxwellian velocity distribution
! This method comes from "The hybrid multiscalesimulation technology"
! by A.S. Lipatov, section 7.2.2.

! Based on 3D speed distribution
       CALL get_Maxwellian_vr(vr) ! obtain Maxwellian vr
       vr = vr * v_thermal

       theta = ACOS(0.99999_dp - 1.9999_dp*dprand(3))
       vz = vr * COS(theta)
       vxy = SQRT(vr*vr - vz*vz)

       phi = 2.0_dp*pi*dprand(4)
       vx = vxy * COS(phi)
       vy = vxy * SIN(phi)

! update particle
       my_particles(my_thread)%particles(ind,4) = vx
       my_particles(my_thread)%particles(ind,5) = vy
       my_particles(my_thread)%particles(ind,6) = vz + LISM_vH

! Maxwellian particles are region 0 (LISM)
       my_particles(my_thread)%particles(ind,7) = 0.1_dp

! LISM particles have unit mass multiplied by OUT-OF-PLANE thickness
       mass0 = x0 * twopi
       my_particles(my_thread)%particles(ind,0) = mass0

    END DO

    my_particles(my_thread)%my_npart = my_particles(my_thread)%my_npart + my_number

  END SUBROUTINE inject_maxwellian
!
!***********************************************************************
!
  SUBROUTINE plasma(xxi,zzi,n_p,T_p,ux,uz,region_p)
    REAL(KIND=dp), INTENT(in)  :: xxi, zzi
    REAL(KIND=dp), INTENT(out) :: n_p, ux, uz, T_p, region_p
    REAL(KIND=dp) :: zi, xyr
    INTEGER :: ii, jj, kkx, kkz, kk, nn
! return plasma number density, velocity components and temp at a given point
! here we DO NOT interpolate, which improves speed

! make sure we don't try to access value OUTSIDE the plasma domain
    xyr = MIN(xxi,0.99999_dp*p_xmax)
    zi = MIN(zzi,0.99999_dp*p_zmax)
    zi = MAX(zi,0.99999_dp*p_zmin)   ! ASSUMES p_zmin < 0

! find which grid we're on and nearest grid point
    kkx = -1
    DO nn = 1, ngrids
       IF (abs(xyr) .GE. grids_boundaries_x(ngrids-nn)) THEN
          kkx = ngrids-nn+1
          EXIT
       END IF
    END DO
    kkz = -1
    DO nn = 1, ngrids
       IF (abs(zi) .GE. grids_boundaries_z(ngrids-nn)) THEN
          kkz = ngrids-nn+1
          EXIT
       END IF
    END DO
    kk = max(kkx,kkz)
    ii = FLOOR(xyr / grids_dx(kk)) + 1
    jj = FLOOR(zi / grids_dx(kk)) + 1

    n_p = p_dens(kk,ii,jj)
    T_p = p_temp(kk,ii,jj)
    ux = p_ux(kk,ii,jj)
    uz = p_uz(kk,ii,jj)
    region_p = p_region(kk,ii,jj)

! for TESTING
!!$    test_dens_max = MAX(test_dens_max,n_p)
!!$    test_temp_max = MAX(test_temp_max,T_p)
!!$    test_ux_max   = MAX(test_ux_max,ux)
!!$    test_uz_max   = MAX(test_uz_max,uz)
!!$    test_dens_min = MIN(test_dens_min,n_p)
!!$    test_temp_min = MIN(test_temp_min,T_p)
!!$    test_ux_min   = MIN(test_ux_min,ux)
!!$    test_uz_min   = MIN(test_uz_min,uz)

  END SUBROUTINE plasma
!
!***********************************************************************
!
  SUBROUTINE exchange(massi,xi,yi,zi,vx,vy,vz,region, &
                      ux,uz,v_p_th,xn,yn,zn,vxn,vyn,vzn)
    REAL(KIND=dp), INTENT (in) :: massi, xi, yi, zi, vx, vy, vz, region
    REAL(KIND=dp), INTENT (in) :: ux, uz, v_p_th
    REAL(KIND=dp), INTENT (out) :: xn, yn, zn, vxn, vyn, vzn
    REAL(KIND=dp) :: xr, vxr, vxy_sq, vyr_sq, dprand
    REAL(KIND=dp) :: vv1, vv2, vv3, vv, vvsq, ww1, ww2, ww3, pp1, pp2, p, vdotw
    REAL(KIND=dp) :: rotvv1, rotvv2, rotvv3, wprojv1, wprojv2, wprojv3
    REAL(KIND=dp) :: r1_1, r1_2, r1_3, rr, r2_1, r2_2, r2_3
    REAL(KIND=dp) :: tp, cp, sp, vvnew1, vvnew2, vvnew3, omega, vr, phi

    my_n_chex = my_n_chex + 1

! axisymmetric projection
    xr = SQRT(xi*xi + yi*yi)
    vxr = (vx*xi + vy*yi) / xr
    vxy_sq = vx*vx + vy*vy
    vyr_sq = ABS(vxy_sq - vxr*vxr)

! NOTE: MASS of new particle same as OLD

! rotate vH about z-axis to bring us onto the frame of the plasma
    vv1 = vxr-ux
    vv2 = SQRT(vyr_sq)          ! sign is not important due to symmetry
    vv3 = vz-uz

    vvsq = vv1*vv1+vv2*vv2+vv3*vv3
    vv = SQRT(vvsq)
! include a selection effect due to thermal protons (OF COURSE we need vyr_sq here)
    CALL get_vr_phi(vv,v_p_th,vr,phi)

    pp1 = -vv2
    pp2 = vv1
    p = SQRT(pp1*pp1+pp2*pp2)
             
    pp1 = pp1/p
    pp2 = pp2/p                 ! normalized rotation axis vector; pp3=0

    cp = COS(phi)
    sp = SIN(phi)
    tp = 1.0_dp-cp
    rotvv1 = (tp*pp1*pp1+cp)*vv1 + (tp*pp1*pp2)   *vv2 + (-sp*pp2)*vv3
    rotvv2 = (tp*pp1*pp2)   *vv1 + (tp*pp2*pp2+cp)*vv2 +  (sp*pp1)*vv3
    rotvv3 =        (sp*pp2)*vv1 +       (-sp*pp1)*vv2 +      (cp)*vv3

    ww1 = vr/vv * rotvv1
    ww2 = vr/vv * rotvv2
    ww3 = vr/vv * rotvv3

    vdotw = vv1*ww1+vv2*ww2+vv3*ww3

    wprojv1 = vdotw/vvsq * vv1
    wprojv2 = vdotw/vvsq * vv2
    wprojv3 = vdotw/vvsq * vv3

    r1_1 = ww1 - wprojv1
    r1_2 = ww2 - wprojv2
    r1_3 = ww3 - wprojv3

    rr = SQRT(r1_1*r1_1 + r1_2*r1_2 + r1_3*r1_3)

    r2_1 = rr * pp1                      ! pp has unit length
    r2_2 = rr * pp2
    r2_3 = 0.0_dp                        ! since pp3=0

    CALL rand_number(dprand)
    omega = 2.0_dp*pi * dprand             ! random phase for scatter by phi

    vvnew1 = wprojv1 + r1_1*COS(omega) + r2_1*SIN(omega)
    vvnew2 = wprojv2 + r1_2*COS(omega) + r2_2*SIN(omega)
    vvnew3 = wprojv3 + r1_3*COS(omega) + r2_3*SIN(omega)

    vxn = ux + vvnew1
    vyn =      vvnew2
    vzn = uz + vvnew3

! return post-chex particle back to x-z plane
    xn = xr
    yn = 0.0_dp
    zn = zi

! add ch-ex event to plasma source term
    CALL source_chex(massi,xi,yi,zi,vx,vy,vz,vxn,vyn,vzn)

  END SUBROUTINE exchange
!
!***********************************************************************
!
  SUBROUTINE source_chex(mass,xi,yi,zi,vxold,vyold,vzold,vxnew,vynew,vznew)
    REAL(KIND=dp), INTENT(in)  :: mass, xi, yi, zi, vxold, vzold, vyold
    REAL(KIND=dp), INTENT(in)  :: vxnew, vznew, vynew
    REAL(KIND=dp) :: dvxr, dvz, dvsq, xyr, vxrold
    INTEGER :: ii, jj, nn, kkx, kkz, kk

! add ch-ex source to CARTESIAN grid

! MAKE SURE that abs(zi) and xi are LESS THAN grids_boundaries(ngrids)
    xyr = SQRT(xi*xi + yi*yi)
    IF ( (xyr .LT. grids_boundaries_x(ngrids)) .AND. &
     (ABS(zi) .LT. grids_boundaries_z(ngrids)) ) THEN

! find projection of vxy into xr using dot product
       vxrold = (vxold*xi + vyold*yi) / xyr

! find which grid we're on and nearest grid point
       kkx = -1
       DO nn = 1, ngrids
          IF (abs(xyr) .GE. grids_boundaries_x(ngrids-nn)) THEN
             kkx = ngrids-nn+1
             EXIT
          END IF
       END DO
       kkz = -1
       DO nn = 1, ngrids
          IF (abs(zi) .GE. grids_boundaries_z(ngrids-nn)) THEN
             kkz = ngrids-nn+1
             EXIT
          END IF
       END DO
       kk = max(kkx,kkz)
       ii = FLOOR(xyr / grids_dx(kk)) + 1
       jj = FLOOR(zi / grids_dx(kk)) + 1

       dvxr = vxrold - vxnew
       dvz  = vzold - vznew
       dvsq = (vxold*vxold + vzold*vzold + vyold*vyold &
             - vxnew*vxnew - vznew*vznew - vynew*vynew)

! no mass source since we assume m_p = m_H
! momentum and energy source terms px, pz, m*u^2

       my_nchex_current = my_nchex_current + 1
       my_chex_events(my_nchex_current)%level = kk
       my_chex_events(my_nchex_current)%i = ii
       my_chex_events(my_nchex_current)%j = jj
       my_chex_events(my_nchex_current)%source_m    = 0.0_dp
       my_chex_events(my_nchex_current)%source_px   = dvxr * mass
       my_chex_events(my_nchex_current)%source_pz   = dvz  * mass
       my_chex_events(my_nchex_current)%source_mvsq = dvsq * mass
              
       if (my_nchex_current>=my_nchex_events) then
        call collect_sources_node        
       endif
       

!       source_px(kk,ii,jj)   = source_px(kk,ii,jj)   + dvxr * mass
!       source_pz(kk,ii,jj)   = source_pz(kk,ii,jj)   + dvz  * mass
!       source_mvsq(kk,ii,jj) = source_mvsq(kk,ii,jj) + dvsq * mass
!       source_nchex(kk,ii,jj) = source_nchex(kk,ii,jj) + 1

    END IF

  END SUBROUTINE source_chex
!
!***********************************************************************
!
  SUBROUTINE source_phi(xi,yi,zi,vx,vy,vz,mass_loss)
    REAL(KIND=dp), INTENT(in) :: xi, yi, zi, vx, vz, vy, mass_loss
! Record mass and momentum source from photoionized neutrals
    REAL(KIND=dp) :: xyr, vxr, vsq
    INTEGER :: ii, jj, nn, kkx, kkz, kk

! add ch-ex source to CARTESIAN grid

! MAKE SURE that abs(zi) and xi are LESS THAN grids_boundaries(ngrids)
    xyr = SQRT(xi*xi + yi*yi)
    IF ( (xyr .LT. grids_boundaries_x(ngrids)) .AND. &
    (ABS(zi) .LT. grids_boundaries_z(ngrids)) ) THEN

! find projection of vxy into xr using dot product
       vxr = (vx*xi + vy*yi) / xyr

! find which grid we're on and nearest grid point
       kkx = -1
       DO nn = 1, ngrids
          IF (abs(xyr) .GE. grids_boundaries_x(ngrids-nn)) THEN
             kkx = ngrids-nn+1
             EXIT
          END IF
       END DO
       kkz = -1
       DO nn = 1, ngrids
          IF (abs(zi) .GE. grids_boundaries_z(ngrids-nn)) THEN
             kkz = ngrids-nn+1
             EXIT
          END IF
       END DO
       kk = max(kkx,kkz)
       ii = FLOOR(xyr / grids_dx(kk)) + 1
       jj = FLOOR(zi / grids_dx(kk)) + 1

       vsq = vx*vx + vy*vy + vz*vz

       source_m(kk,ii,jj)    = source_m(kk,ii,jj)    + mass_loss
       source_px(kk,ii,jj)   = source_px(kk,ii,jj)   + vxr * mass_loss
       source_pz(kk,ii,jj)   = source_pz(kk,ii,jj)   + vz  * mass_loss
       source_mvsq(kk,ii,jj) = source_mvsq(kk,ii,jj) + vsq * mass_loss

    END IF

  END SUBROUTINE source_phi

!
!***********************************************************************
!
  
  SUBROUTINE collect_sources_node
  
  INTEGER :: nn,kk,ii,jj
  CHARACTER (LEN = 100) str_buf
  
  write (str_buf,*) 'collect_sources_node, nch_ex : ',my_nchex_current
  call pout(str_buf)      
  
!  my_nchex_current = 0               
!  return
  
  !$omp critical (update_global_sources)
   do nn = 1, my_nchex_current
    kk =  my_chex_events(nn)%level
    ii =  my_chex_events(nn)%i
    jj =  my_chex_events(nn)%j
    
    source_px(kk,ii,jj)   = source_px(kk,ii,jj) + my_chex_events(nn)%source_px
    source_pz(kk,ii,jj)   = source_pz(kk,ii,jj) + my_chex_events(nn)%source_pz
    source_mvsq(kk,ii,jj) = source_mvsq(kk,ii,jj) + my_chex_events(nn)%source_mvsq
    source_nchex(kk,ii,jj) = source_nchex(kk,ii,jj) + 1
   
   enddo       
   !$omp end critical (update_global_sources)
   
   my_nchex_current = 0               
  
  END SUBROUTINE
  
!
!***********************************************************************
!
  SUBROUTINE collect_sources(rtime)
  REAL(KIND=dp) :: gdx, volume, rtime, xx, zz, dx
  INTEGER :: i, j, k, nx, nz, nn, gnx, ii, jj, kkx, kkz, kk, ierr

! NORMALIZE source with respect to cell AREA (third dimension already done)
  DO kk = 1, ngrids
     gnx = grids_nx(kk)
     gdx = grids_dx(kk)
     DO ii = 1, gnx
        volume = pi*gdx**3 * (2.0_dp*ii - 1.0_dp)
        source_m(   kk,ii,:) = source_m(   kk,ii,:) * macro_mass / volume / rtime
        source_px(  kk,ii,:) = source_px(  kk,ii,:) * macro_mass / volume / rtime
        source_pz(  kk,ii,:) = source_pz(  kk,ii,:) * macro_mass / volume / rtime
        source_mvsq(kk,ii,:) = source_mvsq(kk,ii,:) * macro_mass / volume / rtime
     END DO
  END DO

  g_source_px = 0.0_dp;  g_source_pz = 0.0_dp
  g_source_m  = 0.0_dp;  g_source_mvsq = 0.0_dp
  g_source_nchex = 0

  CALL MPI_Reduce(source_m, g_source_m, total_ndata, &
       MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(source_px, g_source_px, total_ndata, &
       MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(source_pz, g_source_pz, total_ndata, &
       MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(source_mvsq, g_source_mvsq, total_ndata, &
       MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(source_nchex, g_source_nchex, total_ndata, &
       MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)


! master will fill in all grids -- this is actually NOT NECESSARY,
! but helps with plotting.
  IF (my_rank .EQ. 0) THEN

     DO k = 1, ngrids
        nx = grids_nx(k)
        nz = grids_nz(k)
        dx = grids_dx(k)
        DO j = -nz+1, nz
           zz = (j-0.5_dp)*dx
           DO i = 1, nx
              xx = (i-0.5_dp)*dx

! find which grid we're on and nearest grid point
              kkx = -1
              DO nn = 1, ngrids
                 IF (abs(xx) .GE. grids_boundaries_x(ngrids-nn)) THEN
                    kkx = ngrids-nn+1
                    EXIT
                 END IF
              END DO
              kkz = -1
              DO nn = 1, ngrids
                 IF (abs(zz) .GE. grids_boundaries_z(ngrids-nn)) THEN
                    kkz = ngrids-nn+1
                    EXIT
                 END IF
              END DO
              kk = max(kkx,kkz)
              ii = FLOOR(xx / grids_dx(kk)) + 1
              jj = FLOOR(zz / grids_dx(kk)) + 1

              g_source_m(k,i,j)     = g_source_m(kk,ii,jj)
              g_source_px(k,i,j)    = g_source_px(kk,ii,jj)
              g_source_pz(k,i,j)    = g_source_pz(kk,ii,jj)
              g_source_mvsq(k,i,j)  = g_source_mvsq(kk,ii,jj)
              g_source_nchex(k,i,j) = g_source_nchex(kk,ii,jj)

           END DO
        END DO
     END DO

  END IF

END SUBROUTINE collect_sources

END MODULE neutrals_routines1_2D
!
!***********************************************************************
!***********************************************************************
!
MODULE neutrals_routines2_2D
!
! Second level subroutines
!

!
  USE global_2D
  USE neutrals_routines1_2D
  USE random_gen_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pass_plasma_grid, return_source_grid, return_neutral_grid
  PUBLIC :: pass_plasma_grid_mhdam, return_source_grid_mhdam, return_neutral_grid_mhdam
  PUBLIC :: mc_init, move, split, write_neutral_grid, write_source, write_plasma, write_plasma_ascii
  PUBLIC :: create_neutral_grid, output_raw, load_raw, balance_npart, mc_neutrals

CONTAINS


  SUBROUTINE pass_plasma_grid_mhdam(kk, nx, nz, grid, region)
    REAL(KIND=dp), DIMENSION(-nz:nz-1,0:nx-1,1:4), INTENT(in) :: grid
    INTEGER,       DIMENSION(    -nz:nz-1,0:nx-1), INTENT(in) :: region    
    INTEGER :: kk, nx, nz
    
    INTEGER :: jj, ii
    
    DO jj = -nz+1, nz
    DO ii = 1, nx 
      p_dens(  kk,ii,jj) = grid(jj-1,ii-1,1)
      p_ux(    kk,ii,jj) = grid(jj-1,ii-1,3)
      p_uz(    kk,ii,jj) = grid(jj-1,ii-1,2)
      p_temp(  kk,ii,jj) = grid(jj-1,ii-1,4)
      p_region(kk,ii,jj) = region(jj-1,ii-1)
    END DO
    END DO
    
!    DO jj = -nx+1, nx        
!    DO ii = 1, nx 
!      p_dens(  kk,ii,jj) =  grid(1,-jj,ii)
!      p_ux(    kk,ii,jj) =  grid(2,-jj,ii)
!      p_uz(    kk,ii,jj) = -grid(3,-jj,ii)
!      p_temp(  kk,ii,jj) =  grid(4,-jj,ii)
!      p_region(kk,ii,jj) =  region(-jj,ii)
!    END DO
!    END DO
        
  END SUBROUTINE pass_plasma_grid_mhdam
  
  SUBROUTINE pass_plasma_grid(kk, nx, nz, grid, region)
    REAL(KIND=dp), DIMENSION(:,:,:), INTENT(in) :: grid
    INTEGER,       DIMENSION(1:nx,-nz+1:nz), INTENT(in) :: region
    INTEGER :: kk, nx, nz

    p_dens(kk,1:nx,-nz+1:nz) = grid(1,:,:)
    p_ux(  kk,1:nx,-nz+1:nz) = grid(2,:,:)
    p_uz(  kk,1:nx,-nz+1:nz) = grid(3,:,:)
    p_temp(kk,1:nx,-nz+1:nz) = grid(4,:,:)
    p_region(kk,1:nx,-nz+1:nz) = region(:,:)

  END SUBROUTINE pass_plasma_grid
!
!***********************************************************************
!
  SUBROUTINE return_source_grid_mhdam(kk, nx, nz, grid)
    REAL(KIND=dp), DIMENSION(-nz:nz-1,0:nx-1,1:5), INTENT(out) :: grid
    INTEGER :: kk, nx, nz
    INTEGER :: ii,jj
    
    IF (my_rank .NE. 0) THEN  ! only master has complete neutral grid data
      return
    END IF

    DO jj = -nz+1, nz        
    DO ii = 1, nx 
     grid(jj-1,ii-1,1) =  g_source_m(    kk,ii,jj)
     grid(jj-1,ii-1,3) =  g_source_px(   kk,ii,jj)
     grid(jj-1,ii-1,2) =  g_source_pz(   kk,ii,jj)
     grid(jj-1,ii-1,4) =  g_source_mvsq( kk,ii,jj) * 0.5_dp
     grid(jj-1,ii-1,5) =  g_source_nchex(kk,ii,jj)
    END DO
    END DO

    
    !~DO jj = -nx+1, nx        
    !~DO ii = 1, nx 
    !~ grid(1,-jj,ii) =  g_source_m(    kk,ii,jj)
    !~ grid(2,-jj,ii) =  g_source_px(   kk,ii,jj)
    !~ grid(3,-jj,ii) = -g_source_pz(   kk,ii,jj)
    !~ grid(4,-jj,ii) =  g_source_mvsq( kk,ii,jj) * 0.5_dp
    !~END DO
    !~END DO
    

  END SUBROUTINE return_source_grid_mhdam
!
!***********************************************************************
!
  SUBROUTINE return_source_grid(kk, nx, nz, grid, nchex)
    REAL(KIND=dp), DIMENSION(:,:,:), INTENT(out) :: grid
    INTEGER,       DIMENSION(:,:), INTENT(out) :: nchex
    INTEGER :: kk, nx, nz
    
    IF (my_rank .EQ. 0) THEN      ! only master has complete sources
       grid(1,:,:) = g_source_m(    kk,1:nx,-nz+1:nz)
       grid(2,:,:) = g_source_px(   kk,1:nx,-nz+1:nz)
       grid(3,:,:) = g_source_pz(   kk,1:nx,-nz+1:nz)
       grid(4,:,:) = g_source_mvsq( kk,1:nx,-nz+1:nz) / 2.0_dp
       nchex(:,:)  = g_source_nchex(kk,1:nx,-nz+1:nz)
    END IF

  END SUBROUTINE return_source_grid

!
!***********************************************************************
!
  SUBROUTINE return_neutral_grid_mhdam(kk, nx, nz, region, grid)
    REAL(KIND=dp), DIMENSION(-nz:nz-1,0:nx-1,1:4), INTENT(out) :: grid
    INTEGER :: kk, nx, nz, region
    
    INTEGER :: ii,jj
    
    IF (my_rank .NE. 0) THEN  ! only master has complete neutral grid data
      return
    END IF
      
    DO jj = -nz+1, nz        
    DO ii = 1, nx 
       grid(jj-1,ii-1,1) = grid_dens( kk,region,ii,jj)
       grid(jj-1,ii-1,3) = grid_ux(   kk,region,ii,jj)
       grid(jj-1,ii-1,2) = grid_uz(   kk,region,ii,jj)
       grid(jj-1,ii-1,4) = grid_temp( kk,region,ii,jj)       
    END DO
    END DO
      
!    DO jj = -nx+1, nx        
!    DO ii = 1, nx 
!       grid(1,-jj,ii) =  grid_dens( kk,region,ii,jj)
!       grid(2,-jj,ii) =  grid_ux(   kk,region,ii,jj)
!       grid(3,-jj,ii) = -grid_uz(   kk,region,ii,jj)
!       grid(4,-jj,ii) =  grid_temp( kk,region,ii,jj)       
!    END DO
!    END DO

  END SUBROUTINE return_neutral_grid_mhdam

!
!***********************************************************************
!
  SUBROUTINE return_neutral_grid(kk, nx, nz, region, grid, npart)
    REAL(KIND=dp), DIMENSION(:,:,:), INTENT(out) :: grid
    INTEGER,       DIMENSION(:,:), INTENT(out) :: npart
    INTEGER :: kk, nx, nz, region    
            
    IF (my_rank .EQ. 0) THEN      ! only master has complete neutral grid data
       grid(1,:,:) = grid_dens( kk,region,1:nx,-nz+1:nz)
       grid(2,:,:) = grid_ux(   kk,region,1:nx,-nz+1:nz)
       grid(3,:,:) = grid_uz(   kk,region,1:nx,-nz+1:nz)
       grid(4,:,:) = grid_temp( kk,region,1:nx,-nz+1:nz)

       npart(:,:)  = grid_npart(kk,region,1:nx,-nz+1:nz)
    END IF

  END SUBROUTINE return_neutral_grid
!
!***********************************************************************
!
  SUBROUTINE write_plasma
    INTEGER :: kk, nx, nz

    IF (my_rank .EQ. 0) THEN
     OPEN(unit=1,file='plasma.data',form='unformatted',status='unknown',action='write')
     WRITE(1) ngrids
     WRITE(1) grids_nxmax
     WRITE(1) grids_nzmax
     DO kk = ngrids, 1, -1
        nx = grids_nx(kk)
        nz = grids_nz(kk)
        WRITE(1) nx
        WRITE(1) nz
        WRITE(1) grids_boundaries_x(kk)
        WRITE(1) grids_boundaries_z(kk)
        WRITE(1) p_dens(kk,1:nx,-nz+1:nz)
        WRITE(1) p_ux(kk,1:nx,-nz+1:nz)
        WRITE(1) p_uz(kk,1:nx,-nz+1:nz)
        WRITE(1) p_temp(kk,1:nx,-nz+1:nz)
        WRITE(1) p_region(kk,1:nx,-nz+1:nz)
     END DO
     CLOSE(1)

  END IF


  END SUBROUTINE write_plasma
  
  !
!***********************************************************************
!
  SUBROUTINE write_plasma_ascii
    REAL(KIND=dp) :: dx, xx, zz
    INTEGER :: kk, ii, jj, nx, nz

    IF (my_rank .EQ. 0) THEN


       OPEN(unit=1,file='plasma_ascii.data',form='formatted',status='unknown',action='write')
       WRITE(1,*) ngrids
       WRITE(1,*) grids_nxmax
       WRITE(1,*) grids_nzmax
       DO kk = ngrids, 1, -1
          dx = grids_dx(kk)/au
          nx = grids_nx(kk)
          nz = grids_nz(kk)
          WRITE(1,*) nx
          WRITE(1,*) nz
          WRITE(1,*) grids_boundaries_x(kk)
          WRITE(1,*) grids_boundaries_z(kk)
          DO ii = 1, nx
             xx = (ii-0.5_dp)*dx
             DO jj = -nz+1, nz
                zz = (jj-0.5_dp)*dx

                WRITE(1,19) xx, zz, p_dens(kk,ii,jj), p_ux(kk,ii,jj), &
                     p_uz(kk,ii,jj), p_temp(kk,ii,jj),p_region(kk,ii,jj)

             END DO
          END DO
       END DO
       CLOSE(1)

    END IF

19  FORMAT(6F14.4,i6)

  END SUBROUTINE write_plasma_ascii

!
!***********************************************************************
!
  SUBROUTINE write_source
    INTEGER :: kk, nx, nz

    IF (my_rank .EQ. 0) THEN

       OPEN(unit=1,file='source.data',form='unformatted',status='unknown',action='write')
       WRITE(1) ngrids
       WRITE(1) grids_nxmax
       WRITE(1) grids_nzmax
       DO kk = ngrids, 1, -1
          nx = grids_nx(kk)
          nz = grids_nz(kk)
          WRITE(1) nx
          WRITE(1) nz
          WRITE(1) grids_boundaries_x(kk)
          WRITE(1) grids_boundaries_z(kk)
          WRITE(1) g_source_m(kk,1:nx,-nz+1:nz)
          WRITE(1) g_source_px(kk,1:nx,-nz+1:nz)
          WRITE(1) g_source_pz(kk,1:nx,-nz+1:nz)
          WRITE(1) g_source_mvsq(kk,1:nx,-nz+1:nz) / 2.0_dp
          WRITE(1) g_source_nchex(kk,1:nx,-nz+1:nz)
       END DO
       CLOSE(1)

    END IF

  END SUBROUTINE write_source
!
!***********************************************************************
!
  SUBROUTINE write_neutral_grid
    INTEGER :: kk, nx, nz, region

    IF (my_rank .EQ. 0) THEN

       OPEN(unit=1,file='neutral_grid.data',form='unformatted',status='unknown',action='write')
       WRITE(1) ngrids
       WRITE(1) grids_nxmax
       WRITE(1) grids_nzmax
       DO kk = ngrids, 1, -1
          nx = grids_nx(kk)
          nz = grids_nz(kk)
          WRITE(1) nx
          WRITE(1) nz
          WRITE(1) grids_boundaries_x(kk)
          WRITE(1) grids_boundaries_z(kk)
          
          DO region = -1, 3
             WRITE(1) grid_dens( kk,region,1:nx,-nz+1:nz)
             WRITE(1) grid_ux(   kk,region,1:nx,-nz+1:nz)
             WRITE(1) grid_uz(   kk,region,1:nx,-nz+1:nz)
             WRITE(1) grid_temp( kk,region,1:nx,-nz+1:nz)
             WRITE(1) grid_npart(kk,region,1:nx,-nz+1:nz)
          END DO

       END DO
       CLOSE(1)

    END IF

   END SUBROUTINE write_neutral_grid
!
!***********************************************************************
!
  SUBROUTINE mc_init(p_xmax_in, p_zmin_in, p_zmax_in, &
                     ngrids_in, grids_nx_in, grids_nz_in, &
                     grids_boundaries_x_in, grids_boundaries_z_in, &
                     restart_in, load_file_in, au_total_neutrals, &
                     lism_nh_in, lism_vh_in, lism_th_in, &
                     verbosity_in, photoionize_in)
    use hdf5
    REAL(KIND=dp) :: lism_nh_in, lism_vh_in, lism_th_in
    REAL(KIND=dp) :: p_xmax_in, p_zmin_in, p_zmax_in, achieved_density
    REAL(KIND=dp) :: xk, xkm1, zk, zkm1, total_volume2
    REAL(KIND=dp), DIMENSION(:) :: grids_boundaries_x_in, grids_boundaries_z_in
    INTEGER, DIMENSION(:) :: grids_nx_in, grids_nz_in
    INTEGER, DIMENSION(1:4) :: put_seq
    INTEGER, intent(in) :: ngrids_in, restart_in, verbosity_in, photoionize_in
    INTEGER :: kk, n, ierr
    REAL(KIND=dp) :: au_total_neutrals,total_neutrals
    CHARACTER (LEN = lenstr) :: load_file_in
    INTEGER :: ncores,i    
    REAL(8)    :: aux_num    
    CHARACTER (LEN = 100) str_buf
    INTEGER :: status(MPI_STATUS_SIZE), loc_nprocMPI, loc_my_rank
    
    INTEGER(8) :: i64
    INTEGER(HSIZE_T) iHSIZE_T
    integer size1,size2
    
#ifdef CH_OMPCPP               
    INTEGER OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, loc_my_rank, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, loc_nprocMPI, ierr)    ! nprocMPI = # of MPI processors
    
    allocate(mpi_balancing(0:loc_nprocMPI-1))
    
    call h5open_f(ierr)
    
    inquire(iolength=size1)i64
    inquire(iolength=size2)iHSIZE_T
    
    if ((size1.ne.size2).or.(size1.ne.8).or.(size2.ne.8)) then    
      print *,'mc_init Error: INTEGER(8) and HSIZE_T types must be both 64-bit, whereas their lengths are ', size1,' and ', & 
            size2,' bytes respectively'
    !  call mpi_abort()
    endif  
    
        
    nthreads  = 1
    my_thread = 0    
    my_rank   = loc_my_rank
    nprocMPI  = loc_nprocMPI
    
    
    !print *,'Entering omp region'
    !$omp parallel private(put_seq,n,aux_num,i,str_buf)
#ifdef CH_OMPCPP               
      my_thread = OMP_GET_THREAD_NUM()          
      nthreads  = omp_get_num_threads()         
#endif 
      
      my_rank  = loc_my_rank
      nprocMPI = loc_nprocMPI
    
      call init_pout()
  
      
      
      write (str_buf,*) 'Threads allocated : ',nthreads
      call pout(str_buf)      
                  
                
      ! Set random number seed to be different for each processor      
      put_seq = 10*my_thread + my_rank * nprocMPI + 37 * (/ (i - 1, i = 1, 4) /)
      put_seq = my_rank*nthreads + my_thread   + 37 * (/ (i - 1, i = 1, 4) /)
      put_seq = my_rank + my_thread*loc_nprocMPI  + 37 * (/ (i - 1, i = 1, 4) /)
      !put_seq(1) = 50
      !put_seq(2) = put_seq(1) + 10
      !put_seq(3) = put_seq(2) + 10
      !put_seq(4) = put_seq(3) + 10
      CALL rand_seed(put_seq)
      
      CALL rand_number(aux_num)
      my_nchex_events  = 80000 + NINT(30000*aux_num)                   
      my_nchex_current = 0
!      !$omp critical (allocate_my_chex_events)
      allocate(my_chex_events(1:my_nchex_events))     
 !     !$omp end critical (allocate_my_chex_events)

                              
                  
    !$omp end parallel
    
    !print *,'Leaving omp region'
    !print *,'nthreads = ',nthreads
     
    ncores   = nprocMPI*nthreads
    
    p_xmax = p_xmax_in
    p_zmin = p_zmin_in
    p_zmax = p_zmax_in
    p_xmax_sq = p_xmax*p_xmax

    ngrids = ngrids_in
    
    
    call init_vr_phi 
    call init_Maxwellian_vr
        
    
    !print *,'ALLOCATE(grids_nx(ngrids), grids_nz(ngrids), grids_dx(ngrids))'
    
    ALLOCATE(grids_nx(ngrids), grids_nz(ngrids), grids_dx(ngrids))            
    ALLOCATE(grids_boundaries_x(0:ngrids), grids_boundaries_z(0:ngrids))
    grids_boundaries_x(0) = 0.0_dp;  grids_boundaries_z(0) = 0.0_dp

    grids_nx = grids_nx_in
    grids_nz = grids_nz_in
    grids_nxmax = MAXVAL(grids_nx)
    grids_nzmax = MAXVAL(grids_nz)
    grids_boundaries_x = grids_boundaries_x_in
    grids_boundaries_z = grids_boundaries_z_in
    DO kk = 1, ngrids
! ASSUME SQUARE grid cells
       grids_dx(kk) = grids_boundaries_x(kk)/grids_nx(kk)
    END DO
    total_ndata = ngrids*2*grids_nxmax*grids_nzmax

    IF (restart_in .EQ. 0) THEN
       restart = .FALSE.
    ELSE
       restart = .TRUE.
    END IF
    IF (photoionize_in .EQ. 0) THEN
       photoionize = .FALSE.
    ELSE
       photoionize = .TRUE.
    END IF      
    
    load_file = load_file_in
    verbosity = verbosity_in
    
    ! set up global parameters

    ! LISM parameters
    LISM_nH = lism_nh_in               ! per cubic meter
    LISM_vH = lism_vh_in               ! meters per second
    LISM_TH = lism_th_in               ! degrees Kelvin

    k_B = 1.38e-23_dp                  ! Joules per Kelvin (m^2 kg s^{-2} K^{-1})
    neutral_mass = 1.67e-27_dp         ! Kilograms
    LISM_v_thermal = SQRT(2.0_dp*k_B*LISM_TH/neutral_mass) ! 3D thermal speed

    !print *,'allocating source grids '

    ALLOCATE(source_px(ngrids,   grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(source_pz(ngrids,   grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(source_m(ngrids,    grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(source_mvsq(ngrids, grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(source_nchex(ngrids,grids_nxmax,-grids_nzmax+1:grids_nzmax))
    
    !print *,'allocating global source grids '

! output cartesian source grids
    ALLOCATE(g_source_px(ngrids,   grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(g_source_pz(ngrids,   grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(g_source_m(ngrids,    grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(g_source_mvsq(ngrids, grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(g_source_nchex(ngrids,grids_nxmax,-grids_nzmax+1:grids_nzmax))
    
    !print *,'allocating plasma grid '

! allocate plasma grid
    ALLOCATE(  p_dens(ngrids,  grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(  p_temp(ngrids,  grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(  p_ux(ngrids,    grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(  p_uz(ngrids,    grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(  p_region(ngrids,grids_nxmax,-grids_nzmax+1:grids_nzmax))
    p_dens = 0.0_dp;  p_temp = 0.0_dp
    p_ux = 0.0_dp;    p_uz = 0.0_dp
    p_region = 0
    
    !print *,'allocating master plasma grid '

! allocate neutral grid (ONLY on MASTER) -- use regions (-1 is total)
    IF (my_rank .EQ. 0) THEN
       ALLOCATE(  grid_dens(ngrids,  -1:3,grids_nxmax,-grids_nzmax+1:grids_nzmax))
       ALLOCATE(  grid_temp(ngrids,  -1:3,grids_nxmax,-grids_nzmax+1:grids_nzmax))
       ALLOCATE(  grid_ux(ngrids,    -1:3,grids_nxmax,-grids_nzmax+1:grids_nzmax))
       ALLOCATE(  grid_uz(ngrids,    -1:3,grids_nxmax,-grids_nzmax+1:grids_nzmax))
       ALLOCATE(  grid_npart(ngrids, -1:3,grids_nxmax,-grids_nzmax+1:grids_nzmax))
       grid_dens = 0.0_dp;  grid_temp = 0.0_dp
       grid_ux = 0.0_dp;    grid_uz = 0.0_dp
       grid_npart = 0
    END IF
    
    
    !print *,'all shared grids allocated '

! photoionization region outer boundary
    phi_rmax_sq = (400.0_dp*au)**2

    timescale = (p_zmax-p_zmin)/ABS(LISM_vH)
    time = 0.0_dp
    

! set up coarse timestep which moves LISM atoms about 1.5% of a mfp
!    timestep = secs_per_year / 2.0_dp    ! half a year
    timestep = secs_per_year             ! one year


    total_volume = (p_xmax+boundary_layer)*((p_zmax-p_zmin)+2.0_dp*boundary_layer)
! ASSUME uniformly dens GRID-CELLS
    total_volume2 = 0.0_dp
    DO kk = 1, ngrids
       xk = grids_boundaries_x(kk)
       xkm1 = grids_boundaries_x(kk-1)
       zk = grids_boundaries_z(kk)
       zkm1 = grids_boundaries_z(kk-1)
       total_volume2 = total_volume2 + 2.0_dp*4.0_dp**(ngrids-kk)*(xk*zk-xkm1*zkm1)
    END DO

  total_neutrals = 3*au_total_neutrals*(total_volume/total_volume2) / ngrids / 2

! define this in terms of total_neutrals first, then update below
    LISM_density = total_neutrals/total_volume

! choose a boundary ring thickness so that 99.9999998% of particles with
! velocity directed inward which make it into the domain come from this layer
    boundary_layer = 6.0_dp*(LISM_v_thermal+ABS(LISM_vH))*timestep

! number of particles from 0 to p_xmax+boundary_layer
    my_npart_boundary = NINT(LISM_density*boundary_layer* &
                             (p_xmax+boundary_layer)/ncores)
    achieved_density = ncores*my_npart_boundary &
                      / (boundary_layer*(p_xmax+boundary_layer))
    npart = NINT(achieved_density*total_volume,8)

! NOTE: must define macro_mass in terms of EXPECTED # of particles (npart)
    mega_number = LISM_nH*total_volume/npart ! # neutrals in megaparticle
    macro_mass = neutral_mass * mega_number           ! mass of macroparticle
    LISM_density = npart/total_volume
       
    PRINT *,'expected total number of particles'
    PRINT *, npart 

    IF ( (my_rank .eq. 0) .and. (verbosity .ge. 2) ) THEN
       PRINT *,'lengthscale (m)'
       PRINT 9, p_xmax
       PRINT *,'timescale (sec)'
       PRINT 9, timescale
       PRINT *,'timestep (sec)'
       PRINT 9, timestep
       PRINT *,'# of neutrals per megaparticle'
       PRINT 9, mega_number
       PRINT *,'size of boundary layer/xmax'
       PRINT 9, boundary_layer/p_xmax
       PRINT *,'LISM number density'
       PRINT 9, LISM_density
       PRINT *,'achieved number density'
       PRINT 9, achieved_density
       PRINT *,'number of particle per processor per boundary sector'
       PRINT *, my_npart_boundary
       PRINT *,'expected total number of particles'
       PRINT *, npart
    END IF

    my_ntotal = npart/ncores + ncores       ! allow for rounding
! allow for SPLITTING with no recombination for comp1 and comp2 neutrals
    my_ntotal = 5*(ngrids*my_ntotal + 400*NINT(SQRT(REAL(my_ntotal,dp))))
    
    

    IF ( (my_rank .eq. 0) .and. (verbosity .ge. 1) ) THEN
       PRINT *,'maximum number of particles'
       PRINT *, ncores*my_ntotal
    END IF
    
    PRINT *,'maximum number of particles'
    PRINT *, int(ncores,8)*int(my_ntotal,8)
    


! mass, x, y, z, vx, vy, vz, region

    !print *,'allocate  my_particles'
    
    allocate(my_particles(0:nthreads-1))     
    do i=0, nthreads-1
      my_particles(i)%my_npart = 0
      ALLOCATE(my_particles(i)%particles(my_ntotal,0:7))   
      my_particles(i)%particles = 0.0_dp
    enddo
    
    !print *,'my_particles allocated'
    
  !~!$omp parallel
  !~    allocate(my_particles(0:nthreads-1))        
  !~    my_particles(my_thread)%my_npart = 0
  !~    ALLOCATE(my_particles(my_thread)%particles(my_ntotal,0:7))   
  !~    my_particles(my_thread)%particles = 0.0_dp
  !~!$omp end parallel
    

    time = 0.0_dp
    IF ( restart ) THEN
! load particle data from existing file
       CALL load_raw
       CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
    ELSE
! set up initial distribution 
!       CALL inject_maxwellian(0.0_dp,p_xmax,-p_xmax,p_xmax, &
!                              LISM_v_thermal,npart/nproc)
    END IF
    IF ( (my_rank .eq. 0) .and. (verbosity .ge. 1) ) THEN
       PRINT *,'number of timesteps'
       PRINT 9, end_time/timestep
       IF ( restart ) THEN
          PRINT *, 'time in loaded calculation'
          PRINT 9, time
       END IF
    END IF

! run for another 'endtime'
    end_time = end_time + time
    sourcetime = sourcetime + time

9 FORMAT(6e13.4)

  END SUBROUTINE mc_init
!
!***********************************************************************
!
  SUBROUTINE move
    REAL(KIND=dp), PARAMETER :: ideal_chex_prob = 0.02_dp
    REAL(KIND=dp) :: massi, xi, yi, zi, ri, vx, vz, vy, region, xyr_sq, ri_sq
    REAL(KIND=dp) :: n_p, T_p, ux, uz, region_p
    REAL(KIND=dp) :: n_p2, T_p2, ux2, uz2, region_p2
    REAL(KIND=dp) :: xn, yn, zn, vxn, vzn, vyn
    REAL(KIND=dp) :: xr, vxr, vxy_sq, vyr_sq, v_p_th, delta_u, omega_p
    REAL(KIND=dp) :: xr2, vxr2, vyr_sq2, v_p_th2, delta_u2, omega_p2
    REAL(KIND=dp) :: v_rel_p, sigma_ex, subtime, dt, ideal_dt, ds
    REAL(KIND=dp) :: v_rel_p2, ideal_dt2
    REAL(KIND=dp) :: beta_ph, mass_loss, dprand
    REAL(KIND=dp) :: derf  ! need to link in error function
    INTEGER :: my_old_npart, i, count, npart_inject
    
    my_min_dt = timestep
    my_max_dt = 0.0_dp
    my_ave_dt = 0.0_dp
    count = 0
    my_n_chex = 0

    my_tot_mass_loss = 0.0_dp
    my_old_npart = my_particles(my_thread)%my_npart
    i = 0
    DO WHILE (i .LT. my_particles(my_thread)%my_npart)
       i = i + 1

       massi = my_particles(my_thread)%particles(i,0)   
       xi  = my_particles(my_thread)%particles(i,1)
       yi  = my_particles(my_thread)%particles(i,2)
       zi  = my_particles(my_thread)%particles(i,3)
       vx = my_particles(my_thread)%particles(i,4)
       vy = my_particles(my_thread)%particles(i,5)
       vz = my_particles(my_thread)%particles(i,6)
       region = my_particles(my_thread)%particles(i,7)

       ri = sqrt(xi*xi + yi*yi + zi*zi)

       subtime = 0.0_dp
! use adaptive timestep
       DO WHILE (subtime .LT. timestep)

          xr = SQRT(xi*xi + yi*yi)
          CALL plasma(xr,zi,n_p,T_p,ux,uz,region_p)
          v_p_th = SQRT(2.0_dp*k_B*T_p/neutral_mass)

          vxr = (vx*xi + vy*yi)/xr
          vxy_sq = vx*vx + vy*vy
          vyr_sq = ABS(vxy_sq - vxr*vxr)

! from Ripken and Fahr 1983, notation follows Lipatov, Zank & Pauls 1998
          delta_u = SQRT( (vxr-ux)**2 + (vz-uz)**2 + vyr_sq)
          omega_p = delta_u/v_p_th
          v_rel_p = v_p_th * ( EXP(-omega_p*omega_p)/SQRT(pi) + &
                                   (omega_p + 0.5_dp/omega_p)*derf(omega_p) )
   sigma_ex = ( 1.64e-7_dp - 6.95e-9_dp*LOG(v_rel_p/1e-2_dp) )**2 * 1e-4_dp !m^2
   
          ideal_dt = ideal_chex_prob/(n_p*sigma_ex*v_rel_p)
! NOTE: might need to check that dt is small enough over entire space step
! (seems ok, even at 1 AU)

          dt = MIN(ideal_dt,timestep)

          IF ( xr .LT. 4.0_dp*au ) THEN
! improve accuracy near axis
             dt = dt * 0.5_dp
             IF ( xr .LT. 2.0_dp*au ) THEN
                dt = dt * 0.5_dp
                IF ( xr .LT. 1.0_dp*au ) THEN
                   dt = dt * 0.5_dp
                END IF
             END IF
          END IF

! check FUTURE dt
          xr2 = SQRT((xi+vx*dt)**2+(yi+vy*dt)**2)
          CALL plasma(xr2,zi+vz*dt,n_p2,T_p2,ux2,uz2,region_p2)
          v_p_th2 = SQRT(2.0_dp*k_B*T_p2/neutral_mass)

          vxr2 = (vx*xi + vy*yi)/xr2
          vyr_sq2 = ABS(vxy_sq - vxr2*vxr2)
          delta_u2 = SQRT( (vxr2-ux2)**2 + (vz-uz2)**2 + vyr_sq2)
          omega_p2 = delta_u2/v_p_th2
          v_rel_p2 = v_p_th2 * ( EXP(-omega_p2*omega_p2)/SQRT(pi) + &
                                   (omega_p2 + 0.5_dp/omega_p2)*derf(omega_p2) )
          ideal_dt2 = ideal_chex_prob/(n_p2*sigma_ex*v_rel_p2)

          IF ( (dt/ideal_dt2 .GT. 2.0_dp) .and. (ri .gt. 10.0_dp*au) ) THEN
! IMPORTANT: reduce dt to resolve discontinuities in the source terms
             dt = 0.1_dp*min(ideal_dt2,timestep)
          END IF

! DIAGNOSTIC:
          my_max_dt = MAX(my_max_dt,dt)
          my_min_dt = MIN(my_min_dt,dt)

          IF (subtime+dt .gt. timestep) THEN
             dt = timestep - subtime
          END IF

          my_ave_dt = my_ave_dt + dt
          count = count + 1

          ds = v_rel_p * dt  ! use (RELATIVE) distance travelled for Prob(chex)
          CALL rand_number(dprand)
          IF ( dprand .LT. n_p*sigma_ex*ds ) THEN   ! charge exchange occurs

             CALL exchange(massi,xi,yi,zi,vx,vy,vz,region, &
                           ux,uz,v_p_th,xn,yn,zn,vxn,vyn,vzn)

             xi = xn;  yi = yn;  zi = zn
             vx = vxn;  vy = vyn;  vz = vzn

             region = region_p

          END IF

          if (photoionize) then
             ri_sq = xi*xi + yi*yi + zi*zi
! NOTE: need to photoionize each small step for accuracy close to the Sun
! photoionization; assumes beta_ph = beta_ph_E at 1au
             beta_ph = beta_ph_E * r_E_sq / ri_sq
             mass_loss = beta_ph * dt * massi
             my_tot_mass_loss = my_tot_mass_loss + mass_loss
             CALL source_phi(xi,yi,zi,vx,vy,vz,mass_loss)
             massi = massi - mass_loss
          end if

          subtime = subtime + dt
! advance particles in a straight line at current velocity (Euler's method)
          xi = xi + vx*dt
          yi = yi + vy*dt
          zi = zi + vz*dt

       END DO

! update particle specs
       my_particles(my_thread)%particles(i,0) = massi
       my_particles(my_thread)%particles(i,1) = xi
       my_particles(my_thread)%particles(i,2) = yi
       my_particles(my_thread)%particles(i,3) = zi
       my_particles(my_thread)%particles(i,4) = vx
       my_particles(my_thread)%particles(i,5) = vy
       my_particles(my_thread)%particles(i,6) = vz
       my_particles(my_thread)%particles(i,7) = region

       xyr_sq = xi*xi + yi*yi

! remove particles that have moved out of the domain
! photoionization will take care of rmin bound
       IF ( ( xyr_sq .GT. p_xmax_sq) .OR. (zi .GT. p_zmax) .OR. (zi .lT. p_zmin) &
            .OR. (massi .LT. 0.0_dp) )  THEN
          my_particles(my_thread)%particles(i,:) = my_particles(my_thread)%particles(my_particles(my_thread)%my_npart,:)   ! fill gap with last particle
          my_particles(my_thread)%my_npart = my_particles(my_thread)%my_npart - 1    ! stop tracking particles outside domain
          i = i - 1                  ! don't increment if we remove a particle          
       END IF

    END DO

    my_ave_dt = my_ave_dt / count

    my_particles(my_thread)%particles(my_particles(my_thread)%my_npart+1:my_old_npart,:) = 0.0_dp ! zero removed particles
    
    !write(pout_unit,*) 'remove particles (proc,thread,#particles)',my_rank,my_thread,my_old_npart-my_particles(my_thread)%my_npart
    
    !write(pout_unit,*) 'thread: ',my_thread
    !write(pout_unit,*) '# particles ',my_particles(my_thread)%my_npart
    
    !do i=1,my_particles(my_thread)%my_npart
    !  write(pout_unit,*) i,my_particles(my_thread)%particles(i,:)
    !enddo
        
! inject at UPSTREAM boundary
    IF (my_particles(my_thread)%my_npart+my_npart_boundary .LE. my_ntotal) THEN
       CALL inject_maxwellian(0.0_dp,p_xmax+boundary_layer, &
            p_zmax,p_zmax+boundary_layer, LISM_v_thermal,my_npart_boundary)
            
    !   write(pout_unit,*)'inject at UPSTREAM boundary',my_rank,my_thread,my_npart_boundary
       
    END IF

! inject at SIDESTREAM boundary
! NOTE here we only inject particles into the PLASMA domain (p_zmin ... p_zmax)
    npart_inject = my_npart_boundary * nint((p_zmax-p_zmin)/(p_xmax+boundary_layer))
    IF (my_particles(my_thread)%my_npart+npart_inject .LE. my_ntotal) THEN
       CALL inject_maxwellian(p_xmax,p_xmax+boundary_layer, &
                              p_zmin,p_zmax,LISM_v_thermal,npart_inject)
!       write(pout_unit,*)'inject at SIDESTREAM boundary',my_rank,my_thread,npart_inject
    END IF

! inject at DOWNSTREAM boundary
    IF (my_particles(my_thread)%my_npart+my_npart_boundary .LE. my_ntotal) THEN
       CALL inject_maxwellian(0.0_dp,p_xmax+boundary_layer, &
            p_zmin-boundary_layer,p_zmin, LISM_v_thermal,my_npart_boundary)
!       write(pout_unit,*)'inject at DOWNSTREAM boundary',my_rank,my_thread,my_npart_boundary
    END IF
    
!    do i=1,my_particles(my_thread)%my_npart
!      write(pout_unit,*) i,my_particles(my_thread)%particles(i,:)
!    enddo

  END SUBROUTINE move
!
!***********************************************************************
!
  SUBROUTINE split
    REAL(KIND=dp) :: massi, xi, yi, zi, vxi, vyi, vzi, region
    REAL(KIND=dp) :: n_p, T_p, ux, uz, region_p
    REAL(KIND=dp) :: new_massi, ideal_massi
    REAL(KIND=dp) :: dprand, xyr, dxi, ri, xg
    INTEGER :: my_old_npart, part, ii, nn, kkx, kkz, kk, ind

    my_old_npart = my_particles(my_thread)%my_npart
    part = 0
    DO WHILE (part .LT. my_particles(my_thread)%my_npart)
       part = part + 1

       xi = my_particles(my_thread)%particles(part,1)
       yi = my_particles(my_thread)%particles(part,2)
       zi = my_particles(my_thread)%particles(part,3)

       region = my_particles(my_thread)%particles(part,7)
       massi = my_particles(my_thread)%particles(part,0)

       xyr = SQRT(xi*xi + yi*yi)
       ri = SQRT(xyr*xyr + zi*zi)

       IF (xyr .GT. 8.0_dp*au) THEN
          xyr = xyr - 8.0_dp*au     ! split before crossing grid boundary
       END IF
       IF (ri .GT. 8.0_dp*au) THEN
          ri = ri - 8.0_dp*au     ! split before crossing grid boundary
       END IF

! find which grid we're on and nearest grid point
       kkx = -1
       DO nn = 1, ngrids
          IF (abs(xyr) .GE. grids_boundaries_x(ngrids-nn)) THEN
             kkx = ngrids-nn+1
             EXIT
          END IF
       END DO
       kkz = -1
       DO nn = 1, ngrids
          IF (abs(zi) .GE. grids_boundaries_z(ngrids-nn)) THEN
             kkz = ngrids-nn+1
             EXIT
          END IF
       END DO
       kk = max(kkx,kkz)
       dxi = grids_dx(kk)
       ii = floor(xyr / dxi) + 1

! split as though we are 1/4 up grid-cell from the x-axis
       xg = (ii-0.75_dp)*dxi

       ideal_massi = xg*twopi / 4.0_dp**(ngrids-kk)

! extra FIDDLING
       IF ( xyr .GT. p_xmax/8.0_dp ) THEN        ! should exclude termination shock
          ideal_massi = 2.0_dp * ideal_massi
          IF ( xyr .GT. p_xmax/4.0_dp ) THEN     ! should exclude termination shock
             ideal_massi = 2.2_dp * ideal_massi
             IF ( xyr .GT. p_xmax/2.0_dp ) THEN  ! should exclude heliosheath
                ideal_massi = 2.0_dp * ideal_massi
             END IF
          END IF
       END IF

! allow for both positive and negative LISM FLOW
       IF ( ( (xyr .LT. 12.0_dp*dxi) .AND. (zi .LT.  150.0_dp*au) &
            .and. (sign(1.0_dp,LISM_vh) .lt. 0.0_dp) ) .or. &
            ( (xyr .LT. 12.0_dp*dxi) .AND. (zi .gT. -150.0_dp*au) &
            .and. (sign(1.0_dp,LISM_vh) .gt. 0.0_dp) ) ) THEN
! try to improve axis
          ideal_massi = ideal_massi / 2.2_dp
          IF ( (xyr .LT. 6.0_dp*dxi) ) THEN
             ideal_massi = ideal_massi / 2.0_dp
             IF ( (xyr .LT. 3.0_dp*dxi) ) THEN
                ideal_massi = ideal_massi / 2.0_dp
             END IF
          END IF
       END IF

       IF (massi .GT. 2.0_dp*ideal_massi) THEN

          IF ( (my_particles(my_thread)%my_npart .GT. my_ntotal-1) .and. (verbosity .ge. 1) ) THEN
             PRINT *, 'cant split'
          ELSE

             new_massi = massi / 2.0_dp
             my_particles(my_thread)%particles(part,0) = new_massi

             vxi = my_particles(my_thread)%particles(part,4)
             vyi = my_particles(my_thread)%particles(part,5)
             vzi = my_particles(my_thread)%particles(part,6)
             
             ind = my_particles(my_thread)%my_npart + 1

             my_particles(my_thread)%particles(ind,0) = new_massi
             my_particles(my_thread)%particles(ind,1) = xi
             my_particles(my_thread)%particles(ind,2) = yi
             my_particles(my_thread)%particles(ind,3) = zi
             my_particles(my_thread)%particles(ind,4) = vxi
             my_particles(my_thread)%particles(ind,5) = vyi
             my_particles(my_thread)%particles(ind,6) = vzi
             my_particles(my_thread)%particles(ind,7) = region

             my_particles(my_thread)%my_npart = my_particles(my_thread)%my_npart + 1

          END IF
       END IF

       IF ( ( (zi .LT.  120.0_dp*au) .AND. (xyr .LT. p_xmax/2.0_dp) &
            .and. (sign(1.0_dp,LISM_vh) .lt. 0.0_dp) ) .or. & 
            ( (zi .gT. -120.0_dp*au) .AND. (xyr .LT. p_xmax/2.0_dp) &
            .and. (sign(1.0_dp,LISM_vh) .gt. 0.0_dp) ) ) THEN
! don't recombine too much in the TAIL
          CALL plasma(xyr,zi,n_p,T_p,ux,uz,region_p)
          IF ( region_p .GT. 1.9_dp ) THEN
             ideal_massi = 0.2_dp*ideal_massi
          END IF
       END IF

       IF ( (massi .LT. 0.5_dp*ideal_massi) &
            .AND. (region .LT. 1.5_dp) ) THEN
! don't recombine comp2 and comp3 neutrals
          CALL rand_number(dprand)
          IF (dprand .LT. 0.5_dp) THEN
! double mass
             my_particles(my_thread)%particles(part,0) = massi * 2.0_dp
          ELSE
! dissintegrate
             my_particles(my_thread)%particles(part,:) = my_particles(my_thread)%particles(my_particles(my_thread)%my_npart,:)
             my_particles(my_thread)%my_npart = my_particles(my_thread)%my_npart - 1
          END IF
       END IF
       
    END DO

    my_particles(my_thread)%particles(my_particles(my_thread)%my_npart+1:my_old_npart,:) = 0.0_dp

  END SUBROUTINE split
!
!***********************************************************************
!
  SUBROUTINE create_neutral_grid
    REAL(KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: H_dens, H_temp
    REAL(KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: tmp1, H_ux, H_uz
    INTEGER, DIMENSION (:,:,:), ALLOCATABLE :: tmp2, H_npart
    REAL(KIND=dp) :: mass, xi, yi, zi, xyr, xyr_sq, mass_x, area
    REAL(KIND=dp) :: vx, vy, vz, vxr, ux, uz
    REAL(KIND=dp) :: vrel_sq, vxy_sq, vyr_sq, xx, zz, dx
    INTEGER :: i, j, k, ii, jj, kkx, kkz, kk, ierr, region, nn, nx, nz, thread
    
! collate neutral velocity, density and temperature onto grid

    ALLOCATE(H_dens(ngrids,grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(H_temp(ngrids,grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(H_ux(ngrids,grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(H_uz(ngrids,grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(H_npart(ngrids,grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(tmp1(ngrids,grids_nxmax,-grids_nzmax+1:grids_nzmax))
    ALLOCATE(tmp2(ngrids,grids_nxmax,-grids_nzmax+1:grids_nzmax))

    DO region = -1, 3

       H_dens = 0.0_dp; H_temp = 0.0_dp
       H_ux = 0.0_dp; H_uz = 0.0_dp
       H_npart = 0
       
       do thread = 0, nthreads - 1

       DO i = 1, my_particles(thread)%my_npart

! output only particles born in region "region"
          IF ( (region .EQ. -1) .OR. (FLOOR(my_particles(thread)%particles(i,7)) .EQ. region) ) THEN

             xi = my_particles(thread)%particles(i,1)
             yi = my_particles(thread)%particles(i,2)
             zi = my_particles(thread)%particles(i,3)

             xyr_sq = xi*xi + yi*yi
! don't include boundary points
             IF ( (xyr_sq .LT. p_xmax_sq) .AND. (zi .LT. p_zmax) &
                  .and. (zi .gt. p_zmin) ) THEN

                mass = my_particles(thread)%particles(i,0)
                vx = my_particles(thread)%particles(i,4)
                vy = my_particles(thread)%particles(i,5)
                vz = my_particles(thread)%particles(i,6)

                xyr = SQRT(xyr_sq)
! find projection of vxy into xyr using dot product
                vxr = (vx*xi + vy*yi) / xyr
          
! find which grid we're on and nearest grid point
                kkx = -1
                DO nn = 1, ngrids
                   IF (abs(xyr) .GE. grids_boundaries_x(ngrids-nn)) THEN
                      kkx = ngrids-nn+1
                      EXIT
                   END IF
                END DO
                kkz = -1
                DO nn = 1, ngrids
                   IF (abs(zi) .GE. grids_boundaries_z(ngrids-nn)) THEN
                      kkz = ngrids-nn+1
                      EXIT
                   END IF
                END DO
                kk = max(kkx,kkz)
                ii = FLOOR(xyr / grids_dx(kk)) + 1
                jj = FLOOR(zi / grids_dx(kk)) + 1

! must allow for linear variation in mass in x-direction
                mass_x = mass / (xyr*twopi)

                H_npart(kk,ii,jj) = H_npart(kk,ii,jj) + 1
                H_dens( kk,ii,jj) = H_dens( kk,ii,jj) +       mass_x
                H_ux(   kk,ii,jj) = H_ux(   kk,ii,jj) + vxr * mass_x
                H_uz(   kk,ii,jj) = H_uz(   kk,ii,jj) + vz  * mass_x

             END IF
             
          END IF
 
       END DO
       enddo

! temporary storage for REDUCE
       tmp1 = 0.0_dp
       tmp2 = 0

! add together npart and store on MASTER
       CALL MPI_Reduce(H_npart, tmp2, total_ndata, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       H_npart = tmp2

! add together all density data
       CALL MPI_ALLReduce(H_dens, tmp1, total_ndata, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       H_dens = tmp1

! add together all velocity data
       CALL MPI_ALLReduce(H_ux, tmp1, total_ndata, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       H_ux = tmp1

       CALL MPI_ALLReduce(H_uz, tmp1, total_ndata, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       H_uz = tmp1
       

!       WHERE (H_dens /= 0.0_dp)
!          H_ux = H_ux / H_dens
!          H_uz = H_uz / H_dens
!       END WHERE

        WHERE (H_dens > 1d-10)
          H_ux = H_ux / H_dens
          H_uz = H_uz / H_dens
       END WHERE


       H_temp = 0.0_dp
       
       do thread = 0, nthreads - 1
       DO i = 1, my_particles(thread)%my_npart

! output only particles born in region "region"
          IF ( (region .EQ. -1) .OR. (FLOOR(my_particles(thread)%particles(i,7)) .EQ. region) ) THEN

             xi = my_particles(thread)%particles(i,1)
             yi = my_particles(thread)%particles(i,2)
             zi = my_particles(thread)%particles(i,3)

             xyr_sq = xi*xi + yi*yi
! don't include boundary points
             IF ( (xyr_sq .LT. p_xmax_sq) .AND. (zi .LT. p_zmax) &
                  .and. (zi .gt. p_zmin) ) THEN

                mass = my_particles(thread)%particles(i,0)
                vx = my_particles(thread)%particles(i,4)
                vy = my_particles(thread)%particles(i,5)
                vz = my_particles(thread)%particles(i,6)

                xyr = SQRT(xyr_sq)
! find projection of vxy into xyr using dot product
                vxr = (vx*xi + vy*yi) / xyr

                vxy_sq = vx*vx + vy*vy
                vyr_sq = vxy_sq - vxr*vxr

! find which grid we're on and nearest grid point
                kkx = -1
                DO nn = 1, ngrids
                   IF (abs(xyr) .GE. grids_boundaries_x(ngrids-nn)) THEN
                      kkx = ngrids-nn+1
                      EXIT
                   END IF
                END DO
                kkz = -1
                DO nn = 1, ngrids
                   IF (abs(zi) .GE. grids_boundaries_z(ngrids-nn)) THEN
                      kkz = ngrids-nn+1
                      EXIT
                   END IF
                END DO
                kk = max(kkx,kkz)
                ii = FLOOR(xyr / grids_dx(kk)) + 1
                jj = FLOOR(zi / grids_dx(kk)) + 1

! must allow for linear variation in mass in x-direction
                mass_x = mass / (xyr*twopi)

                ux = H_ux(kk,ii,jj)
                uz = H_uz(kk,ii,jj)

! calculate temperature
                vrel_sq = (vxr-ux)**2 + (vz-uz)**2 + vyr_sq

                H_temp(kk,ii,jj) = H_temp(kk,ii,jj) + vrel_sq * mass_x

             END IF
             
          END IF

       END DO
       enddo

! add together all vrel_sq data
       CALL MPI_Reduce(H_temp, tmp1, total_ndata, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       H_temp = tmp1

       IF (my_rank .EQ. 0) THEN
          WHERE (H_dens /= 0.0_dp)
             H_temp = neutral_mass * H_temp / (3.0_dp*k_B*H_dens)
          END WHERE
! total values are only stored on MASTER
          DO kk = 1, ngrids
             nx = grids_nx(kk)
             nz = grids_nz(kk)
             area = grids_dx(kk)**2
             grid_dens(kk,region,1:nx,-nz+1:nz) = H_dens(kk,1:nx,-nz+1:nz) &
                  * mega_number / area
             grid_ux(kk,region,1:nx,-nz+1:nz) = H_ux(kk,1:nx,-nz+1:nz)
             grid_uz(kk,region,1:nx,-nz+1:nz) = H_uz(kk,1:nx,-nz+1:nz)
             grid_temp(kk,region,1:nx,-nz+1:nz) = H_temp(kk,1:nx,-nz+1:nz)
             grid_npart(kk,region,1:nx,-nz+1:nz) = H_npart(kk,1:nx,-nz+1:nz)
          END DO
       END IF

    END DO  ! region loop


    DEALLOCATE(H_dens, H_temp, H_ux, H_uz, H_npart, tmp1, tmp2)
    

! master will fill in all grids -- this is actually NOT NECESSARY,
! but helps with plotting.
    IF (my_rank .EQ. 0) THEN
       
       DO region = -1, 3

          DO k = 1, ngrids
             nx = grids_nx(k)
             nz = grids_nz(k)
             dx = grids_dx(k)
             DO j = -nz+1, nz
                zz = (j-0.5_dp)*dx
                DO i = 1, nx
                   xx = (i-0.5_dp)*dx
                   
! find which grid we're on and nearest grid point
                   kkx = -1
                   DO nn = 1, ngrids
                      IF (abs(xx) .GE. grids_boundaries_x(ngrids-nn)) THEN
                         kkx = ngrids-nn+1
                         EXIT
                      END IF
                   END DO
                   kkz = -1
                   DO nn = 1, ngrids
                      IF (abs(zz) .GE. grids_boundaries_z(ngrids-nn)) THEN
                         kkz = ngrids-nn+1
                         EXIT
                      END IF
                   END DO
                   kk = max(kkx,kkz)
                   ii = FLOOR(xx / grids_dx(kk)) + 1
                   jj = FLOOR(zz / grids_dx(kk)) + 1

                   grid_dens(k,region,i,j)  = grid_dens(kk,region,ii,jj)
                   grid_ux(k,region,i,j)    = grid_ux(kk,region,ii,jj)
                   grid_uz(k,region,i,j)    = grid_uz(kk,region,ii,jj)
                   grid_temp(k,region,i,j)  = grid_temp(kk,region,ii,jj)
                   grid_npart(k,region,i,j) = grid_npart(kk,region,i,jj)

                END DO
             END DO
          END DO

       END DO

  END IF

  

  END SUBROUTINE create_neutral_grid
!
!***********************************************************************
!
  SUBROUTINE output_raw(rawfile)
    use hdf5
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: rank, ierr, i
    integer(8) :: my_npart,total_npart
    integer(8), DIMENSION (0:nprocMPI-1) :: npart
    CHARACTER (LEN=60) :: rawfile
    
    ! hdf5 IO
    INTEGER :: info
    INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension    
    INTEGER(HID_T) :: file_id       ! File identifier 
    INTEGER(HID_T) :: plist_id      ! Property list identifier 	
    INTEGER(HID_T) :: attr_id       ! Attribute identifier 
    INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier 
    INTEGER(HID_T) :: atype_id
    INTEGER(HID_T) :: gr_id
    INTEGER(HID_T) :: dset_id       ! Dataset identifier 
    INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
    INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
    
    INTEGER(HSIZE_T),  DIMENSION(2) :: dimsf
    INTEGER(HSIZE_T),  DIMENSION(2) :: my_count  
    INTEGER(HSIZE_T), DIMENSION(2)  :: my_offset             
    
    REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: tmp_particles
    
    !INTEGER(HID_T) :: int64
    
    !call H5Tcopy_f(H5T_NATIVE_INTEGER,int64,ierr)
    !call H5Tset_precision_f(int64, 64_8, ierr)
    
        
    
    my_npart = 0
    do i = 0, nthreads - 1
      my_npart = my_npart + my_particles(i)%my_npart
    enddo
      
    call MPI_allgather(my_npart, 1, MPI_LONG_LONG_INT, npart, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD, ierr)
    
    total_npart = sum(npart)
    
    if (total_npart == 0) then
      return
    endif
    
    info = MPI_INFO_NULL
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, ierr)
    
    CALL h5fcreate_f(rawfile, H5F_ACC_TRUNC_F, file_id, ierr, access_prp = plist_id)
    CALL h5pclose_f(plist_id, ierr)
    
    call h5gopen_f(file_id, '/', gr_id, ierr)
    
    CALL h5screate_f(H5S_SCALAR_F,aspace_id,ierr)    
    
    
    
    
    
    ! Create the data space for the  dataset. 
    !dimsf(1) = 1
    !dimsf(2) = 0
    !CALL h5screate_simple_f(1, dimsf, filespace, ierr)          
    !CALL h5dcreate_f(file_id, 'total_npart', int64, filespace, dset_id, ierr)
    !CALL h5sclose_f(filespace, ierr)
    !CALL h5dwrite_f(dset_id, int64, total_npart, dimsf, ierr)
    !CALL h5dclose_f(dset_id, ierr)    
        
    !CALL h5acreate_f(gr_id, 'total_npart', int64, aspace_id, attr_id, ierr)
    !CALL h5awrite_f(attr_id, int64, total_npart, adims, ierr)     
    !CALL h5aclose_f(attr_id, ierr)
    
    CALL h5acreate_f(gr_id, 'p_xmax', H5T_NATIVE_DOUBLE, aspace_id, attr_id, ierr)
    CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, p_xmax, adims, ierr) 
    CALL h5aclose_f(attr_id, ierr)
        
    CALL h5acreate_f(gr_id, 'boundary_layer', H5T_NATIVE_DOUBLE, aspace_id, attr_id, ierr)
    CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, boundary_layer, adims, ierr) 
    CALL h5aclose_f(attr_id, ierr)
    
    CALL h5acreate_f(gr_id, 'LISM_v_thermal', H5T_NATIVE_DOUBLE, aspace_id, attr_id, ierr)
    CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, LISM_v_thermal, adims, ierr) 
    CALL h5aclose_f(attr_id, ierr)
    
    CALL h5acreate_f(gr_id, 'time', H5T_NATIVE_DOUBLE, aspace_id, attr_id, ierr)
    CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, time, adims, ierr) 
    CALL h5aclose_f(attr_id, ierr)
    
    CALL h5sclose_f(aspace_id, ierr)
    call h5gclose_f(gr_id, ierr)
    
    
    dimsf(1) = total_npart
    dimsf(2) = 8
        
    
    ! Create the data space for the  dataset. 
    CALL h5screate_simple_f(2, dimsf, filespace, ierr)      
      ! Create the dataset with default properties.      
    CALL h5dcreate_f(file_id, 'particles', H5T_NATIVE_DOUBLE, filespace, dset_id, ierr)
    CALL h5sclose_f(filespace, ierr)
        
    my_count(2)  = dimsf(2)
    my_offset(2) = 0
    
    my_offset(1) = 0
    if (my_rank>0) then
      my_offset(1) = sum(npart(0:my_rank-1))
    endif
    
    ! Create property list for collective dataset write      
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
    
    print *,'output_raw, total number of particles',total_npart
    
    do i = 0, nthreads - 1
      my_count(1) = my_particles(i)%my_npart
      if (my_count(1).eq.0) then
        cycle
      endif
      
      print *,'         thread',i,'number of particles ', my_count(1)
      
      
      
      !allocate(tmp_particles(1:my_count(1),0:7))
      !tmp_particles(1:my_count(1),0:7) = my_particles(i)%particles(1:my_count(1),0:7)
      
      
      
      CALL h5screate_simple_f(2, my_count, memspace, ierr) 
      if (ierr < 0) then
        print *, 'h5screate_simple_f error'
      endif
             
      ! Select hyperslab in the file.      
      CALL h5dget_space_f(dset_id, filespace, ierr)
      if (ierr < 0) then
        print *, 'h5dget_space_f error'
      endif
      
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,  my_offset, my_count, ierr)
      if (ierr < 0) then
        print *, 'h5sselect_hyperslab_f error'
      endif
                             
      ! Write the dataset collectively.       
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, my_particles(i)%particles(1:my_count(1),0:7), dimsf, ierr, &
                     file_space_id = filespace, &
                     mem_space_id = memspace)
                     
      !CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tmp_particles, dimsf, ierr, &
      !               file_space_id = filespace, &
      !               mem_space_id = memspace, xfer_prp = plist_id)
                            
          
      ! Close dataspaces.     
      CALL h5sclose_f(filespace, ierr)
      CALL h5sclose_f(memspace, ierr)            
      
      !deallocate(tmp_particles)
      
      my_offset(1) = my_offset(1) + my_particles(i)%my_npart
    
    enddo
        
    
    ! Close the dataset.      
    CALL h5dclose_f(dset_id, ierr)    
    
    ! Close property list.      
    CALL h5pclose_f(plist_id,  ierr)
    
    CALL h5fclose_f(file_id,  ierr)
    
    
    
!~
!~    
!~
!~    
!~  CALL MPI_Reduce(my_npart,total_npart,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!~
!~    IF (my_rank .EQ. 0) THEN
!~
!~       OPEN(unit=1,file=rawfile,form='unformatted',status='unknown',action='write')
!~       WRITE(1) total_npart
!~       WRITE(1) p_xmax
!~       WRITE(1) boundary_layer
!~       WRITE(1) LISM_v_thermal
!~       WRITE(1) time
!~       WRITE(1) nprocMPI
!~       WRITE(1) 1.0_dp*my_npart + 0.1_dp
!~       IF ( verbosity .ge. 3 ) THEN
!~          PRINT *, 'writing my own data, ',my_npart
!~       end IF
!~       WRITE(1) particle(1:my_npart,0:7)
!~
!~       DO rank = 1, nprocMPI-1
!~
!~          CALL MPI_RECV(npart, 1, MPI_INTEGER, rank, &
!~                           1000+rank, MPI_COMM_WORLD, status, ierr)
!~          CALL MPI_RECV(particle(1:npart,0:7), 8*npart, MPI_DOUBLE_PRECISION, &
!~                           rank, 2000+rank, MPI_COMM_WORLD, status, ierr)
!~
!~          WRITE(1) 1.0_dp*npart + 0.1_dp
!~          IF ( verbosity .ge. 3 ) THEN
!~             PRINT *, 'writing data from proc ',rank,', ',npart
!~          end IF
!~          WRITE(1) particle(1:npart,0:7)
!~
!~          CALL release_cache()
!~
!~       END DO
!~       CLOSE(1)
!~
!~    ELSE
!~
!~       CALL MPI_SEND(my_npart, 1, MPI_INTEGER, 0, 1000+my_rank, &
!~                     MPI_COMM_WORLD, ierr)
!~       CALL MPI_SEND(particle(1:my_npart,0:7), 8*my_npart, &
!~                     MPI_DOUBLE_PRECISION, 0, 2000+my_rank, MPI_COMM_WORLD, ierr)
!~
!~    END IF
!~
!~9 FORMAT(6e13.4)

  END SUBROUTINE output_raw
!
!***********************************************************************
!
  SUBROUTINE load_raw
    use hdf5
    REAL(KIND=dp) :: load_rmax, load_boundary_layer, load_LISM_v_thermal
    REAL(KIND=dp) :: real_my_npart
    REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: temp
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER(8) :: load_total_npart, new_npart
    INTEGER :: blocksize, read_npart
    INTEGER :: rank, ierr, i
    CHARACTER (LEN = 100) str_buf
    
    ! hdf5 IO
    INTEGER :: info
    INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension    
    INTEGER(HID_T) :: file_id       ! File identifier 
    INTEGER(HID_T) :: plist_id      ! Property list identifier 	
    INTEGER(HID_T) :: attr_id       ! Attribute identifier 
    INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier 
    INTEGER(HID_T) :: atype_id
    INTEGER(HID_T) :: gr_id
    INTEGER(HID_T) :: dset_id       ! Dataset identifier 
    INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
    INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
    
    INTEGER(HSIZE_T),  DIMENSION(2) :: dimsf
    INTEGER(HSIZE_T),  DIMENSION(2) :: maxdimsf
    INTEGER(HSIZE_T),  DIMENSION(2) :: my_count  
    INTEGER(HSIZE_T), DIMENSION(2)  :: my_offset             
    
    
    !info = MPI_INFO_NULL
    !CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    !CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, ierr)
     
    ! Open the file collectively.       
    !CALL h5fopen_f(load_file, H5F_ACC_RDONLY_F, file_id, ierr, access_prp = plist_id)    
    !CALL h5pclose_f(plist_id, ierr)
    
    CALL h5fopen_f(load_file, H5F_ACC_RDONLY_F, file_id, ierr)

    call h5gopen_f(file_id, '/', gr_id, ierr)

    CALL h5aopen_name_f(gr_id, 'time', attr_id, ierr)
    CALL h5aread_f(attr_id, H5T_NATIVE_DOUBLE, time, adims, ierr) 
    CALL h5aclose_f(attr_id, ierr)

    call h5gclose_f(gr_id, ierr)
    
    CALL h5dopen_f(file_id, 'particles', dset_id, ierr)            
    CALL h5dget_space_f(dset_id, filespace, ierr)
    call H5Sget_simple_extent_dims_f(filespace, dimsf, maxdimsf, ierr)
    
    load_total_npart = dimsf(1)
    new_npart = load_total_npart/(nprocMPI*nthreads)
    
    my_count(2)  = dimsf(2)
    my_offset(2) = 0
    
    
    my_offset(1) = my_rank*nthreads*new_npart
    
    ! Create property list for collective dataset write      
    !CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr) 
    !CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    
    do i = 0, nthreads - 1
      my_count(1) = new_npart
      if ((my_rank.eq.(nprocMPI-1)).and.(i.eq.(nthreads - 1))) then
        my_count(1) = load_total_npart - (nthreads*nprocMPI - 1)*new_npart 
        !my_count(1) = load_total_npart - (nthreads*(nprocMPI - 1)*new_npart + (nthreads - 1)*new_npart)
      endif
      
      !allocate(tmp_particles(1:my_count(1),0:7))
      !tmp_particles(:,:) = my_particles(i)%particles(1:my_count(1),0:7)
            
      
      write (str_buf,*) 'read [',my_offset(1),'..',my_offset(1)+my_count(1) - 1,'] particles'
      call pout(str_buf)      
      
      print *, 'read [',my_offset(1),'..',my_offset(1)+my_count(1)-1,'] particles'
      
      CALL h5screate_simple_f(2, my_count, memspace, ierr) 
             
      ! Select hyperslab in the file.      
      CALL h5dget_space_f(dset_id, filespace, ierr)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,  my_offset, my_count, ierr)
                             
      
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, my_particles(i)%particles(1:my_count(1),0:7), dimsf, ierr, &
                     file_space_id = filespace, &
                     mem_space_id = memspace)
                     
      !CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, my_particles(i)%particles(1:my_count(1),0:7), dimsf, ierr, &
      !               file_space_id = filespace, &
      !               mem_space_id = memspace, xfer_prp = plist_id)
                     
                            
          
      ! Close dataspaces.     
      CALL h5sclose_f(filespace, ierr)
      CALL h5sclose_f(memspace, ierr)            
      
      my_offset(1) = my_offset(1) + new_npart
      
      my_particles(i)%my_npart = my_count(1)
      
      !deallocate(tmp_particles)
    
    enddo
    
    ! Close property list.      
    !CALL h5pclose_f(plist_id,  ierr)
    
    ! Close the dataset.      
    CALL h5dclose_f(dset_id, ierr)    
    
    CALL h5fclose_f(file_id,  ierr)
    
    
    

!~    IF (my_rank .EQ. 0) THEN
!~       OPEN(unit=1,file=load_file,form='unformatted',status='old',action='read')
!~       READ(1) load_total_npart
!~
!~       IF ( verbosity .ge. 2 ) THEN
!~          PRINT *,'loading ',load_total_npart,' particles'
!~       end IF
!~       IF (load_total_npart .GT. my_ntotal*nprocMPI) THEN
!~          PRINT *,'too many particles, STOPPING'
!~          STOP
!~       END IF
!~
!~       READ(1) load_rmax
!~       READ(1) load_boundary_layer
!~       READ(1) load_LISM_v_thermal
!~       READ(1) time
!~       READ(1) load_nproc
!~
!~       blocksize = 3*MAX(load_total_npart/load_nproc,my_ntotal)
!~       ALLOCATE(temp(1:blocksize,0:7))
!~
!~       new_npart = load_total_npart/nprocMPI   ! how many each proc will receive
!~       read_npart = 0
!~       DO rank = 1, nprocMPI-1
!~
!~          DO WHILE (read_npart .LT. new_npart)
!~             READ(1) real_my_npart
!~             my_npart = FLOOR(real_my_npart)
!~             IF ( verbosity .ge. 3 ) THEN
!~                PRINT *, 'reading ', my_npart, ' particles'
!~             end IF
!~             READ(1) temp(read_npart+1:read_npart+my_npart,0:7)
!~             read_npart = read_npart + my_npart
!~          END DO
!~ 
!~          CALL MPI_SEND(new_npart, 1, MPI_INTEGER, rank, 1000+rank, &
!~               MPI_COMM_WORLD, ierr)
!~          IF ( verbosity .ge. 3 ) THEN
!~             PRINT *, 'sending ',new_npart,' particles to proc # ',rank
!~          end IF
!~          CALL MPI_SEND(temp(1:new_npart,0:7), 8*new_npart, &
!~               MPI_DOUBLE_PRECISION, rank, 2000+rank, MPI_COMM_WORLD, ierr)
!~ 
!~          temp(1:read_npart-new_npart,0:7) = temp(new_npart+1:read_npart,0:7)
!~          read_npart = read_npart - new_npart
!~
!~          CALL release_cache()
!~
!~       END DO
!~
!~! master takes particles left over from rounding
!~       new_npart = new_npart + load_total_npart - nprocMPI*new_npart
!~       DO WHILE (read_npart .LT. new_npart)
!~          READ(1) real_my_npart
!~          my_npart = FLOOR(real_my_npart)
!~         
!~          READ(1) temp(read_npart+1:read_npart+my_npart,0:7)
!~          read_npart = read_npart + my_npart
!~       END DO
!~
!~       IF ( verbosity .ge. 3 ) THEN
!~          PRINT *, 'keeping remaining ', new_npart, ' particles for myself'
!~       end IF
!~       my_npart = new_npart
!~       particle(1:my_npart,0:7) = temp(1:new_npart,0:7)
!~
!~       CLOSE(1)
!~       DEALLOCATE(temp)
!~       CALL release_cache()
!~
!~    ELSE
!~
!~       DO rank = 1, nprocMPI-1
!~
!~          IF (my_rank .EQ. rank) THEN
!~             CALL MPI_RECV(my_npart, 1, MPI_INTEGER, 0, 1000+rank, &
!~                  MPI_COMM_WORLD, status, ierr)
!~             CALL MPI_RECV(particle(1:my_npart,0:7), 8*my_npart, &
!~                  MPI_DOUBLE_PRECISION, 0, 2000+rank, MPI_COMM_WORLD, status, ierr)
!~             
!~          END IF
!~       END DO
!~
!~    END IF
!~
!~    CALL MPI_BCAST(time, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!~
!~    particle(my_npart+1:my_ntotal,:) = 0.0_dp
!~
  END SUBROUTINE load_raw
!
!***********************************************************************
!





 subroutine partion_qs(ba,size_ba,p,r,pivot_ind)
  integer, intent(in) :: size_ba
  type(node_workload), DIMENSION (0:size_ba-1) :: ba  
  integer, intent(in)  :: p,r
  integer, intent(out) :: pivot_ind
  integer :: i,j
  
  type(node_workload) :: x, tmp
  
  x = ba(r)
  i = p - 1
  do j=p, r-1
    if (ba(j)%part_disbalance <= x%part_disbalance) then
      i = i + 1
      tmp = ba(i)
      ba(i) = ba(j)
      ba(j) = tmp
    endif  
  enddo
  
  tmp = ba(i+1)
  ba(i+1) = ba(r)
  ba(r) = tmp
  
  pivot_ind = i + 1

end subroutine partion_qs

recursive subroutine quicksort(ba,size_ba,p,r)
  integer, intent(in) :: size_ba
  type(node_workload), DIMENSION (0:size_ba-1) :: ba  
  integer, intent(in)  :: p,r
  
  integer :: q
  
  if (p<r) then
    call partion_qs(ba,size_ba,p,r,q)
    call quicksort(ba,size_ba,p,q-1)
    call quicksort(ba,size_ba,q+1,r)
  endif

end subroutine quicksort

! sort mpi_balancing structure. Disblanced nodes will appear at the end of the array
  subroutine sort_balancing_array(ba,size_ba)
    integer, intent(in) :: size_ba
    type(node_workload), DIMENSION (0:size_ba-1) :: ba    
  
    integer :: i,ind1
    type(node_workload) :: tmp
    CHARACTER (LEN = 100) str_buf
    
    if (size_ba .eq. 1) then
      return
    endif
    
!    write (str_buf,*) '    sort_balancing_array, non sorted array:'
!    call pout(str_buf) 
!    do i = 0,size_ba-1
!      write (str_buf,*) 'procID = ',ba(i)%procID,' disbalance = ',ba(i)%part_disbalance
!      call pout(str_buf)            
!    enddo
    
    
    call quicksort(ba,size_ba,0,size_ba-1)
    
    ! reverse the array
    do i = 0,size_ba/2-1
      ind1 = size_ba - 1 - i
      tmp = ba(ind1)
      ba(ind1) = ba(i)
      ba(i) = tmp
    enddo
    
!    write (str_buf,*) '    sort_balancing_array, sorted array:'
!    call pout(str_buf) 
!    do i = 0,size_ba-1
!      write (str_buf,*) 'procID = ',ba(i)%procID,' disbalance = ',ba(i)%part_disbalance
!      call pout(str_buf)            
!    enddo
  
  end subroutine sort_balancing_array
  
  subroutine add_particles_node(buffer, bnpart)
    REAL(KIND=dp), intent(in) :: buffer(1:,0:)
    integer,       intent(in) :: bnpart
    
    integer :: i, sP, iSenShift, avail
    
    sP        = bnpart
    iSenShift = 1
        
    do i = 0, nthreads - 1
    
      avail = my_ntotal - my_particles(i)%my_npart
      
      if (avail <= 0) then
        cycle
      endif
      
      if (sP<=avail) then
        my_particles(i)%particles(my_particles(i)%my_npart+1:my_particles(i)%my_npart+sP,:) = & 
                buffer(iSenShift:iSenShift+sP-1,:)        
        my_particles(i)%my_npart = my_particles(i)%my_npart + sP
        exit
      else
        my_particles(i)%particles(my_particles(i)%my_npart+1:my_ntotal,:) = buffer(iSenShift:iSenShift+avail-1,:)        
                
        iSenShift  = iSenShift + avail
        sP         = sP        - avail
        
        my_particles(i)%my_npart =  my_ntotal
      endif
    enddo         
  
  end subroutine add_particles_node
  
  subroutine balance_node()
    integer(8) :: my_npart
    integer    :: npart_average, ithread, jthread   
    type(node_workload), DIMENSION (0:nthreads-1) :: node_balancing
        
    integer :: iSenShift,sP,iRecShift,rP,iR,iS,tmp
    
    
    
    if (nthreads == 1) then
      return
    endif
  
    my_npart = 0
    do ithread = 0, nthreads - 1
      my_npart = my_npart + my_particles(ithread)%my_npart
    enddo      
                
    npart_average = my_npart/nthreads + 1    
    
    !print *,'particles on node: ',my_npart,'average: ',npart_average 
            
    do ithread = 0, nthreads - 1
      node_balancing(ithread)%procID = ithread;
      node_balancing(ithread)%part_disbalance = my_particles(ithread)%my_npart - npart_average   

     ! print *,'thread ',ithread, 'particles ',my_particles(ithread)%my_npart,'disbalance ', node_balancing(ithread)%part_disbalance
    enddo
    
    call sort_balancing_array(node_balancing,nthreads)
    
    jthread    = nthreads - 1
    iR         = node_balancing(jthread)%procID
    rP         = min(abs(node_balancing(jthread)%part_disbalance), my_ntotal - my_particles(iR)%my_npart)
    iRecShift  = my_particles(iR)%my_npart+1 ! shift in receiver buffer        
    
    do ithread = 0, nthreads - 1    
    
      if (node_balancing(ithread)%part_disbalance <= 0)  then
        exit
      endif

      iS = node_balancing(ithread)%procID          ! sender id
      sP = node_balancing(ithread)%part_disbalance ! sender number of particles
      iSenShift = my_particles(node_balancing(ithread)%procID)%my_npart
            
      do while (.true.)
      
        if (rP.eq.0) then ! go to the next thread, this is full
          jthread    = jthread  - 1
          
          if (jthread <= 0) then
            exit
          endif
          
          if (node_balancing(jthread)%part_disbalance >= 0)  then
            exit
          endif
          
          iR         = node_balancing(jthread)%procID
          rP         = min(abs(node_balancing(jthread)%part_disbalance), my_ntotal - my_particles(iR)%my_npart)
          if (rP.eq.0) then
            cycle
          endif
          iRecShift  = my_particles(iR)%my_npart+1 
        endif
              
        
        if (sP <= rP) then          
          
          my_particles(iR)%particles(iRecShift:iRecShift+sP-1,:) = my_particles(iS)%particles(iSenShift-sP+1:iSenShift,:)           
          
          
          !print *,'thread',iS,'sends [',iSenShift-sP+1,iSenShift,'] to ',iR,' [',iRecShift,iRecShift+sP-1,']'
          
          my_particles(iS)%my_npart = my_particles(iS)%my_npart-sP
          my_particles(iR)%my_npart = my_particles(iR)%my_npart+sP
          
          iRecShift  = my_particles(iR)%my_npart+1
                    
          rP = rP - sP            
          
          

          
          
          exit
        else                  
            
          my_particles(iR)%particles(iRecShift:iRecShift+rP-1,:) = my_particles(iS)%particles(iSenShift-rP+1:iSenShift,:)
          
          !print *,'thread',iS,'sends [',iSenShift-rP+1,iSenShift,'] to ',iR,' [',iRecShift,iRecShift+rP-1,']'
            

          my_particles(iS)%my_npart = my_particles(iS)%my_npart-rP
          my_particles(iR)%my_npart = my_particles(iR)%my_npart+rP
          
          sP        = sP - rP
                    
                  
          iSenShift = iSenShift - rP
          jthread   = jthread  - 1
          iR         = node_balancing(jthread)%procID
          rP         = min(abs(node_balancing(jthread)%part_disbalance), my_ntotal - my_particles(iR)%my_npart)
          iRecShift  = my_particles(iR)%my_npart+1 
        endif        
      enddo  
      
      if (node_balancing(jthread)%part_disbalance >= 0)  then
        exit
      endif

    
        
    enddo
    
    !my_npart = 0
    !do ithread = 0, nthreads - 1      
    !  my_npart = my_npart + my_particles(ithread)%my_npart
    !  print *,'thread ',ithread, 'particles ',my_particles(ithread)%my_npart
    !enddo
    !print *,'particles on node: ',my_npart
     
  end subroutine balance_node
  

  SUBROUTINE balance_npart()    
    real(KIND=dp) :: npart_ratio
    real, parameter :: max_ratio = 0.1d0    
    INTEGER :: tmp, min_rank, max_rank, tot_min, tot_max, ierr
    
    integer(8) :: my_npart ! number of particles on node
    integer(8) :: total_npart,npart_average
    integer(8), DIMENSION (0:nprocMPI-1) :: npart
    
    integer :: my_disbalance, nrecv_part
    integer :: iproc,jproc,i
    integer :: iSenShift,sP,iRecShift,rP,iR,iS
    
    REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: buffer
    CHARACTER (LEN = 100) str_buf
    
        
        
    my_npart = 0
    do i = 0, nthreads - 1
      my_npart = my_npart + my_particles(i)%my_npart
    enddo      
        
    
    call MPI_allgather(my_npart, 1, MPI_LONG_LONG_INT, npart, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD, ierr)    
        
    write (str_buf,*) '             , my particles before: ', my_npart
    call pout(str_buf) 
    
    total_npart   = sum(npart)    
    npart_average = total_npart/nprocMPI + 1
    
    IF (my_rank .eq. 0) then
      print*,'total number of particles ', total_npart
    endif
        
    do iproc = 0, nprocMPI - 1
      mpi_balancing(iproc)%procID = iproc;
      mpi_balancing(iproc)%part_disbalance = npart(iproc) - npart_average            
    enddo
    
    my_disbalance     = mpi_balancing(my_rank)%part_disbalance
    
    ! allocate sender/receiver buffer
    if (my_disbalance .ne.  0) then
      allocate(buffer(1:abs(my_disbalance),0:7))
      buffer = 0d0
    endif
    
    npart_ratio = real(abs(my_disbalance))/real(npart_average)
      
    !if ((my_disbalance > 0).and.(npart_ratio>max_ratio)) then
    if (my_disbalance > 0) then
      sP        = my_disbalance
      iSenShift = 1
      
      do i = 0, nthreads - 1
        if (sP<=my_particles(i)%my_npart) then
          buffer(iSenShift:iSenShift+sP-1,:) = my_particles(i)%particles(my_particles(i)%my_npart-sP+1:my_particles(i)%my_npart,:)
          my_particles(i)%my_npart = my_particles(i)%my_npart - sP
          exit
        else
          buffer(iSenShift:iSenShift+my_particles(i)%my_npart-1,:) = my_particles(i)%particles(1:my_particles(i)%my_npart,:)
          iSenShift  = iSenShift + my_particles(i)%my_npart
          sP         = sP        - my_particles(i)%my_npart
          
          my_particles(i)%my_npart =  0
        endif
      enddo          
    endif
    
    call sort_balancing_array(mpi_balancing,nprocMPI)
    
    jproc      = nprocMPI - 1
    iR         = mpi_balancing(jproc)%procID 
    rP         = abs(mpi_balancing(jproc)%part_disbalance)
    iRecShift  = 1 ! shift in receiver buffer    
    nrecv_part = 0
    
    do iproc = 0, nprocMPI - 1
      npart_ratio = real(mpi_balancing(iproc)%part_disbalance)/real(npart_average)
      
      if (npart_ratio<=0) then
        exit
      endif
      
!      if (npart_ratio<=max_ratio) then
!        exit
!      endif
      
      iS = mpi_balancing(iproc)%procID          ! sender id
      sP = mpi_balancing(iproc)%part_disbalance ! sender number of particles
      iSenShift = 1
      
      do while (.true.)
        
        if (mpi_balancing(jproc)%part_disbalance >= 0) then
          write (str_buf,*) 'node balancing error'
          call pout(str_buf) 
          call mpi_abort()
        endif
        
        if (sP <= rP) then
          if (my_rank == iS) then        
            CALL MPI_SEND(buffer(iSenShift:iSenShift+sP-1,:), &
                     8*sP,MPI_DOUBLE_PRECISION,iR, &
                     0,MPI_COMM_WORLD,ierr)
          elseif (my_rank == iR) then  
            CALL MPI_RECV(buffer(iRecShift:iRecShift+sP-1,:), &
                    8*sP,MPI_DOUBLE_PRECISION,iS, &
                    0,MPI_COMM_WORLD,status,ierr)
            nrecv_part = nrecv_part + sP 
          endif
          rP = rP - sP
          iRecShift = iRecShift + sP
          exit
        else        
          if (my_rank == iS) then        
            CALL MPI_SEND(buffer(iSenShift:iSenShift+rP-1,:), &
                     8*rP,MPI_DOUBLE_PRECISION,iR, &
                     0,MPI_COMM_WORLD,ierr)
          elseif (my_rank == iR) then  
            CALL MPI_RECV(buffer(iRecShift:iRecShift+rP-1,:), &
                    8*rP,MPI_DOUBLE_PRECISION,iS, &
                    0,MPI_COMM_WORLD,status,ierr)
            nrecv_part = nrecv_part + rP 
          endif
          sP        = sP - rP
          
          iSenShift = iSenShift + rP
          jproc     = jproc - 1
          iR        = mpi_balancing(jproc)%procID 
          rP        = abs(mpi_balancing(jproc)%part_disbalance)
          iRecShift = 1         
        endif
        
      enddo      
          
    enddo
    
    if (nrecv_part > 0) then
      call add_particles_node(buffer,nrecv_part)
    endif
    
    if (my_disbalance .ne.  0) then
      deallocate(buffer)
    endif
    
    call balance_node()
    
    my_npart = 0
    do i = 0, nthreads - 1
      my_npart = my_npart + my_particles(i)%my_npart
    enddo      
    
    write (str_buf,*) 'balance_npart, my particles after: ', my_npart
    call pout(str_buf) 
!~
!~    IF (nprocMPI .gt. 1) THEN   ! obviously....
!~
!~    CALL MPI_ALLReduce(my_npart,min_npart,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
!~    CALL MPI_ALLReduce(my_npart,max_npart,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
!~    npart_ratio = real(max_npart,dp) / real(min_npart,dp)
!~
!~    loop_count = 1
!~! keep send max_npart's excess particles to min_npart until acceptable ratio
!~    do while ( (npart_ratio .gt. max_ratio) .and. (loop_count .le. nprocMPI/2) )
!~       if ( (my_rank .eq. 0) .and. (verbosity .ge. 4) ) then
!~          print *, 'max_npart/min_npart = ', npart_ratio
!~       end if
!~
!~       IF (my_npart .EQ. min_npart) THEN
!~          tmp = 1
!~       ELSE
!~          tmp = 0
!~       END IF
!~
!~       CALL MPI_ALLReduce(tmp,tot_min,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!~
!~       IF (tot_min .gt. 1) THEN
!~! more than 1 node has the same min number of particles -- should resolve itself..
!~          exit
!~       else
!~! there is exactly one proc with the minimum number of particles
!~          CALL MPI_ALLReduce(my_npart,max_npart,1, &
!~               MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
!~       
!~          IF (my_npart .EQ. max_npart) THEN
!~             tmp = 1
!~          ELSE
!~             tmp = 0
!~          END IF
!~
!~          CALL MPI_ALLReduce(tmp,tot_max,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!~       
!~          IF (tot_max .gt. 1) THEN
!~! more than 1 node has the same max number of particles -- should resolve itself..
!~             exit
!~          else
!~! there is exactly one proc with the maximum number of particles
!~
!~             CALL MPI_ALLReduce(my_npart,tot_npart,1, &
!~                  MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!~             mean_npart = tot_npart/nprocMPI
!~
!~             IF (my_npart .EQ. min_npart) THEN
!~                tmp = my_rank
!~             ELSE
!~                tmp = 0
!~             END IF
!~             CALL MPI_ALLReduce(tmp,min_rank,1, &
!~                  MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)    
!~
!~             IF (my_npart .EQ. max_npart) THEN
!~                tmp = my_rank
!~             ELSE
!~                tmp = 0
!~             END IF
!~             CALL MPI_ALLReduce(tmp,max_rank,1, &
!~                  MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)    
!~
!~
!~             IF (my_rank .EQ. max_rank) THEN
!~                CALL MPI_SEND(particle(mean_npart+1:my_npart,0:7), &
!~                     8*(my_npart-mean_npart),MPI_DOUBLE_PRECISION,min_rank, &
!~                     1000+my_rank,MPI_COMM_WORLD,ierr)
!~                my_npart = mean_npart
!~                particle(mean_npart+1:max_npart,:) = 0.0_dp
!~             END IF
!~
!~             IF (my_rank .EQ. min_rank) THEN
!~           CALL MPI_RECV(particle(my_npart+1:my_npart+max_npart-mean_npart,:), &
!~                    8*(max_npart-mean_npart),MPI_DOUBLE_PRECISION,max_rank, &
!~                    1000+max_rank,MPI_COMM_WORLD,status,ierr)
!~                my_npart = min_npart+max_npart-mean_npart
!~             END IF
!~
!~          END IF
!~
!~       END IF
!~
!~  CALL MPI_ALLReduce(my_npart,min_npart,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
!~  CALL MPI_ALLReduce(my_npart,max_npart,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
!~       npart_ratio = real(max_npart,dp) / real(min_npart,dp)
!~
!~       loop_count = loop_count + 1
!~    end do
!~
!~ END IF
!~
!~  if ( (my_rank .eq. 0) .and. (verbosity .ge. 3) ) then
!~     print *, 'number of loops for balancing npart = ', loop_count-1
!~  end if

  END SUBROUTINE balance_npart
!
!***********************************************************************
!
  SUBROUTINE mc_neutrals(run_time_in, min_nchex)
    REAL(KIND=dp), intent(in) :: run_time_in
    REAL(KIND=dp) :: run_time, mc_time, old_timestep, skiptime, ntime, min_nchex
    REAL(KIND=dp) :: total_mass_sum, tot_mass_loss, min_dt, max_dt, ave_dt
    INTEGER :: total_n_chex, total_my_npart, ierr,loop,i
    REAL(KIND=dp) p_mc_time,p_timestep ! private copies of coresponding variables
    


!    call write_plasma_ascii
!    call write_plasma
    my_n_chex = 0

! convert times from years to seconds
    run_time = run_time_in * secs_per_year

    mc_time = 0.0_dp
! initialize source term arrays
    source_px = 0.0_dp;  source_pz = 0.0_dp
    source_m  = 0.0_dp;  source_mvsq = 0.0_dp
    source_nchex = 0

    old_timestep = timestep

    skiptime = 5.0_dp * secs_per_year
    ntime = mc_time + skiptime
    
    IF (my_rank .eq. 0)  then
      PRINT *,'begin parallel loop '
    endif
    
    !$omp parallel private(p_mc_time,p_timestep,loop)
    p_mc_time  = mc_time
    p_timestep = timestep
    loop = 0
    DO WHILE (p_mc_time .LT. run_time)
      
      IF ((my_rank .eq. 0) .and.(my_thread .eq. 0)) then
        PRINT *,'loop ',loop        
      endif

!~       IF ((mc_time .GE. ntime).and.(.false.) )THEN          ! output DIAGNOSTICS
!~          CALL MPI_Reduce(my_n_chex,total_n_chex,1, &
!~               MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!~          CALL MPI_Reduce(my_npart,total_my_npart,1, &
!~               MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!~          CALL MPI_Reduce(SUM(particle(1:my_npart,0)),total_mass_sum,1, &
!~               MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!~          CALL MPI_Reduce(my_tot_mass_loss,tot_mass_loss,1, &
!~               MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!~          CALL MPI_Reduce(my_min_dt,min_dt,1, &
!~               MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
!~          CALL MPI_Reduce(my_max_dt,max_dt,1, &
!~               MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
!~          CALL MPI_Reduce(my_ave_dt,ave_dt,1, &
!~               MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!~          IF ( (my_rank .eq. 0) .and. (verbosity .ge. 1) ) THEN
!~             IF ( verbosity .ge. 2 ) THEN
!~                PRINT 8, mc_time / secs_per_year
!~             PRINT *, '# charge exchanged this step ', total_n_chex
!~             end IF
!~             PRINT 10, total_mass_sum*macro_mass, total_my_npart
!~             IF ( verbosity .ge. 3 ) THEN
!~                PRINT 12, tot_mass_loss/timestep
!~                PRINT 14, min_dt/secs_per_year, max_dt/secs_per_year, &
!~                     ave_dt/nprocMPI/secs_per_year
!~                PRINT *
!~             end IF
!~8            FORMAT('time (yrs) = ',F13.6)
!~10           FORMAT('sum of masses = ',E13.6,' (kg), # particles = ',I14)
!~12           FORMAT('total mass loss due to photoionization (kg/sec) ',E13.6)
!~14           FORMAT('min, max, mean dt this step (yrs)',3E13.6)
!~          END IF
!~          ntime = ntime + skiptime
!~       END IF
       
      
       CALL move
       CALL split     
       
!       mc_time = mc_time + timestep
!       IF (run_time - mc_time .LT. timestep) THEN 
!          timestep = run_time - mc_time  ! adjust timestep so we finish on time
!       END IF

      ! balance every 25th step
      if (mod(loop,1).eq.0) then
       !$omp barrier
       !$omp master
                     
       CALL balance_npart()
       
       !call balance_node()
       
!       do i = 0, nthreads - 1
!        PRINT *,'  thread',i,my_particles(i)%my_npart
!       enddo
       
       !$omp end master      
       !$omp barrier
      endif

       p_mc_time = p_mc_time + p_timestep
       IF (run_time - p_mc_time .LT. p_timestep) THEN 
          p_timestep = run_time - p_mc_time  ! adjust timestep so we finish on time
       END IF
       
       loop = loop + 1

!       CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

    END DO
    timestep = old_timestep              ! put timestep back for next iteration
            
    call collect_sources_node
    !$omp end parallel
    
    IF (my_rank .eq. 0)  then
      PRINT *,'end parallel loop '
    endif
    
    !return
    
!    CALL MPI_Barrier(MPI_COMM_WORLD, ierr)    
    
    CALL collect_sources(run_time)

    CALL create_neutral_grid

    if (my_rank .eq. 0) then
! EXCLUDE outer grid since we may not be calculating out to the boundaries
       min_nchex = minval(g_source_nchex(2:ngrids,:,:))
    end if
    
    call write_neutral_grid
    call write_source

  END SUBROUTINE mc_neutrals

END MODULE neutrals_routines2_2D

