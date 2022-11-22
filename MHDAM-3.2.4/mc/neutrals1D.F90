! This code follows the trajectories of a number of mega neutral particles
! and allows the motion to be altered through charge exchange events.
!
! This is the 1D version in CARTESIAN coordinates. Particles have 3D velocity,
! but only 1D position. 
!
! This version runs as a subroutine, so that if may be plugged into CHOMBO.
! Cumulative source terms and grid values are ALWAYS calculated
!
! particle array: massi, zi, vxi, vyi, vzi, region
!  
! SI units are employed throughout this code.
!
! Here we do not use the Lipatov et al 1998 approach to counting down the
! charge-exchange time. Instead we use Prob(chex) = n_p sigma_ex ds.
!
!***********************************************************************
!***********************************************************************
!
MODULE mpi_defs_1D
  IMPLICIT NONE
  INCLUDE "mpif.h"
  INTEGER :: status(MPI_STATUS_SIZE), nproc, my_rank
END MODULE mpi_defs_1D
!
!***********************************************************************
!
MODULE global_1D
!
  USE mpi_defs_1D
!
! Global variable declarations
!
  IMPLICIT NONE
!
!  INCLUDE "mpif.h"
!
  INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND(14,307)
  INTEGER, PARAMETER, PUBLIC :: lenstr = 60
  REAL(KIND=dp), PARAMETER, PUBLIC :: au = 1.5e11_dp            ! m
  REAL(KIND=dp), PARAMETER, PUBLIC :: pi = &
       3.1415926535897932384626433832795_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: secs_per_year = 31536000.0_dp ! 365 days
  REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: particle
  REAL(KIND=dp), PUBLIC :: LISM_nH, LISM_vH, LISM_TH, k_B
  REAL(KIND=dp), PUBLIC :: LISM_v_thermal, LISM_density, boundary_layer
  REAL(KIND=dp), PUBLIC :: total_volume, neutral_mass, mega_number, macro_mass
  REAL(KIND=dp), PUBLIC :: time, timescale, timestep, my_maxchexprob
  REAL(KIND=dp), PUBLIC :: test_dens_max, test_temp_max, test_ux_max, test_uz_max
  REAL(KIND=dp), PUBLIC :: test_dens_min, test_temp_min, test_ux_min, test_uz_min
  REAL(KIND=dp), PUBLIC :: my_min_dt, my_max_dt, my_ave_dt
  INTEGER,       PUBLIC :: my_ntotal, my_npart, my_n_chex
  INTEGER,       PUBLIC :: npart, my_npart_boundary
!  LOGICAL, PUBLIC :: restart
  integer, PUBLIC :: restart, verbosity
!
! mpi defs
!  INTEGER :: status(MPI_STATUS_SIZE), nproc, my_rank
!
! plasma grid parameters
  REAL(KIND=dp), DIMENSION (:), ALLOCATABLE :: p_dens, p_temp
  REAL(KIND=dp), DIMENSION (:), ALLOCATABLE :: p_ux, p_uy, p_uz, p_region
  REAL(KIND=dp), PUBLIC :: p_zmax, p_zmin, p_dz
  INTEGER,       PUBLIC :: p_nz
!
! source grid variables
  REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: source_px, source_py
  REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: source_pz, source_mvsq
  REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: source_nchex
  REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: grid_px, grid_py, grid_pz
  REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: grid_mvsq, grid_nchex
!  REAL(KIND=dp), DIMENSION (:), ALLOCATABLE :: src_momx, src_momy
!  REAL(KIND=dp), DIMENSION (:), ALLOCATABLE :: src_momz, src_energy
!
! neutral grid parameters
  REAL(KIND=dp), DIMENSION (:), ALLOCATABLE :: grid_dens, grid_ux, grid_uy
  REAL(KIND=dp), DIMENSION (:), ALLOCATABLE :: grid_uz, grid_temp
!

END MODULE global_1D
!
!***********************************************************************
!***********************************************************************
!
MODULE neutrals_routines1_1D
!
! Lowest level subroutines
!
  USE global_1D

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: get_vr_phi, get_Maxwellian_vr, inject_maxwellian, source_chex
  PUBLIC :: plasma, exchange

CONTAINS

  SUBROUTINE get_vr_phi(vH,vT,vr,phi)
! Here we return the speed and direction (in the plasma velocity frame)
! of the new neutral after ch-ex. Follows Section 7.4 of Lipatov's book (2002)
    REAL(KIND=dp), INTENT(in)  :: vH, vT
    REAL(KIND=dp), INTENT(out) :: vr, phi
    REAL(KIND=dp), SAVE :: int_r(101,101), int_phi(101,101,101), drH, drp, dphi
    REAL(KIND=dp) :: rH, rp, int_p
    REAL(KIND=dp) :: j1, k1, rand
    INTEGER :: i, j, k, dj, dk, jlow, jhigh, klow, khigh
    LOGICAL, SAVE :: firstcall = .TRUE.

    IF (firstcall) THEN
       firstcall = .FALSE.
       
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

    END IF

    i = FLOOR(vH/vT * (100.0_dp*drH)) + 1
    i = MIN(i,101)
    i = MAX(i,1)

    CALL RANDOM_NUMBER(rand)

! use bisection to find speed
    jlow = 1
    jhigh = 101
    dj = jhigh - jlow
    j = 25
    DO WHILE (dj .GT. 1)
       j1 = int_r(i,j)
       IF (j1 .GT. rand) THEN
          jhigh = j
       ELSE
          jlow = j
       END IF
       j = FLOOR((jhigh+jlow)/2.0_dp)
       dj = jhigh - jlow
    END DO
    j = jhigh

!    vr = (j-1.5_dp)*vT*drp
    vr = (j-1.3_dp)*vT*drp
!    vr = (j-1.0_dp)*vT*drp
!    vr = (j-0.5_dp)*vT*drp
!    CALL RANDOM_NUMBER(rand)        ! try to randomize speed a bit more
!    vr = (j-rand-0.5_dp)*vT*drp

    CALL RANDOM_NUMBER(rand)

    IF (vT/vH * (100.0_dp*drH) .LT. 0.5_dp) THEN

! use ISOTROPIC speed distribution for very hot particles
       phi = rand * pi

    ELSE

! use bisection to find angle
       klow = 1
       khigh = 101
       dk = khigh - klow
       k = 50
       DO WHILE (dk .GT. 1)
          k1 = int_phi(i,j,k)
          IF (k1 .GT. rand) THEN
             khigh = k
          ELSE
             klow = k
          END IF
          k = FLOOR((khigh+klow)/2.0_dp)
          dk = khigh - klow
       END DO
       k = khigh
       phi = (k-1.5_dp)*dphi
!       phi = (k-2.0_dp)*dphi

    END IF  ! check for hot particle

! combinations tried:
!                    -1, -1  --  6300, 27.2
!                    -1.5, -2  --  6330, 27.0
!                    -1.5, -1.5  --  

  END SUBROUTINE get_vr_phi
!
!***********************************************************************
!
  SUBROUTINE get_Maxwellian_vr(vr)
! Here we randomly select a speed from a Maxwellian distribution.
! Follows Section 7.3 of Lipatov's book (2002)
    REAL(KIND=dp), INTENT(out) :: vr
    REAL(KIND=dp), SAVE :: int_Maxwellian(0:400)
    REAL(KIND=dp) :: dvr, tot_int, i1, rand
    INTEGER :: i, di, ilow, ihigh
    LOGICAL, SAVE :: firstcall = .TRUE.

    dvr = 0.01_dp
    IF (firstcall) THEN
       firstcall = .FALSE.
       
       int_Maxwellian(0) = 0.0_dp
       DO i = 1, 400
          vr = (i-0.5_dp)*dvr

          int_Maxwellian(i) = int_Maxwellian(i-1) + vr*vr * EXP(-vr*vr) * dvr

       END DO

       tot_int = int_Maxwellian(400)
       int_Maxwellian = int_Maxwellian / tot_int

    END IF

    CALL RANDOM_NUMBER(rand)

! use bisection to find speed
    ilow = 0
    ihigh = 400
    di = ihigh - ilow
    i = 100
    DO WHILE (di .GT. 1)
       i1 = int_Maxwellian(i)
       IF (i1 .GT. rand) THEN
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
  SUBROUTINE inject_maxwellian(v_thermal,my_number,zmin,zmax)
! inject a Maxwellian distributions into a prescibed subsector
    REAL(KIND=dp), INTENT (in) :: v_thermal, zmin, zmax
    INTEGER, INTENT (in)  :: my_number
    REAL(KIND=dp) :: z0, vr, theta, phi
    REAL(KIND=dp), DIMENSION(5) :: rand
    INTEGER :: i
    REAL(KIND=dp) :: vx, vy, vxy, vz
!    REAL(KIND=dp) :: xn, yn, zn, vxn, vyn, vzn
 
    DO i = 1, my_number

! unit mass for now
       particle(my_npart+i,0) = 1.0_dp

       CALL RANDOM_NUMBER(rand)
       z0 = zmin + rand(1)*(zmax-zmin)

       particle(my_npart+i,1) = z0

! Inject neutrals to have a Maxwellian velocity distribution
! This method comes from "The hybrid multiscalesimulation technology"
! by A.S. Lipatov, section 7.2.2.

! Based on 3D speed distribution
       CALL get_Maxwellian_vr(vr) ! obtain Maxwellian vr
       vr = vr * v_thermal

! must AVOID making vx=vy=0
       theta = ACOS(0.999999_dp - 1.99999_dp*rand(2))
       vz = vr * COS(theta)
       vxy = SQRT(vr*vr - vz*vz)

       phi = 2.0_dp*pi*rand(3)
       vx = vxy * COS(phi)
       vy = vxy * SIN(phi)

! update particle
       particle(my_npart+i,2) = vx
       particle(my_npart+i,3) = vy
       particle(my_npart+i,4) = vz + LISM_vH

! Maxwellian particles are region 1 (LISM)
       particle(my_npart+i,5) = 1.1_dp

    END DO

    my_npart = my_npart + my_number

  END SUBROUTINE inject_maxwellian
!
!***********************************************************************
!
  SUBROUTINE plasma(zi_in,n_p,T_p,ux_p,uy_p,uz_p,region_p)
    REAL(KIND=dp), INTENT(in)  :: zi_in
    REAL(KIND=dp), INTENT(out) :: n_p, ux_p, uy_p, uz_p, T_p, region_p
    REAL(KIND=dp) :: zi
    INTEGER :: ii
! return plasma number density, velocity components and temp at a given point

    zi = MIN(zi_in,0.999_dp*p_zmax)
    zi = MAX(zi,1.001_dp*p_zmin)

    ii = FLOOR((zi-p_zmin)/p_dz) + 1

    n_p = p_dens(ii)
    T_p = p_temp(ii)
    ux_p = p_ux(ii)
    uy_p = p_uy(ii)
    uz_p = p_uz(ii)
    region_p = p_region(ii)

  END SUBROUTINE plasma
!
!***********************************************************************
!
  SUBROUTINE exchange(massi,zi,vx,vy,vz,region,vxn,vyn,vzn)
    REAL(KIND=dp), INTENT (in)  :: massi, zi, vx, vy, vz, region
    REAL(KIND=dp), INTENT (out) :: vxn, vyn, vzn
    REAL(KIND=dp) :: rand, n_p, T_p, ux, uy, uz, region_p, v_p_th
    REAL(KIND=dp) :: vv1, vv2, vv3, vv, vvsq, ww1, ww2, ww3, pp1, pp2, p, vdotw
    REAL(KIND=dp) :: rotvv1, rotvv2, rotvv3, wprojv1, wprojv2, wprojv3
    REAL(KIND=dp) :: r1_1, r1_2, r1_3, rr, r2_1, r2_2, r2_3
    REAL(KIND=dp) :: tp, cp, sp, vvnew1, vvnew2, vvnew3, omega, vr, phi

    CALL plasma(zi,n_p,T_p,ux,uy,uz,region_p)
    v_p_th = SQRT(2.0_dp*k_B*T_p/neutral_mass)

! NOTE: MASS of new particle same as OLD
    vv1 = vx-ux
    vv2 = vy-uy
    vv3 = vz-uz
    vvsq = vv1*vv1+vv2*vv2+vv3*vv3
    vv = SQRT(vvsq)
! include a selection effect due to thermal protons
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

    CALL RANDOM_NUMBER(rand)
    omega = 2.0_dp*pi * rand             ! random phase for scatter by phi

    vvnew1 = wprojv1 + r1_1*COS(omega) + r2_1*SIN(omega)
    vvnew2 = wprojv2 + r1_2*COS(omega) + r2_2*SIN(omega)
    vvnew3 = wprojv3 + r1_3*COS(omega) + r2_3*SIN(omega)

    vxn = vvnew1 + ux
    vyn = vvnew2 + uy
    vzn = vvnew3 + uz

! add ch-ex event to plasma source term
    CALL source_chex(massi,zi,vx,vy,vz,vxn,vyn,vzn,FLOOR(region))

  END SUBROUTINE exchange
!
!***********************************************************************
!
  SUBROUTINE source_chex(mass,zi,vxold,vyold,vzold,vxnew,vynew,vznew,region)
    REAL(KIND=dp), INTENT(in)  :: mass, zi, vxold, vzold, vyold
    REAL(KIND=dp), INTENT(in)  :: vxnew, vznew, vynew
    INTEGER, INTENT(in) ::  region
    REAL(KIND=dp) :: dvx, dvy, dvz, dvsq
    INTEGER :: ii

! add ch-ex source to CARTESIAN grid

! MAKE SURE that p_zmin < zi < p_zmax
    IF ( (zi .GT. p_zmin) .AND. (zi .LT. p_zmax) ) THEN
       
       ii = FLOOR((zi-p_zmin)/p_dz) + 1

       dvx  = vxold - vxnew
       dvy  = vyold - vynew
       dvz  = vzold - vznew
       dvsq = (vxold*vxold + vzold*vzold + vyold*vyold &
             - vxnew*vxnew - vznew*vznew - vynew*vynew)

! no mass source since we assume m_p = m_H
! momentum and energy source terms pz, m*u^2
       source_px(region,ii)   = source_px(region,ii) + dvx * mass
       source_py(region,ii)   = source_py(region,ii) + dvy * mass
       source_pz(region,ii)   = source_pz(region,ii) + dvz * mass
       source_mvsq(region,ii) = source_mvsq(region,ii) + dvsq * mass

! count number of ch-ex per grid point as an ACCURACY DIAGNOSTIC
       source_nchex(region,ii) = source_nchex(region,ii) + mass

       my_n_chex = my_n_chex + 1

    END IF

  END SUBROUTINE source_chex

END MODULE neutrals_routines1_1D
!
!***********************************************************************
!***********************************************************************
!
MODULE neutrals_routines2_1D
!
! Second level subroutines
!
#ifdef F_ABSOFT
!DIR$ NAME (release_cache="_f90a_free_all")
#endif
!
  USE global_1D
  USE neutrals_routines1_1D

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: mc_init, move, collect_source, collect_grid
  PUBLIC :: output_source, output_grid, output_raw, load_raw
  PUBLIC :: mc_neutrals

CONTAINS

  SUBROUTINE mc_init(p_nz_in, p_zmin_in, p_zmax_in, restart_in, loadfile, &
               total_neutrals, lism_nh_in, lism_vh_in, lism_th_in, verbosity_in)
    REAL(KIND=dp) :: lism_nh_in, lism_vh_in, lism_th_in
    REAL(KIND=dp) :: p_zmin_in, p_zmax_in, achieved_density
    INTEGER, DIMENSION(:), ALLOCATABLE :: put_seq
    INTEGER :: n, ierr, p_nz_in, total_neutrals
    CHARACTER (LEN = lenstr) :: loadfile
!    LOGICAL :: restart_in
    integer :: restart_in, verbosity_in

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)    ! nproc = # of processors

    p_nz = p_nz_in
    p_zmin = p_zmin_in
    p_zmax = p_zmax_in
    restart = restart_in
    verbosity = verbosity_in

! Set random number seed to be different for each processor
    CALL RANDOM_SEED(size=n)
    ALLOCATE(put_seq(n))
    put_seq = my_rank * my_rank
    CALL RANDOM_seed(put=put_seq)

! set up global parameters

! LISM parameters
    LISM_nH = lism_nh_in               ! per cubic meter
    LISM_vH = lism_vh_in               ! meters per second
    LISM_TH = lism_th_in               ! degrees Kelvin

    k_B = 1.38e-23_dp                  ! Joules per Kelvin (m^2 kg s^{-2} K^{-1})
    neutral_mass = 1.67e-27_dp         ! Kilograms
    LISM_v_thermal = SQRT(2.0_dp*k_B*LISM_TH/neutral_mass) ! 3D thermal speed

    ALLOCATE(source_px(1:2,1:p_nz),source_py(1:2,1:p_nz))
    ALLOCATE(source_pz(1:2,1:p_nz),source_mvsq(1:2,1:p_nz),source_nchex(1:2,1:p_nz))
    source_px = 0.0_dp;  source_py = 0.0_dp
    source_pz = 0.0_dp;  source_mvsq = 0.0_dp;  source_nchex = 0.0_dp
    ALLOCATE(grid_px(1:2,1:p_nz), grid_py(1:2,1:p_nz))
    ALLOCATE(grid_pz(1:2,1:p_nz), grid_mvsq(1:2,1:p_nz), grid_nchex(1:2,1:p_nz))
    grid_px = 0.0_dp;  grid_py = 0.0_dp;  grid_pz = 0.0_dp
    grid_mvsq = 0.0_dp;  grid_nchex = 0
    ALLOCATE(grid_dens(p_nz),grid_ux(p_nz),grid_uy(p_nz))
    allocate(grid_uz(p_nz),grid_temp(p_nz))
    grid_dens = 0.0_dp;  grid_ux = 0.0_dp; grid_uy = 0.0_dp
    grid_uz = 0.0_dp;  grid_temp = 0.0_dp
    ALLOCATE(p_dens(1:p_nz), p_ux(1:p_nz), p_uy(1:p_nz), p_uz(1:p_nz))
    allocate(p_temp(1:p_nz), p_region(1:p_nz))
!    ALLOCATE(src_momx(1:p_nz),src_momy(1:p_nz),src_momz(1:p_nz),src_energy(1:p_nz))

    timescale = (p_zmax-p_zmin)/ABS(LISM_vH)
    time = 0.0_dp
! set up coarse timestep which moves LISM atoms about 1.5% of a mfp
!    timestep = secs_per_year / 2.0_dp    ! half a year
    timestep = secs_per_year             ! one year

! for true 1D we assume that dx=dy=1 m
    total_volume = p_zmax - p_zmin
    LISM_density = total_neutrals/total_volume

! choose a boundary layer thickness so that 99.9999998% of particles with
! velocity directed inward which make it into the domain come from this layer
    boundary_layer = 6.0_dp*(LISM_v_thermal+ABS(LISM_vH))*timestep

    my_npart_boundary = FLOOR( LISM_density*boundary_layer/nproc + 0.5_dp)

    achieved_density = nproc*my_npart_boundary / boundary_layer
    npart = FLOOR(achieved_density*total_volume)

! NOTE: must define macro_mass in terms of EXPECTED # of particles (npart)
    mega_number = LISM_nH*total_volume/npart ! # neutrals in megaparticle
    macro_mass = neutral_mass * mega_number           ! mass of macroparticle
    LISM_density = npart/total_volume

    my_ntotal = npart/nproc + nproc       ! allow for rounding
! create extra space for random fluctuations in my_npart -- Poisson process
    my_ntotal = my_ntotal + &
         FLOOR(MAX(0.2_dp*my_ntotal,200_dp*SQRT(1.0_dp*my_ntotal))) + 1
! particle array: mass (degraded for photoionization), x, z, vx, vz, region

    my_ntotal = 3*my_ntotal

! mass, z, vx, vy, vz, region --- MAX_MASS removed
    ALLOCATE(particle(my_ntotal,0:5))
    particle = 0.0_dp
    my_npart = 0

    IF ( restart .eq. 1 ) THEN
! load particle data from existing file
       CALL load_raw(loadfile)
       CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
    ELSE
! set up initial distribution 
       CALL inject_maxwellian(LISM_v_thermal, &
            nint(npart/nproc*120.0_dp*au/p_zmax),0.0_dp*au,120.0_dp*au)
    END IF

9 FORMAT(6e13.4)

  END SUBROUTINE mc_init
!
!***********************************************************************
!
  SUBROUTINE move
    REAL(KIND=dp), PARAMETER :: ideal_chex_prob = 0.01_dp
    REAL(KIND=dp) :: massi, zi, vx, vz, vy, region
    REAL(KIND=dp) :: n_p, T_p, ux, uy, uz, region_p
    REAL(KIND=dp) :: vxn, vzn, vyn
    REAL(KIND=dp) :: v_p_th, delta_u, omega_p
    REAL(KIND=dp) :: v_rel_p, sigma_ex, subtime, dt, ideal_dt, ds
    REAL(KIND=dp) :: n_p2, T_p2, ux2, uy2, uz2, region_p2, v_rel_p2, ideal_dt2
    REAL(KIND=dp) :: delta_u2, rand, dz
    REAL(KIND=dp) :: derf  ! need to link in error function
    INTEGER :: my_old_npart, i, count
    
    my_min_dt = timestep
    my_max_dt = 0.0_dp
    my_ave_dt = 0.0_dp
    count = 0

    my_old_npart = my_npart
    i = 0
    DO WHILE (i .LT. my_npart)
       i = i + 1

       massi = particle(i,0)   
       zi  = particle(i,1)
       vx = particle(i,2)
       vy = particle(i,3)
       vz = particle(i,4)
       region = particle(i,5)

       subtime = 0.0_dp
! use adaptive timestep
       DO WHILE (subtime .LT. timestep)

          CALL plasma(zi,n_p,T_p,ux,uy,uz,region_p)
          v_p_th = SQRT(2.0_dp*k_B*T_p/neutral_mass)

! from Ripken and Fahr 1983, notation follows Lipatov, Zank & Pauls 1998
          delta_u = SQRT( (vx-ux)**2 + (vy-uy)**2 + (vz-uz)**2)
          omega_p = delta_u/v_p_th
          v_rel_p = v_p_th * ( EXP(-omega_p*omega_p)/SQRT(pi) + &
                                   (omega_p + 0.5_dp/omega_p)*derf(omega_p) )
   sigma_ex = ( 1.64e-7_dp - 6.95e-9_dp*LOG(v_rel_p/1e-2_dp) )**2 * 1e-4_dp !m^2

          ideal_dt = ideal_chex_prob/(n_p*sigma_ex*v_rel_p)
! NOTE: might need to check that dt is small enough over entire space step
! (seems ok, even at 1 AU)

          dt = MIN(ideal_dt,timestep)


! check FUTURE dt
          CALL plasma(zi+vz*dt,n_p2,T_p2,ux2,uy2,uz2,region_p2)
          delta_u2 = SQRT( (vx-ux2)**2 + (vy-uy2)**2 + (vz-uz2)**2)
          v_rel_p2 = v_rel_p * SQRT(T_p2/T_p) * delta_u2/delta_u  ! ESTIMATE

          ideal_dt2 = ideal_chex_prob/(n_p2*sigma_ex*v_rel_p2)

          IF (ideal_dt/ideal_dt2 .GT. 2.0_dp) THEN
             dt = 0.01_dp*ideal_dt2
          END IF

! TESTING:
          my_max_dt = MAX(my_max_dt,dt)
          my_min_dt = MIN(my_min_dt,dt)

          IF (subtime+dt .GT. timestep) THEN
             dt = timestep - subtime
          END IF

          my_ave_dt = my_ave_dt + dt
          count = count + 1

          subtime = subtime + dt

          ds = v_rel_p * dt  ! use (RELATIVE) distance travelled for Prob(chex)
          CALL RANDOM_NUMBER(rand)
          IF ( rand .LT. n_p*sigma_ex*ds ) THEN   ! charge exchange occurs

             CALL RANDOM_NUMBER(rand)
             dz = rand * dt * vz

             CALL exchange(massi,zi+dz,vx,vy,vz,region,vxn,vyn,vzn)

             vx = vxn;  vy = vyn;  vz = vzn
             region = region_p

          END IF

! advance particles in a straight line at current velocity (Euler's method)
          zi = zi + vz*dt

       END DO

! update particle specs
       particle(i,0) = massi
       particle(i,1) = zi
       particle(i,2) = vx
       particle(i,3) = vy
       particle(i,4) = vz
       particle(i,5) = region

! remove particles that have moved out of the domain
       IF ( (zi .GT. p_zmax) .OR. (zi .LT. p_zmin) ) THEN
          particle(i,:) = particle(my_npart,:)   ! fill gap with last particle
          my_npart = my_npart - 1    ! stop tracking particles outside domain
          i = i - 1                  ! don't increment if we remove a particle
       END IF

    END DO

    my_ave_dt = my_ave_dt / count

    particle(my_npart+1:my_old_npart,:) = 0.0_dp ! zero removed particles

! INJECT new particles at inflow boundary 
    if (my_npart+my_npart_boundary .gt. my_ntotal) then
       print *, 'maximum number of particles exceeded by proc ', my_rank
    else
       CALL inject_maxwellian(LISM_v_thermal,my_npart_boundary,p_zmin-boundary_layer,p_zmin)
    end if

  END SUBROUTINE move
!
!***********************************************************************
!
  subroutine collect_grid
    REAL(KIND=dp), DIMENSION (:), ALLOCATABLE :: H_dens, H_ux, H_uy, H_uz, H_temp
    REAL(KIND=dp) :: massi, zi, vx, vy, vz, ux, uy, uz
    REAL(KIND=dp) :: vrel_sq
    INTEGER :: i, iz, ierr
    
! create grid on each node
    ALLOCATE(H_dens(p_nz),H_ux(p_nz),H_uy(p_nz),H_uz(p_nz),H_temp(p_nz))

    H_dens = 0.0_dp;  H_temp = 0.0_dp
    H_ux = 0.0_dp;  H_uz = 0.0_dp;  H_uz = 0.0_dp

    DO i = 1, my_npart

       zi = particle(i,1)

! don't include boundary points
       IF ( (zi .LT. p_zmax) .AND. (zi .GT. p_zmin) ) THEN

          massi = particle(i,0)
          vx = particle(i,2)
          vy = particle(i,3)
          vz = particle(i,4)
          
          iz = FLOOR((zi-p_zmin)/p_dz) + 1
                
! distribute number density to adjacent nodes
          H_dens(iz) = H_dens(iz) + massi
! calculate bulk velocity
          H_ux(iz)   = H_ux(iz) + vx * massi
          H_uy(iz)   = H_uy(iz) + vy * massi
          H_uz(iz)   = H_uz(iz) + vz * massi

       END IF
 
    END DO

    grid_ux = 0.0_dp;  grid_uy = 0.0_dp;  grid_uz = 0.0_dp; 
    grid_dens = 0.0_dp; grid_temp = 0.0_dp

! add together all density data
    CALL MPI_Allreduce(H_dens, grid_dens, p_nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
! add together all velocity data
    CALL MPI_Allreduce(H_ux, grid_ux, p_nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(H_uy, grid_uy, p_nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(H_uz, grid_uz, p_nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    WHERE (grid_dens /= 0.0_dp)
       grid_ux = grid_ux / grid_dens
       grid_uy = grid_uy / grid_dens
       grid_uz = grid_uz / grid_dens
    END WHERE

! calculate temperature
    H_temp = 0.0_dp
    DO i = 1, my_npart

       zi = particle(i,1)

! don't include boundary points
       IF ( (zi .LT. p_zmax) .AND. (zi .GT. p_zmin) ) THEN

          massi = particle(i,0)
          vx = particle(i,2)
          vy = particle(i,3)
          vz = particle(i,4)

          iz = FLOOR((zi-p_zmin)/p_dz) + 1

          ux = grid_ux(iz)
          uy = grid_uy(iz)
          uz = grid_uz(iz)

! calculate temperature
          vrel_sq = (vx-ux)**2 + (vy-uy)**2 + (vz-uz)**2
          H_temp(iz) = H_temp(iz) + vrel_sq * massi

       END IF

    END DO

! add together all vrel_sq data
    CALL MPI_reduce(H_temp, grid_temp, p_nz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (my_rank .eq. 0) then

       WHERE (grid_dens /= 0.0_dp)
! H temp must be normalized
          grid_temp = grid_temp * neutral_mass / (3.0_dp*k_B*grid_dens)
       END WHERE

! normalise grid density
       grid_dens = grid_dens * mega_number / p_dz ! per m^3

    end if
    
    DEALLOCATE(H_dens, H_temp, H_ux, H_uy, H_uz)
#ifdef F_ABSOFT
    CALL release_cache()
#endif

  end subroutine collect_grid
!
!***********************************************************************
!
  SUBROUTINE output_grid(gridfile)
    INTEGER :: kk, ierr
    CHARACTER (LEN=60) :: gridfile

    IF (my_rank .EQ. 0) THEN
!       OPEN(unit=1,file=trim(gridfile),form='formatted',status='unknown',action='write')
       OPEN(unit=1,file='output.grid',form='formatted',status='unknown',action='write')
       WRITE(1,*) p_nz
       WRITE(1,*) p_zmin/au
       WRITE(1,*) p_zmax/au
       DO kk = 1, p_nz
          WRITE(1,88) (p_zmin + (kk-0.5_dp)*p_dz)/au, &
               grid_dens(kk), grid_ux(kk), grid_uy(kk), grid_uz(kk), grid_temp(kk)
       END DO

       CLOSE(1)

       if (verbosity .ge. 3) then
          print *
          print *, 'average grid density ', sum(grid_dens)/p_nz
          print *, 'average grid temperature ', sum(grid_temp)/p_nz
       end if

    END IF

    CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

88 FORMAT(6e15.5)

  END SUBROUTINE output_grid
!
!***********************************************************************
!
  SUBROUTINE output_raw(rawfile)
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: rank, npart, total_npart, ierr
    CHARACTER (LEN=60) :: rawfile

    total_npart = 0
  CALL MPI_Reduce(my_npart,total_npart,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    IF (my_rank .EQ. 0) THEN

       OPEN(unit=1,file=rawfile,form='unformatted',status='unknown',action='write')
       WRITE(1) total_npart
       WRITE(1) p_zmin
       WRITE(1) p_zmax
       WRITE(1) LISM_v_thermal
       WRITE(1) time
       WRITE(1) nproc
       WRITE(1) 1.0_dp*my_npart + 0.1_dp
       WRITE(1) particle(1:my_npart,0:5)

       DO rank = 1, nproc-1
          CALL MPI_RECV(npart, 1, MPI_INTEGER, rank, &
                           1000+rank, MPI_COMM_WORLD, status, ierr)
          CALL MPI_RECV(particle(1:npart,0:5), 6*npart, MPI_DOUBLE_PRECISION, &
                           rank, 2000+rank, MPI_COMM_WORLD, status, ierr)

          WRITE(1) 1.0_dp*npart + 0.1_dp
          WRITE(1) particle(1:npart,0:5)

#ifdef F_ABSOFT
          CALL release_cache()
#endif

       END DO
       CLOSE(1)

    ELSE

       DO rank = 1, nproc-1

          IF (my_rank .EQ. rank) THEN
             CALL MPI_SEND(my_npart, 1, MPI_INTEGER, 0, 1000+rank, &
                           MPI_COMM_WORLD, ierr)
             CALL MPI_SEND(particle(1:my_npart,0:5), 6*my_npart, &
                    MPI_DOUBLE_PRECISION, 0, 2000+rank, MPI_COMM_WORLD, ierr)
          END IF
       END DO

    END IF

9 FORMAT(6e13.4)

  END SUBROUTINE output_raw
!
!***********************************************************************
!
  SUBROUTINE load_raw(loadfile)
    REAL(KIND=dp) :: load_p_zmin, load_p_zmax, load_LISM_v_thermal
    REAL(KIND=dp) :: real_my_npart
    REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: temp
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: load_total_npart, load_nproc
    INTEGER :: blocksize, new_npart, read_npart
    INTEGER :: rank, ierr
    CHARACTER (LEN = lenstr) :: loadfile

    IF (my_rank .EQ. 0) THEN
       OPEN(unit=1,file=loadfile,form='unformatted',status='old',action='read')
       READ(1) load_total_npart

       if (verbosity .ge. 3) then
          PRINT *,'loading ',load_total_npart,' particles'
       end if
       IF (load_total_npart .GT. my_ntotal*nproc) THEN
          PRINT *,'too many particles, STOPPING'
          STOP
       END IF

       READ(1) load_p_zmin
       READ(1) load_p_zmax
       READ(1) load_LISM_v_thermal
       READ(1) time
       READ(1) load_nproc

       blocksize = 2*MAX(load_total_npart/load_nproc,my_ntotal)
       ALLOCATE(temp(1:blocksize,0:5))

       new_npart = load_total_npart/nproc   ! how many each proc will receive
       read_npart = 0
       DO rank = 1, nproc-1

          DO WHILE (read_npart .LT. new_npart)
             READ(1) real_my_npart
             my_npart = FLOOR(real_my_npart)
             if (verbosity .ge. 3) then
                PRINT *, 'reading ', my_npart, ' particles'
             end if
             READ(1) temp(read_npart+1:read_npart+my_npart,0:5)
             read_npart = read_npart + my_npart
          END DO
 
          CALL MPI_SEND(new_npart, 1, MPI_INTEGER, rank, 1000+rank, &
               MPI_COMM_WORLD, ierr)
          if (verbosity .ge. 3) then
             PRINT *, 'sending ',new_npart,' particles to proc # ',rank
          end if
          CALL MPI_SEND(temp(1:new_npart,0:5), 6*new_npart, &
               MPI_DOUBLE_PRECISION, rank, 2000+rank, MPI_COMM_WORLD, ierr)
 
          temp(1:read_npart-new_npart,0:5) = temp(new_npart+1:read_npart,0:5)
          read_npart = read_npart - new_npart

#ifdef F_ABSOFT
          CALL release_cache()
#endif

       END DO

! master takes particles left over from rounding
       new_npart = new_npart + load_total_npart - nproc*new_npart
       DO WHILE (read_npart .LT. new_npart)
          READ(1) real_my_npart
          my_npart = FLOOR(real_my_npart)
         
          READ(1) temp(read_npart+1:read_npart+my_npart,0:5)
          read_npart = read_npart + my_npart
       END DO

       if (verbosity .ge. 3) then
          PRINT *, 'keeping remaining ', new_npart, ' particles for myself'
       end if
       my_npart = new_npart
       particle(1:my_npart,0:5) = temp(1:new_npart,0:5)

       CLOSE(1)
       DEALLOCATE(temp)

#ifdef F_ABSOFT
       CALL release_cache()
#endif

    ELSE

       DO rank = 1, nproc-1

          IF (my_rank .EQ. rank) THEN
             CALL MPI_RECV(my_npart, 1, MPI_INTEGER, 0, 1000+rank, &
                  MPI_COMM_WORLD, status, ierr)
             CALL MPI_RECV(particle(1:my_npart,0:5), 6*my_npart, &
                  MPI_DOUBLE_PRECISION, 0, 2000+rank, MPI_COMM_WORLD, status, ierr)
             
          END IF
       END DO

    END IF

    CALL MPI_BCAST(time, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    particle(my_npart+1:my_ntotal,:) = 0.0_dp

  END SUBROUTINE load_raw
!
!***********************************************************************
!
  SUBROUTINE collect_source(run_time)
    REAL(KIND=dp) :: run_time
    INTEGER :: ierr

! add together all source data
    CALL MPI_Reduce(source_px,    grid_px,  2*p_nz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce(source_py,    grid_py,  2*p_nz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce(source_pz,    grid_pz,  2*p_nz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce(source_mvsq, grid_mvsq, 2*p_nz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce(source_nchex, grid_nchex, 2*p_nz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (my_rank .eq. 0) then

! normalize with respect to local VOLUME (dx=dy=1 m)
       grid_px = grid_px / p_dz
       grid_py = grid_py / p_dz
       grid_pz = grid_pz / p_dz
       grid_mvsq = grid_mvsq / p_dz

! normalize with respect to time and set units
       grid_px   = macro_mass * grid_px   / run_time  ! kg m/s /m^3 /s
       grid_py   = macro_mass * grid_py   / run_time  ! kg m/s /m^3 /s
       grid_pz   = macro_mass * grid_pz   / run_time  ! kg m/s /m^3 /s
       grid_mvsq = macro_mass * grid_mvsq / run_time  ! kg m^2/s^2 /m^3 /s

    end if

  END SUBROUTINE collect_source
!
!***********************************************************************
!
  SUBROUTINE output_source(sourcefile)
    INTEGER :: i, reg
    CHARACTER (LEN=60) :: sourcefile, totfile
    CHARACTER (LEN=6) :: suffix

    suffix = '_r1_r2'

    IF (my_rank .EQ. 0) THEN

       DO reg = 1, 2

          if (reg .eq. 1) then
             totfile = 'output.source_r1'
          else
             totfile = 'output.source_r2'
          end if

          OPEN(unit=1,file=totfile,status='unknown',action='write')
          WRITE(1,*) p_nz
          WRITE(1,*) p_zmin/au
          WRITE(1,*) p_zmax/au
          DO i = 1, p_nz

! remember to take half mvsq for KINETIC ENERGY
             WRITE(1,FMT='(1F9.3,5E16.6)') (p_zmin+(i-0.5_dp)*p_dz)/au, &
                  grid_px(reg,i), grid_py(reg,i), &
                  grid_pz(reg,i), grid_mvsq(reg,i)/2.0_dp, grid_nchex(reg,i)
          END DO
          CLOSE(1)

       END DO   ! reg loop

! write total source to data file
!       OPEN(unit=1,file=TRIM(sourcefile)//'_tot',status='unknown',action='write')
       OPEN(unit=1,file='output.source_tot',status='unknown',action='write')
       WRITE(1,*) p_nz
       WRITE(1,*) p_zmin/au
       WRITE(1,*) p_zmax/au
       DO i = 1, p_nz

! remember to take half mvsq for KINETIC ENERGY
          WRITE(1,FMT='(1F9.3,5E16.6)') (p_zmin+(i-0.5_dp)*p_dz)/au, &
               SUM(grid_px(:,i)), SUM(grid_py(:,i)), &
               SUM(grid_pz(:,i)), SUM(grid_mvsq(:,i))/2.0_dp, SUM(grid_nchex(:,i))
       END DO
       CLOSE(1)

    END IF

  END SUBROUTINE output_source
!
!***********************************************************************
!
  SUBROUTINE mc_neutrals(run_time, min_nchex, p_dens_in, p_ux_in, &
       p_uy_in, p_uz_in, p_temp_in, p_region_in, &
       src_momx, src_momy, src_momz, src_energy, &
       grid_dens1, grid_ux1, grid_uy1, grid_uz1, grid_temp1)
    REAL(KIND=dp), DIMENSION (:) :: p_dens_in, p_temp_in, p_ux_in
    REAL(KIND=dp), DIMENSION (:) :: p_uy_in, p_uz_in, p_region_in
    REAL(KIND=dp), DIMENSION (:) :: src_momx, src_momy, src_momz, src_energy
    REAL(KIND=dp), DIMENSION (:) :: grid_dens1, grid_ux1, grid_uy1
    REAL(KIND=dp), dimension (:) :: grid_uz1, grid_temp1
    REAL(KIND=dp) :: run_time, mc_time, old_timestep, min_nchex
    INTEGER :: ii, ierr

    p_dens = p_dens_in
    p_temp = p_temp_in
    p_ux = p_ux_in
    p_uy = p_uy_in
    p_uz = p_uz_in
    p_region = p_region_in

! BROADCAST plasma data
    CALL MPI_BCAST(p_dens, p_nz, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(p_ux, p_nz, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(p_uy, p_nz, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(p_uz, p_nz, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(p_temp, p_nz, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(p_region, p_nz, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#ifdef F_ABSOFT
    CALL release_cache()
#endif

    p_dz = (p_zmax-p_zmin)/p_nz

    my_n_chex = 0

! convert times from years to seconds
    run_time = run_time * secs_per_year

    mc_time = 0.0_dp
! initialize source term arrays
    source_px = 0.0_dp
    source_py = 0.0_dp
    source_pz = 0.0_dp
    source_mvsq = 0.0_dp
    source_nchex = 0.0_dp
    old_timestep = timestep

    DO WHILE (mc_time .LT. run_time)

       CALL move

       mc_time = mc_time + timestep
       IF (run_time - mc_time .LT. timestep) THEN 
          timestep = run_time - mc_time  ! adjust timestep so we finish on time
       END IF

       CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

    END DO

    timestep = old_timestep              ! put timestep back for next iteration

    CALL collect_source(run_time)
    min_nchex = SUM(grid_nchex)
    DO ii = 1, p_nz
       min_nchex = MIN(min_nchex,SUM(grid_nchex(:,ii)))
    END DO

    src_momx = grid_px(1,:) + grid_px(2,:)
    src_momy = grid_py(1,:) + grid_py(2,:)
    src_momz = grid_pz(1,:) + grid_pz(2,:)
    src_energy = 0.5_dp * ( grid_mvsq(1,:) + grid_mvsq(2,:) )

    call collect_grid
    grid_dens1 = grid_dens
    grid_ux1 = grid_ux
    grid_uy1 = grid_uy
    grid_uz1 = grid_uz
    grid_temp1 = grid_temp

!!!! source term diagnostics!!!!
    if ( (my_rank .eq. 0) .and. (verbosity .ge. 3) ) then
       print *
       print *, 'min nchex', min_nchex
       print *, 'min energy source', minval(src_energy)
       print *, 'max energy source', maxval(src_energy)
       print *, 'mean energy source', sum(src_energy)/p_nz
    end if

  END SUBROUTINE mc_neutrals

END MODULE neutrals_routines2_1D
