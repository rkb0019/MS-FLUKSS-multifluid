! This code follows the trajectories of 9a number of mega neutral particles
! and allows the motion to be altered through charge exchange events.
!
! This version uses runs on a Cartesian plasma grid for long heliotails.
! Sources are collected on nested grids with different levels of refinement.
!  
! SI units are employed throughout this code.
!
!***********************************************************************
!***********************************************************************
!

!pgi$g noconcur

MODULE global
!
! Global variable declarations
!
  IMPLICIT NONE
!
  INCLUDE 'mpif.h'
!
  
  
  INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND(14,307)
  INTEGER, PARAMETER, PUBLIC :: sp = SELECTED_REAL_KIND(6)
  
  INTEGER, PARAMETER, PUBLIC :: WCOMP = 6 ! number of components stored in g_plasma
  INTEGER, PARAMETER, PUBLIC :: WRHO  = 0
  INTEGER, PARAMETER, PUBLIC :: WVELX = 1
  INTEGER, PARAMETER, PUBLIC :: WVELY = 2
  INTEGER, PARAMETER, PUBLIC :: WVELZ = 3
  INTEGER, PARAMETER, PUBLIC :: WTEMP = 4
  INTEGER, PARAMETER, PUBLIC :: WREG  = 5
  
  
  INTEGER, PARAMETER, PUBLIC :: URHO  = 0
  INTEGER, PARAMETER, PUBLIC :: UMOMX = 1
  INTEGER, PARAMETER, PUBLIC :: UMOMY = 2
  INTEGER, PARAMETER, PUBLIC :: UMOMZ = 3
  INTEGER, PARAMETER, PUBLIC :: UENG  = 4
  
  INTEGER, PARAMETER, PUBLIC :: PNUM  = 9 ! number of components stored for each particle
  
  
  type kinetic_level
    REAL(KIND=dp) :: p_xmin, p_xmax, p_ymin, p_ymax, p_zmin, p_zmax
    REAL(KIND=dp) :: sizex , sizey , sizez
    REAL(KIND=dp) :: p_dx, p_dy, p_dz
    INTEGER       :: p_nx, p_ny, p_nz
    
    INTEGER       :: m_ref_ratio ! ref ratio from previous level. 1 for the base level.
        
    REAL(KIND=dp), DIMENSION (:,:,:),   ALLOCATABLE :: source_m,   source_px,     source_py,   source_pz,   source_mvsq    
    REAL(KIND=dp), DIMENSION (:,:,:),   POINTER     :: g_source_m, g_source_px, g_source_py, g_source_pz, g_source_mvsq    
    REAL(KIND=dp), DIMENSION (:,:,:,:), POINTER     :: g_plasma
    
    ! real becuase we need to store values in FArrayBox later
    REAL(KIND=dp), DIMENSION (:,:,:), ALLOCATABLE   :: source_nchex
    REAL(KIND=dp), DIMENSION (:,:,:), POINTER       :: g_source_nchex
    
    ! neutrals
    REAL(KIND=dp), DIMENSION (:,:,:),   POINTER         :: g_H_dens, g_H_ux, g_H_uy, g_H_uz, g_H_temp
    REAL(KIND=dp), DIMENSION (:,:,:),   ALLOCATABLE     :: H_dens, H_ux, H_uy, H_uz, H_temp
    
    INTEGER, PUBLIC :: num_samples
    INTEGER, PUBLIC :: sampling_dt  ! how frequently should be sampling for the finest level
    
    
    
  end type kinetic_level
  
  
   type charge_exchange_event
    INTEGER :: level
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    REAL(KIND=dp) :: source_m, source_px, source_py, source_pz, source_mvsq    
  end type charge_exchange_event
  
  type particle_storage
    INTEGER :: my_npart ! number of particles actually stored
    REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: particles__
  end type particle_storage
  
  type node_workload
    INTEGER :: procID
    INTEGER :: part_disbalance
  end type node_workload
  
  
  INTEGER, PARAMETER, PUBLIC :: lenstr = 60
  REAL(KIND=dp), PARAMETER, PUBLIC :: au = 1.5e11_dp            ! m
  REAL(KIND=dp), PARAMETER, PUBLIC :: r_E_sq = au*au            ! m^2
  REAL(KIND=dp), PARAMETER, PUBLIC :: beta_ph_E = 8e-8_dp       ! s^{-1}
  REAL(KIND=dp), PARAMETER, PUBLIC :: phi_rmax_sq = (0.0_dp*au)**2 ! photoionization source term outer boundary
  REAL(KIND=dp), PARAMETER, PUBLIC :: pi = &
       3.1415926535897932384626433832795_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: secs_per_year = 31536000.0_dp  
  
  REAL(KIND=dp), PUBLIC :: LISM_nH, LISM_vx, LISM_vy, LISM_vz, LISM_TH, k_B
  REAL(KIND=dp), PUBLIC :: LISM_vH, LISM_v_thermal, LISM_density
  REAL(KIND=dp), PUBLIC :: boundary_thickness, total_volume
  REAL(KIND=dp), PUBLIC :: neutral_mass, mega_number, macro_mass  
  REAL(KIND=dp), PUBLIC :: timestep, my_tot_mass_loss
  REAL(KIND=dp), PUBLIC :: my_min_dt, my_max_dt, my_ave_dt
  INTEGER,       PUBLIC :: my_ntotal, my_n_chex
  INTEGER,       PUBLIC :: boundary_npartx, boundary_nparty, boundary_npartz  
  integer(8),    PUBLIC :: total_npart_verbosity
  
  
  type(kinetic_level), DIMENSION (:), ALLOCATABLE, TARGET, PUBLIC :: m_levels
  INTEGER, PUBLIC :: p_nlevels
  
  
  INTEGER,       PUBLIC :: pout_unit
  CHARACTER (LEN = lenstr), PUBLIC ::pout_name
  
  INTEGER,       PUBLIC :: my_nchex_events, my_nchex_current
  type(charge_exchange_event), DIMENSION (:), ALLOCATABLE, PUBLIC         :: my_chex_events
  type(particle_storage),      DIMENSION (:), ALLOCATABLE, TARGET, PUBLIC :: my_particles
  type(node_workload),         DIMENSION (:), ALLOCATABLE, PUBLIC         :: mpi_balancing
    
  INTEGER :: balance_interval
  
  INTEGER :: nprocMPI, my_rank, nthreads, my_thread
  
  INTEGER, DIMENSION(:),   ALLOCATABLE, PUBLIC :: mpi_requests    
  INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: mpi_statuses
  
  REAL(KIND=dp), PUBLIC :: last_time_step_sec
  
  ! accumulate neutral distribution statistics during calculation (requires extra memory)
  LOGICAL, PUBLIC :: aggregate_neutral_data
  

!$OMP threadprivate (my_nchex_events, my_nchex_current, my_chex_events)
!$OMP threadprivate (my_min_dt, my_max_dt, my_ave_dt)
!$OMP threadprivate (my_n_chex,my_tot_mass_loss)
!!$OMP threadprivate (pout_name,pout_unit)
!$OMP threadprivate (nprocMPI, my_rank, nthreads, my_thread)
  
  
  REAL(KIND=dp) :: int_r_maxwel(10,36,201), int_r_kappa(10,36,101), int_Maxwellian(0:400)
  REAL(KIND=dp) :: int_phi_maxwel(10,36,201,51), int_phi_kappa(10,36,101,51)        
  
  
!  REAL(KIND=dp), DIMENSION (:,:,:),   ALLOCATABLE, PUBLIC :: int_r_maxwel, int_r_kappa
!  REAL(KIND=dp), DIMENSION (:),       ALLOCATABLE, PUBLIC :: int_Maxwellian
!  REAL(KIND=dp), DIMENSION (:,:,:,:), ALLOCATABLE, PUBLIC :: int_phi_maxwel, int_phi_kappa
!!$OMP threadprivate (int_r_maxwel, int_r_kappa,int_phi_maxwel,int_phi_kappa, int_Maxwellian)
  

!
! parameters for kappa distribution in inner heliosheath
  REAL(KIND=dp), PUBLIC :: kappa, kappa_v_th_factor
  REAL(KIND=dp), PUBLIC :: kappa_v_th_factor_sq, gamma_vrel_factor
  REAL(KIND=dp), PUBLIC :: sigma_vrel_maxwel(10,36), sigma_vrel_kappa(10,36)
  
!  REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: sigma_vrel_maxwel, sigma_vrel_kappa  
!!$OMP threadprivate (sigma_vrel_maxwel, sigma_vrel_kappa)
  
  REAL(KIND=dp), PUBLIC :: g_xmin, g_xmax, g_ymin, g_ymax, g_zmin, g_zmax
  
  INTEGER, PARAMETER, PUBLIC :: nsplit_grids = 4
  REAL(KIND=dp), PUBLIC :: split_xmin(1:nsplit_grids) 
  REAL(KIND=dp), PUBLIC :: split_ymin(1:nsplit_grids) 
  REAL(KIND=dp), PUBLIC :: split_zmin(1:nsplit_grids) 
  REAL(KIND=dp), PUBLIC :: split_xmax(1:nsplit_grids) 
  REAL(KIND=dp), PUBLIC :: split_ymax(1:nsplit_grids) 
  REAL(KIND=dp), PUBLIC :: split_zmax(1:nsplit_grids) 
  
  REAL(KIND=dp), PUBLIC :: cell_ratios(0:nsplit_grids) 
  
  
END MODULE global
!
!***********************************************************************
!***********************************************************************
!
MODULE level1_subroutines
!
! Lowest level subroutines
!
  USE global
  USE random_gen_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: gamma_func, initialize_gammas, get_sigma_vrelp
  PUBLIC :: init_vr_phi_kappa, init_vr_phi_maxwel
  PUBLIC :: get_vr_phi_kappa , get_vr_phi_maxwel, get_Maxwellian_vr
  PUBLIC :: boundary_particles, plasma, source_chex, source_phi
  PUBLIC :: initial_particles_distribution
  PUBLIC :: collect_sources_node
  PUBLIC :: pout,init_pout

CONTAINS

  SUBROUTINE init_pout

    ! 1 log file per node
    
    pout_unit = 1000 + my_rank
    
    if (my_rank>0) then
      pout_unit = -1
      return
    endif
    
    write  ( pout_name, '(a6,i3.3,a2,i2.2)' ) 'fort.n', my_rank
    
    if (my_thread == 0) then      
      open(unit = pout_unit,file=pout_name,form='formatted',status='REPLACE' ,action='write')
      !close(pout_unit)          
    endif
    
    return
    
    ! 1 log file per thread
    
!    pout_unit = 1000 + my_rank*nprocMPI + my_thread
!    write  ( pout_name, '(a6,i3.3,a2,i2.2)' ) 'fort.n', my_rank, '.t', my_thread
!    open(unit = pout_unit,file=pout_name,form='formatted',status='REPLACE' ,action='write')
!    close(pout_unit)          
    
    
    
!    close(unit = pout_unit,status='delete')
!    open(unit = pout_unit,file=pout_name,status='unknown')
      !rewind(unit = pout_unit)
  end SUBROUTINE init_pout

  SUBROUTINE pout(str)
    character (len=*), INTENT(in) :: str
    
    !character (len=100), INTENT(in) :: str_final
    
    if (pout_unit < 0) then
      return
    endif
            
    !$omp critical (update_pout)
    !open(unit = pout_unit,file=pout_name,form='formatted',status='unknown', position = 'APPEND')
    write (pout_unit,'(a,i2,a2,a)') '(thread ', my_thread, '): ', str
    flush(pout_unit)
    !close(pout_unit)          
    !$omp end critical (update_pout)
  end SUBROUTINE pout

  SUBROUTINE init_vr_phi_maxwel
  
    INTEGER :: s,i,j,k  
    REAL(KIND=dp) :: rH, rp, int_p, int_phi_a, int_phi_b, phi_a, phi_b
    REAL(KIND=dp) :: ee, sigma_ex, vrel, norm_factor, int_vrel_phi
    REAL(KIND=dp) :: int_vrel_r, local_vT, local_vp, local_vH
    REAL(KIND=dp) :: drH, drp, dphi, dvr
    REAL(KIND=dp) :: vr, tot_int
    
    !allocate(int_Maxwellian(0:400))
    
    !allocate(int_r_maxwel(10,36,201),int_phi_maxwel(10,36,201,51))
    
    !allocate(sigma_vrel_maxwel(10,36))                 
        
        
    dvr = 0.01_dp       
    int_Maxwellian(0) = 0.0_dp
    DO i = 1, 400
      vr = (i-0.5_dp)*dvr
      int_Maxwellian(i) = int_Maxwellian(i-1) + vr*vr * EXP(-vr*vr) * dvr
    END DO

    tot_int = int_Maxwellian(400)
    int_Maxwellian = int_Maxwellian / tot_int
            


    sigma_vrel_maxwel = 0.0_dp

    drH = 1.0_dp/6.0_dp
    drp = 0.02_dp
    dphi = pi/50.0_dp

   !$omp parallel do default(none)  &
   !$omp& private(s,i,j,k, rH, rp, int_p, int_phi_a, int_phi_b, phi_a, phi_b, &
   !$omp& ee, sigma_ex, vrel, norm_factor, int_vrel_phi, &
   !$omp& int_vrel_r, local_vT, local_vp, local_vH, dvr,  vr, tot_int) &
   !$omp& shared (drH,drp,dphi,int_r_maxwel,sigma_vrel_maxwel,int_phi_maxwel)
    DO s = 1, 10
    ! 10 -- 302.5 km/s (6,000K -- 5,500,000K) use this to estimate cross-section
      local_vT = 2.5*(s+1)**2 * 1e3_dp

      DO i = 1, 36
         rH = ( i*drH )**2

         int_phi_a = 0.0_dp
         int_vrel_r = 0.0_dp
         int_r_maxwel(s,i,1) = 0.0_dp
         DO j = 2, 201
            rp = (j-1.0_dp)*drp

            int_vrel_phi = 0.0_dp
            int_p = 0.0_dp
            DO k = 2, 51
               phi_a = (k-2)*dphi
               phi_b = (k-1)*dphi

               local_vH = rH*local_vT
               local_vp = rp*local_vT
               vrel = SQRT(local_vH**2+local_vp**2-2.0_dp*local_vH*local_vp*COS(phi_b))
               ee = vrel*vrel * 5.21e-12_dp   ! energy in keV per nucleon
               sigma_ex = (4.15_dp-0.531_dp*LOG(ee))**2.0_dp * (1.0_dp-EXP(-67.3_dp/ee))**4.5_dp  ! cross-section in units of 10^-20 m
        

    !sigma_ex = 1.0_dp

    ! Trapezoid rule for integrating over phi -- uniform dphi assumed
               int_p = int_p + 0.5_dp*dphi*sigma_ex * rp*rp * EXP(-rp*rp) &
                    * (SQRT(rH*rH + rp*rp - 2.0_dp*rH*rp*COS(phi_a))  &
                     + SQRT(rH*rH + rp*rp - 2.0_dp*rH*rp*COS(phi_b)))

               int_vrel_phi = int_vrel_phi &
                    + sigma_ex * (local_vp)**2.0_dp * EXP(-rp*rp) &
                    * 0.5_dp*(SIN(phi_a)+SIN(phi_b)) * dphi * vrel

            END DO
            int_r_maxwel(s,i,j) = int_r_maxwel(s,i,j-1) + int_p*drp

            int_phi_b = int_vrel_phi
            int_vrel_r = int_vrel_r + 0.5_dp*(int_phi_a + int_phi_b) * local_vT * drp
            int_phi_a = int_phi_b

         END DO
         int_r_maxwel(s,i,:) = int_r_maxwel(s,i,:)/int_r_maxwel(s,i,201)

         norm_factor = 2.0_dp*pi/( SQRT(pi)*local_vT )**3.0_dp
         sigma_vrel_maxwel(s,i) = int_vrel_r * norm_factor * 1e-20_dp

      END DO


      DO i = 1, 36
         rH = ( i*drH )**2
      
         DO j = 1, 201
            rp = (j-0.5_dp)*drp

            int_phi_maxwel(s,i,j,1) = 0.0
            DO k = 2, 51
               phi_a = (k-2)*dphi
               phi_b = (k-1)*dphi

               local_vH = rH*local_vT
               local_vp = rp*local_vT
               vrel = SQRT(local_vH**2+local_vp**2-2.0_dp*local_vH*local_vp*COS(0.5_dp*(phi_a+phi_b)))
               ee = vrel*vrel * 5.21e-12_dp   ! energy in keV per nucleon
               sigma_ex = (4.15_dp-0.531_dp*LOG(ee))**2.0_dp * (1.0_dp-EXP(-67.3_dp/ee))**4.5_dp  ! cross-section in units of 10^-20 m

    !sigma_ex = 1.0_dp

               int_phi_maxwel(s,i,j,k) = int_phi_maxwel(s,i,j,k-1) + 0.5*dphi*sigma_ex * rp*rp &
           * (SQRT(rH*rH + rp*rp - 2.0_dp*rH*rp*COS(phi_a)) &
            + SQRT(rH*rH + rp*rp - 2.0_dp*rH*rp*COS(phi_b))) * EXP(-rp*rp)

            END DO

            int_phi_maxwel(s,i,j,:) = int_phi_maxwel(s,i,j,:)/int_phi_maxwel(s,i,j,51)
            
         END DO
      END DO

    END DO
   !$omp end parallel do

    
  END SUBROUTINE init_vr_phi_maxwel

  SUBROUTINE get_vr_phi_maxwel(vH,vT,vr,phi)
! Here we return the speed and direction (in the plasma velocity frame)
! of the new neutral after ch-ex. Follows Section 7.4 of Lipatov's book (2002)
    REAL(KIND=dp), INTENT(in)  :: vH, vT
    REAL(KIND=dp), INTENT(out) :: vr, phi    
    REAL(KIND=dp) :: j1, k1,dprand, a, fa, fb    
    INTEGER :: s, i, j, k, dj, dk, jlow, jhigh, klow, khigh    
    REAL(KIND=dp) :: drH, drp, dphi
    
    drH = 1.0_dp/6.0_dp
    drp = 0.02_dp
    dphi = pi/50.0_dp
    

    s = NINT( SQRT(vT/2.5e3_dp) - 1.0_dp )     ! vT in km/s
    s = MIN(s,10);  s = MAX(s,1)

    i = NINT( SQRT(vH/vT) /drH )
    i = MIN(i,36);  i = MAX(i,1)

    CALL rand_number(dprand)

! use bisection to find speed
    jlow = 1
    jhigh = 201
    dj = jhigh - jlow
    j = 50
    DO WHILE (dj .GT. 1)
       j1 = int_r_maxwel(s,i,j)
       IF (j1 .GT. dprand) THEN
          jhigh = j
       ELSE
          jlow = j
       END IF
       j = FLOOR((jhigh+jlow)/2.0_dp)
       dj = jhigh - jlow
    END DO
    j = jhigh
    vr = (j-1)*vT*drp

    CALL rand_number(dprand)
! WARNING: need to make sure we pick uniformly in 3D
    dprand = ACOS(0.999999 - 1.99999*dprand)/pi   ! ensure uniform 3D distribution

! use bisection to find angle
    klow = 1
    khigh = 51
    dk = khigh - klow
    k = 24
    DO WHILE (dk .GT. 1)
       k1 = int_phi_maxwel(s,i,j,k)
       IF (k1 .GT. dprand) THEN
          khigh = k
       ELSE
          klow = k
       END IF
       k = FLOOR((khigh+klow)/2.0_dp)
       dk = khigh - klow
    END DO
    k = khigh
! interpolate since we used Trapezoid rule
    fb = int_phi_maxwel(s,i,j,k);  fa = int_phi_maxwel(s,i,j,k-1)
    a = (k-2)*dphi
    phi = a + dphi/(fb-fa) * (dprand - fa)

  END SUBROUTINE get_vr_phi_maxwel
  
  
  SUBROUTINE init_vr_phi_kappa
    INTEGER :: s,i,j,k
    REAL(KIND=dp) :: rH, rp, rp1, rp2, int_p, j1, k1
    REAL(KIND=dp) :: vrel, ee, sigma_ex
    REAL(KIND=dp) :: phi_a, phi_b, local_vT, local_vH, local_vp
    REAL(KIND=dp) :: gamma_kp1, gamma_kmh, norm_factor
    REAL(KIND=dp) :: int_phi_a, int_phi_b, int_vrel_phi, int_vrel_r
    REAL(KIND=dp) :: drH, drp, dphi


    drH = 1.0_dp/6.0_dp
    drp = 0.1_dp
    dphi = pi/50.0_dp
    
    !allocate(int_r_kappa(10,36,101), int_phi_kappa(10,36,101,51))
    
    
!    allocate(sigma_vrel_kappa(10,36))

    sigma_vrel_kappa = 0.0_dp

    
    CALL gamma_func(kappa+1.0_dp,gamma_kp1)
    CALL gamma_func(kappa-0.5_dp,gamma_kmh)

   !$omp  parallel do default(private) &
   !$omp& shared(drH,drp,dphi,gamma_kp1,gamma_kmh,kappa, kappa_v_th_factor, &
   !$omp& int_r_kappa, sigma_vrel_kappa, int_phi_kappa)
    DO s = 1, 10
    ! 10 -- 302.5 km/s (6,000K -- 5,500,000K) use this to estimate cross-section
      local_vT = 2.5*(s+1)**2 * 1e3_dp

      DO i = 1, 36
         rH = ( i*drH )**2

         int_phi_a = 0.0_dp
         int_vrel_r = 0.0_dp
         int_r_kappa(s,i,1) = 0.0_dp
         DO j = 2, 101
            rp1 = ( (j-2)*drp )**2
            rp2 = ( (j-1)*drp )**2

            int_vrel_phi = 0.0_dp
            int_p = 0.0_dp
            DO k = 2, 51
               phi_a = (k-2)*dphi
               phi_b = (k-1)*dphi

               local_vH = rH*local_vT
               local_vp = rp2*local_vT
               vrel = SQRT(local_vH**2+local_vp**2-2.0_dp*local_vH*local_vp*COS(phi_b))
               ee = vrel*vrel * 5.21e-12_dp   ! energy in keV per nucleon
               sigma_ex = (4.15_dp-0.531_dp*LOG(ee))**2.0_dp * (1.0_dp-EXP(-67.3_dp/ee))**4.5_dp  ! cross-section in units of 10^-20 m

    ! Trapezoid rule for integrating over phi -- uniform dphi assumed
               int_p = int_p + sigma_ex * rp2*rp2 &
    * (SQRT(rH*rH + rp2*rp2 - 2.0_dp*rH*rp2*COS(phi_a)) &
    + SQRT(rH*rH + rp2*rp2 - 2.0_dp*rH*rp2*COS(phi_b))) &
    * (1.0_dp+rp2*rp2/kappa / kappa_v_th_factor**2.0_dp)**(-kappa-1.0_dp)

               int_vrel_phi = int_vrel_phi &
                    + (local_vp)**2.0_dp*(1.0_dp + 1.0_dp/kappa * (rp2/kappa_v_th_factor)**2.0_dp)**(-kappa-1.0_dp) &
                    * 0.5_dp*(SIN(phi_a)+SIN(phi_b)) * dphi * vrel * sigma_ex

            END DO
            int_r_kappa(s,i,j) = int_r_kappa(s,i,j-1) + int_p * (rp2-rp1)  ! nonuniform drp


            int_phi_b = int_vrel_phi
            int_vrel_r = int_vrel_r + 0.5_dp*(int_phi_a + int_phi_b) * local_vT*(rp2-rp1)  ! nonuniform drp
            int_phi_a = int_phi_b

         END DO
         int_r_kappa(s,i,:) = int_r_kappa(s,i,:)/int_r_kappa(s,i,101)

         norm_factor = 2.0_dp*pi/( SQRT(pi*kappa)*local_vT*kappa_v_th_factor )**3.0_dp * gamma_kp1/gamma_kmh
         sigma_vrel_kappa(s,i) = int_vrel_r * norm_factor * 1e-20_dp

      END DO


      DO i = 1, 36
         rH = ( i*drH )**2
         
         DO j = 1, 101
            rp = ( (j-1)*drp )**2
         
            int_phi_kappa(s,i,j,1) = 0.0
            DO k = 2, 51
               phi_a = (k-2)*dphi
               phi_b = (k-1)*dphi

               local_vH = rH*local_vT
               local_vp = rp*local_vT
               vrel = sqrt(local_vH**2+local_vp**2-2.0_dp*local_vH*local_vp*COS(phi_b))
               ee = vrel*vrel * 5.21e-12_dp   ! energy in keV per nucleon
               sigma_ex = (4.15_dp-0.531_dp*LOG(ee))**2.0_dp * (1.0_dp-EXP(-67.3_dp/ee))**4.5_dp  ! cross-section in units of 10^-20 m

    ! Trapezoid rule for integrating over phi -- uniform dphi assumed
               int_phi_kappa(s,i,j,k) = int_phi_kappa(s,i,j,k-1) + sigma_ex * rp*rp &
    * (SQRT(rH*rH + rp*rp - 2.0_dp*rH*rp*COS(phi_a)) &
    + SQRT(rH*rH + rp*rp - 2.0_dp*rH*rp*COS(phi_b))) &
    * (1.0_dp+rp*rp/kappa / kappa_v_th_factor**2.0_dp)**(-kappa-1.0_dp)
            END DO

            int_phi_kappa(s,i,j,:) = int_phi_kappa(s,i,j,:)/int_phi_kappa(s,i,j,51)

         END DO
      END DO
    END DO
   !$omp end parallel do

  END SUBROUTINE init_vr_phi_kappa
  
!
!***********************************************************************
!
  SUBROUTINE get_vr_phi_kappa(vH,vT,vr,phi)
! Here we return the speed and direction (in the plasma velocity frame)
! of the new neutral after ch-ex. Follows Section 7.4 of Lipatov's book (2002)
    REAL(KIND=dp), INTENT(in)  :: vH, vT
    REAL(KIND=dp), INTENT(out) :: vr, phi            
    REAL(KIND=dp) :: a, fa, fb, ra, rb, dprand, j1, k1
    INTEGER :: i, j, k, dj, dk, jlow, jhigh, klow, khigh, s
    REAL(KIND=dp) :: drH, drp, dphi
    
    drH = 1.0_dp/6.0_dp
    drp = 0.1_dp
    dphi = pi/50.0_dp

   
   
    s = NINT( SQRT(vT/2.5e3_dp) - 1.0_dp )     ! vT in km/s
    s = MIN(s,10);  s = MAX(s,1)

    i = NINT( SQRT(vH/vT) /drH )
    i = MIN(i,36);  i = MAX(i,1)

    CALL rand_number(dprand)

! use bisection to find speed
    jlow = 1
    jhigh = 101
    dj = jhigh - jlow
    j = 10
    DO WHILE (dj .GT. 1)
       j1 = int_r_kappa(s,i,j)
       IF (j1 .GT. dprand) THEN
          jhigh = j
       ELSE
          jlow = j
       END IF
       j = FLOOR((jhigh+jlow)/2.0_dp)
       dj = jhigh - jlow
    END DO
    j = jhigh
! interpolate since we used Trapezoid rule
    rb = ( (j-1)*drp )**2;  ra = ( (j-2)*drp )**2
    fb = int_r_kappa(s,i,j);  fa = int_r_kappa(s,i,j-1)
    vr = ( ra + (rb-ra)/(fb-fa) * (dprand - fa) ) * vT
    vr = MAX(vr,1e-4_dp*vT)       ! very small vr creates trouble...

    CALL rand_number(dprand)
! WARNING: need to make sure we pick uniformly in 3D
    dprand = ACOS(0.999999 - 1.99999*dprand)/pi   ! ensure uniform 3D distribution
! use bisection to find angle
    klow = 1
    khigh = 51
    dk = khigh - klow
    k = 24
    DO WHILE (dk .GT. 1)
       k1 = int_phi_kappa(s,i,j,k)
       IF (k1 .GT. dprand) THEN
          khigh = k
       ELSE
          klow = k
       END IF
       k = FLOOR((khigh+klow)/2.0_dp)
       dk = khigh - klow
    END DO
    k = khigh
! interpolate since we used Trapezoid rule
    fb = int_phi_kappa(s,i,j,k);  fa = int_phi_kappa(s,i,j,k-1)
    a = (k-2)*dphi
    phi = a + dphi/(fb-fa) * (dprand - fa)

  END SUBROUTINE get_vr_phi_kappa
!
!***********************************************************************
!
  SUBROUTINE get_Maxwellian_vr(vr)
! Here we randomly select a speed from a Maxwellian distribution.
! Follows Section 7.3 of Lipatov's book (2002)
    REAL(KIND=dp), INTENT(out) :: vr    
    REAL(KIND=dp) :: dvr, i1, dprand
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
  SUBROUTINE initial_particles_distribution
! inject a Maxwellian distributions into the whole domain (for debugging purposes only)
  
    INTEGER, PARAMETER :: sample_size = 20000
    REAL(KIND=dp) :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(KIND=dp) :: x0, x1, y0, y1, z0, z1, vx, vy, vz
    REAL(KIND=dp) :: phi, theta, vr, vxy
    REAL(KIND=dp), DIMENSION(5) :: dprand
    INTEGER :: i,ind,ub
    INTEGER(8) :: tot_part,a1,a2,a3
    
    REAL(KIND=dp), DIMENSION(sample_size) ::vx_v,vy_v,vz_v
    REAL(KIND=dp) :: vxm,vym,vzm,H_temp,vrel_sq,vrelx_sq,vrely_sq,vrelz_sq
    
    
    CHARACTER (LEN = 200) str_buf
  
    xmin = g_xmin;     xmax = g_xmax
    ymin = g_ymin;     ymax = g_ymax
    zmin = g_zmin;     zmax = g_zmax
    
    ub      = ubound(dprand, dim = 1)
        
    if (my_rank .EQ. 0) then
      
      
      open(unit = 668,file='vxvyvz.txt',form='formatted',status='REPLACE' ,action='write')      
      
      vxm = 0.0_dp
      vym = 0.0_dp
      vzm = 0.0_dp
      
      do i = 1, sample_size
        do ind = 1, ub
          CALL rand_number(dprand(ind))
        enddo
      
         CALL get_Maxwellian_vr(vr)
         vr = vr * LISM_v_thermal

         theta = ACOS(0.99999_dp - 1.9999_dp*dprand(4))
         vz = vr * COS(theta)
         vxy = SQRT(vr*vr - vz*vz)

         phi = 2.0_dp*pi*dprand(5)
         vy = vxy * COS(phi)
         vx = vxy * SIN(phi)

         vx_v(i) = (vx + LISM_vx)
         vy_v(i) = (vy + LISM_vy)
         vz_v(i) = (vz + LISM_vz)
         
              
         write(668,*) (1e-3_dp)*vx_v(i),(1e-3_dp)*vy_v(i),(1e-3_dp)*vz_v(i)
         
         vxm = vxm + vx_v(i)
         vym = vym + vy_v(i)
         vzm = vzm + vz_v(i)
                  
      enddo
      close(668)
        
    
      vxm = vxm/sample_size
      vym = vym/sample_size
      vzm = vzm/sample_size
      
      vrel_sq = 0.0_dp
      vrelx_sq = 0.0_dp
      vrely_sq = 0.0_dp
      vrelz_sq = 0.0_dp
      do i = 1, sample_size
        vrel_sq = vrel_sq + (vx_v(i) - vxm )**2 &
                + (vy_v(i) - vym )**2 &
                + (vz_v(i) - vzm )**2
                
        vrelx_sq = vrelx_sq + (vx_v(i) - vxm )**2
        vrely_sq = vrely_sq + (vy_v(i) - vym )**2
        vrelz_sq = vrelz_sq + (vz_v(i) - vzm )**2
      enddo
      vrel_sq  = vrel_sq/sample_size
      vrelx_sq = vrelx_sq/sample_size
      vrely_sq = vrely_sq/sample_size
      vrelz_sq = vrelz_sq/sample_size
      
      H_temp = vrel_sq * neutral_mass / (3.0*k_B)
      
      print *, vrelx_sq*neutral_mass/k_B, vrely_sq*neutral_mass/k_B, vrelz_sq*neutral_mass/k_B
      print *,'H_temp ',H_temp
    endif
   
    return
    

    DO i = 1, boundary_npartx
    
       do ind = 1, ub
        CALL rand_number(dprand(ind))
       enddo
               
       
! fill x faces from yzmin-boundary_thickness to yzmax+boundary_thickness
       x0 = xmin + dprand(1)*(xmax-xmin)
       y0 = ymin + dprand(2)*(ymax-ymin)
       z0 = zmin + dprand(3)*(zmax-zmin)
                        
       ind = i

       my_particles(my_thread)%particles__(1,ind) = x0
       my_particles(my_thread)%particles__(2,ind) = y0
       my_particles(my_thread)%particles__(3,ind) = z0
       CALL get_Maxwellian_vr(vr)
       vr = vr * LISM_v_thermal
! must AVOID setting vx=vy=0
       theta = ACOS(0.99999_dp - 1.9999_dp*dprand(4))
       vz = vr * COS(theta)
       vxy = SQRT(vr*vr - vz*vz)

       phi = 2.0_dp*pi*dprand(5)
       vy = vxy * COS(phi)
       vx = vxy * SIN(phi)

       my_particles(my_thread)%particles__(4,ind) = vx + LISM_vx
       my_particles(my_thread)%particles__(5,ind) = vy + LISM_vy
       my_particles(my_thread)%particles__(6,ind) = vz + LISM_vz

       my_particles(my_thread)%particles__(7,ind) = 0.1  ! region 0 (LISM)
       my_particles(my_thread)%particles__(0,ind) = 1.0
       my_particles(my_thread)%particles__(8,ind) = 1.0

    enddo
    
    my_particles(my_thread)%my_npart = my_ntotal
    
    a1 = my_ntotal
    a2 = nprocMPI
    a3 = nthreads
    tot_part = a1*a2*a3
    if (tot_part < 0) then
      write (str_buf,*) 'OMG!!! tot_part < 0 ', a1,a2,a3
      call pout(str_buf)   
    endif
    
    !write (str_buf,*) 'number of particles per thread: ', my_ntotal, ' total particles: ', tot_part
    !call pout(str_buf)   
    !write (str_buf,*) ' total particles: ', tot_part
    !call pout(str_buf)   
    
    
    
    
  end SUBROUTINE initial_particles_distribution
  
!
!***********************************************************************
!
  SUBROUTINE boundary_particles
! inject a Maxwellian distributions into the boundary layer
    REAL(KIND=dp) :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(KIND=dp) :: x0, x1, y0, y1, z0, z1, vx, vy, vz
    REAL(KIND=dp) :: phi, theta, vr, vxy
    REAL(KIND=dp), DIMENSION(10) :: dprand
    INTEGER :: i,ind,ub
    
    REAL(KIND=dp), DIMENSION (:,:), POINTER :: particles
        
  
    xmin = g_xmin;     xmax = g_xmax
    ymin = g_ymin;     ymax = g_ymax
    zmin = g_zmin;     zmax = g_zmax
    
    ub = ubound(dprand, dim = 1)
    
    particles => my_particles(my_thread)%particles__

    DO i = 1, boundary_npartx
    
       do ind = 1, ub
        CALL rand_number(dprand(ind))
       enddo
               
       
! fill x faces from yzmin-boundary_thickness to yzmax+boundary_thickness
       x0 = xmax + dprand(1)*boundary_thickness
       y0 = (ymin-boundary_thickness)*(1.0-dprand(2)) &
          + (ymax+boundary_thickness)*dprand(2)
       z0 = (zmin-boundary_thickness)*(1.0-dprand(3)) &
          + (zmax+boundary_thickness)*dprand(3)
          
       ind = my_particles(my_thread)%my_npart+i

       particles(1,ind) = x0
       particles(2,ind) = y0
       particles(3,ind) = z0
       
       CALL get_Maxwellian_vr(vr)
       vr = vr * LISM_v_thermal
! must AVOID setting vx=vy=0
       theta = ACOS(0.99999_dp - 1.9999_dp*dprand(4))
       vz = vr * COS(theta)
       vxy = SQRT(vr*vr - vz*vz)

       phi = 2.0_dp*pi*dprand(5)
       vy = vxy * COS(phi)
       vx = vxy * SIN(phi)

       particles(4,ind) = vx + LISM_vx
       particles(5,ind) = vy + LISM_vy
       particles(6,ind) = vz + LISM_vz

       particles(7,ind) = 0.1  ! region 0 (LISM)
       particles(0,ind) = 1.0
       particles(8,ind) = 1.0

       x1 = xmin - dprand(6)*boundary_thickness
       y1 = (ymin-boundary_thickness)*(1.0-dprand(7)) &
          + (ymax+boundary_thickness)*dprand(7)
       z1 = (zmin-boundary_thickness)*(1.0-dprand(8)) &
          + (zmax+boundary_thickness)*dprand(8)
          
          
       ind = ind+boundary_npartx

       particles(1,ind) = x1
       particles(2,ind) = y1
       particles(3,ind) = z1

       CALL get_Maxwellian_vr(vr)
       vr = vr * LISM_v_thermal
! must AVOID setting vx=vy=0
       theta = ACOS(0.99999_dp - 1.9999_dp*dprand(9))
       vz = vr * COS(theta)
       vxy = SQRT(vr*vr - vz*vz)

       phi = 2.0_dp*pi*dprand(10)
       vy = vxy * COS(phi)
       vx = vxy * SIN(phi)

       particles(4,ind) = vx + LISM_vx
       particles(5,ind) = vy + LISM_vy
       particles(6,ind) = vz + LISM_vz

       particles(7,ind) = 0.1  ! region 0 (LISM)
       particles(0,ind) = 1.0
       particles(8,ind) = 1.0

    END DO
    my_particles(my_thread)%my_npart = my_particles(my_thread)%my_npart + 2*boundary_npartx

    DO i = 1, boundary_nparty

       do ind = 1, ub
        CALL rand_number(dprand(ind))
       enddo
       
       ind = my_particles(my_thread)%my_npart+i
       
! fill y faces from xmin to xmax & zmin to zmax
       y0 = ymax + dprand(1)*boundary_thickness
       x0 = xmin*(1.0-dprand(2)) + xmax*dprand(2)
       z0 = zmin*(1.0-dprand(3)) + zmax*dprand(3)

       particles(1,ind) = x0
       particles(2,ind) = y0
       particles(3,ind) = z0
       
       
       CALL get_Maxwellian_vr(vr)
       vr = vr * LISM_v_thermal
! must AVOID setting vx=vy=0
       theta = ACOS(0.99999_dp - 1.9999_dp*dprand(4))
       vz = vr * COS(theta)
       vxy = SQRT(vr*vr - vz*vz)

       phi = 2.0_dp*pi*dprand(5)
       vy = vxy * COS(phi)
       vx = vxy * SIN(phi)

       particles(4,ind) = vx + LISM_vx
       particles(5,ind) = vy + LISM_vy
       particles(6,ind) = vz + LISM_vz

       particles(7,ind) = 0.1  ! region 0 (LISM)
       particles(0,ind) = 1.0
       particles(8,ind) = 1.0

       y1 = ymin - dprand(6)*boundary_thickness
       x1 = xmin*(1.0-dprand(7)) + xmax*dprand(7)
       z1 = zmin*(1.0-dprand(8)) + zmax*dprand(8)
       
       ind = ind+boundary_nparty

       particles(1,ind) = x1
       particles(2,ind) = y1
       particles(3,ind) = z1

       CALL get_Maxwellian_vr(vr)
       vr = vr * LISM_v_thermal
! must AVOID setting vx=vy=0
       theta = ACOS(0.99999_dp - 1.9999_dp*dprand(9))
       vz = vr * COS(theta)
       vxy = SQRT(vr*vr - vz*vz)

       phi = 2.0_dp*pi*dprand(10)
       vy = vxy * COS(phi)
       vx = vxy * SIN(phi)

       particles(4,ind) = vx + LISM_vx
       particles(5,ind) = vy + LISM_vy
       particles(6,ind) = vz + LISM_vz

       particles(7,ind) = 0.1  ! region 0 (LISM)
       particles(0,ind) = 1.0
       particles(8,ind) = 1.0

    END DO
    my_particles(my_thread)%my_npart = my_particles(my_thread)%my_npart + 2*boundary_nparty

    DO i = 1, boundary_npartz

       do ind = 1, ub
        CALL rand_number(dprand(ind))
       enddo
       
       ind = my_particles(my_thread)%my_npart+i
       
! fill z faces from xmin to xmax & ymin to ymax
       z0 = zmax + dprand(1)*boundary_thickness
       y0 = ymin*(1.0-dprand(2)) + ymax*dprand(2)
       x0 = xmin*(1.0-dprand(3)) + xmax*dprand(3)

       particles(1,ind) = x0
       particles(2,ind) = y0
       particles(3,ind) = z0
       
       CALL get_Maxwellian_vr(vr)
       vr = vr * LISM_v_thermal
! must AVOID setting vx=vy=0
       theta = ACOS(0.99999_dp - 1.9999_dp*dprand(4))
       vz = vr * COS(theta)
       vxy = SQRT(vr*vr - vz*vz)

       phi = 2.0_dp*pi*dprand(5)
       vy = vxy * COS(phi)
       vx = vxy * SIN(phi)

       particles(4,ind) = vx + LISM_vx
       particles(5,ind) = vy + LISM_vy
       particles(6,ind) = vz + LISM_vz

       particles(7,ind) = 0.1  ! region 0 (LISM)
       particles(0,ind) = 1.0
       particles(8,ind) = 1.0

       z1 = zmin - dprand(6)*boundary_thickness
       y1 = ymin*(1.0-dprand(7)) + ymax*dprand(7)
       x1 = xmin*(1.0-dprand(8)) + xmax*dprand(8)
       
       ind = ind+boundary_npartz

       particles(1,ind) = x1
       particles(2,ind) = y1
       particles(3,ind) = z1

       CALL get_Maxwellian_vr(vr)
       vr = vr * LISM_v_thermal
! must AVOID setting vx=vy=0
       theta = ACOS(0.99999_dp - 1.9999_dp*dprand(9))
       vz = vr * COS(theta)
       vxy = SQRT(vr*vr - vz*vz)

       phi = 2.0_dp*pi*dprand(10)
       vy = vxy * COS(phi)
       vx = vxy * SIN(phi)

       particles(4,ind) = vx + LISM_vx
       particles(5,ind) = vy + LISM_vy
       particles(6,ind) = vz + LISM_vz

       particles(7,ind) = 0.1  ! region 0 (LISM)
       particles(0,ind) = 1.0
       particles(8,ind) = 1.0

    END DO
    my_particles(my_thread)%my_npart = my_particles(my_thread)%my_npart + 2*boundary_npartz

  END SUBROUTINE boundary_particles
!
!***********************************************************************
!
  SUBROUTINE GAMMA_func(X,GA)
!!$
!!$       ==================================================
!!$       Purpose: Compute the gamma function â(x)
!!$       Input :  x  --- Argument of â(x)
!!$                       ( x is not equal to 0,-1,-2,úúú )
!!$       Output:  GA --- â(x)
!!$       ==================================================
!!$
    REAL(KIND=dp) :: X, GA, G(26), Z, R, GR
    INTEGER :: K, M1, M

    IF (X.EQ.INT(X)) THEN
       IF (X.GT.0.0_dp) THEN
          GA=1.0_dp
          M1=NINT(X)-1
          DO K=2,M1
             GA=GA*K
          END DO
       ELSE
          GA=1e300_dp
       ENDIF
    ELSE
       IF (DABS(X).GT.1.0_dp) THEN
          Z=DABS(X)
          M=NINT(Z)
          R=1.0_dp
          DO K=1,M
             R=R*(Z-K)
          END DO
          Z=Z-M
       ELSE
          Z=X
       ENDIF
       DATA G/1.0D0,0.5772156649015329D0, &
            -0.6558780715202538D0, -0.420026350340952D-1, &
            0.1665386113822915D0,-.421977345555443D-1, &
            -.96219715278770D-2, .72189432466630D-2, &
            -.11651675918591D-2, -.2152416741149D-3, &
            .1280502823882D-3, -.201348547807D-4, &
            -.12504934821D-5, .11330272320D-5, &
            -.2056338417D-6, .61160950D-8, &
            .50020075D-8, -.11812746D-8, &
            .1043427D-9, .77823D-11, &
            -.36968D-11, .51D-12, &
            -.206D-13, -.54D-14, .14D-14, .1D-15/
       GR=G(26)
       DO K=25,1,-1
          GR=GR*Z+G(K)
       END DO

       GA=1.0_dp/(GR*Z)
       IF (DABS(X).GT.1.0_dp) THEN
          GA=GA*R
          IF (X.LT.0.0_dp) GA=-pi/(X*GA*DSIN(pi*X))
       ENDIF
    ENDIF
    RETURN

  END SUBROUTINE GAMMA_FUNC
!
!***********************************************************************
!
  SUBROUTINE initialize_gammas
    REAL(KIND=dp) :: gamma_kp1, gamma_kmh

    CALL gamma_func(kappa+1.0_dp,gamma_kp1)
    CALL gamma_func(kappa-0.5_dp,gamma_kmh)

    gamma_vrel_factor = gamma_kp1*gamma_kp1 / &
         ( kappa*(kappa-1.0_dp)*(kappa-1.0_dp) * gamma_kmh*gamma_kmh )

  END SUBROUTINE initialize_gammas
!
!***********************************************************************
!
  SUBROUTINE get_sigma_vrelp(vx,vy,vz,ux,uy,uz,v_p_th_sq,sigma_vrelp, use_kappa)
    REAL(KIND=dp), INTENT (in) :: vx, vy, vz, ux, uy, uz, v_p_th_sq
    REAL(KIND=dp), INTENT (out) :: sigma_vrelp
    REAL(KIND=dp) :: delta_u_sq, s, sigma_vrel1, sigma_vrel2
    INTEGER, INTENT (in) :: use_kappa
    INTEGER :: s1, s2, i

    delta_u_sq = (vx-ux)**2 + (vy-uy)**2 + (vz-uz)**2

! use computed integrals that include variable cross-section
! linearly interpolate thermal speed
    s = SQRT(SQRT(v_p_th_sq)/2.5e3_dp) - 1.0_dp      ! vT in km/s
    s = MIN(s,9.999999_dp);  s = MAX(s,1.0000001_dp)
    s1 = FLOOR(s);  s2 = s1 + 1

    i = NINT( SQRT( SQRT(delta_u_sq/v_p_th_sq) ) *6.0_dp )
    i = MIN(i,36);  i = MAX(i,1)

    IF (use_kappa .EQ. 1) THEN

       sigma_vrel1 = sigma_vrel_kappa(s1,i)
       sigma_vrel2 = sigma_vrel_kappa(s2,i)

    ELSE

       sigma_vrel1 = sigma_vrel_maxwel(s1,i)
       sigma_vrel2 = sigma_vrel_maxwel(s2,i)

    END IF

    sigma_vrelp = sigma_vrel1 + (sigma_vrel2-sigma_vrel1) * (s - s1) !/(s2-s1)

  END SUBROUTINE get_sigma_vrelp
!
!***********************************************************************
!
  SUBROUTINE plasma(xi,yi,zi,n_p,T_p,ux,uy,uz,region_p)
    REAL(KIND=dp), INTENT(in)  :: xi, yi, zi
    REAL(KIND=dp), INTENT(out) :: n_p, ux, uy, uz, T_p, region_p
!    REAL(KIND=dp) :: ri, risq, mod_u
    INTEGER :: i,j,k,ilev
       
    
! return plasma number density, velocity components and temp at a given point
! here we DO NOT interpolate, which improves speed

 !   risq = xi*xi + yi*yi + zi*zi
 !   ri = SQRT(risq)
    
    do ilev = p_nlevels - 1, 0, -1
      if ((m_levels(ilev)%p_xmin<=xi).and.(xi<=m_levels(ilev)%p_xmax).and. &
          (m_levels(ilev)%p_ymin<=yi).and.(yi<=m_levels(ilev)%p_ymax).and. &
          (m_levels(ilev)%p_zmin<=zi).and.(zi<=m_levels(ilev)%p_zmax)) then
          
        i = floor((xi - m_levels(ilev)%p_xmin)/m_levels(ilev)%p_dx)
        j = floor((yi - m_levels(ilev)%p_ymin)/m_levels(ilev)%p_dy)
        k = floor((zi - m_levels(ilev)%p_zmin)/m_levels(ilev)%p_dz)
        
        i = min(max(0,i),m_levels(ilev)%p_nx-1)
        j = min(max(0,j),m_levels(ilev)%p_ny-1)
        k = min(max(0,k),m_levels(ilev)%p_nz-1)
        
        n_p = m_levels(ilev)%g_plasma(i,j,k,WRHO)
        T_p = m_levels(ilev)%g_plasma(i,j,k,WTEMP)
        ux  = m_levels(ilev)%g_plasma(i,j,k,WVELX)
        uy  = m_levels(ilev)%g_plasma(i,j,k,WVELY)
        uz  = m_levels(ilev)%g_plasma(i,j,k,WVELZ)
        region_p = m_levels(ilev)%g_plasma(i,j,k,WREG)                         
        
        return
        exit ! should not be here
      endif          
    enddo

    
! LISM
    region_p = 0.1
    n_p = LISM_nH
    ux  = LISM_vx
    uy  = LISM_vy
    uz  = LISM_vz
    T_p = LISM_TH
    
    return

  END SUBROUTINE plasma
!
!***********************************************************************
!
  SUBROUTINE source_chex(mass,xi,yi,zi,vxold,vyold,vzold,vxnew,vynew,vznew)
    REAL(KIND=dp), INTENT(in)  :: mass, xi, yi, zi, vxold, vzold, vyold
    REAL(KIND=dp), INTENT(in)  :: vxnew, vznew, vynew
    REAL(KIND=dp) :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(KIND=dp) :: dx, dy, dz, dvx, dvy, dvz, dvsq
    INTEGER :: ilev, i, j, k, level
    
          
    level = -1
    
    do ilev = 0, p_nlevels - 1 
      if ((m_levels(ilev)%p_xmin<=xi).and.(xi<=m_levels(ilev)%p_xmax).and. &
          (m_levels(ilev)%p_ymin<=yi).and.(yi<=m_levels(ilev)%p_ymax).and. &
          (m_levels(ilev)%p_zmin<=zi).and.(zi<=m_levels(ilev)%p_zmax)) then
        
        level = ilev  
        cycle
      else
        level = ilev-1                       
        exit 
      endif          
    enddo
    
    if (level < 0) then
      return
    endif
        
    
    dvx = vxold - vxnew
    dvy = vyold - vynew
    dvz = vzold - vznew
    dvsq = (vxold*vxold + vzold*vzold + vyold*vyold &
          - vxnew*vxnew - vznew*vznew - vynew*vynew)

! add source to all levels up to the one we're on
    DO ilev = 0, level

       i = floor((xi - m_levels(ilev)%p_xmin)/m_levels(ilev)%p_dx)
       j = floor((yi - m_levels(ilev)%p_ymin)/m_levels(ilev)%p_dy)
       k = floor((zi - m_levels(ilev)%p_zmin)/m_levels(ilev)%p_dz)
        
       i = min(max(0,i),m_levels(ilev)%p_nx-1)
       j = min(max(0,j),m_levels(ilev)%p_ny-1)
       k = min(max(0,k),m_levels(ilev)%p_nz-1)
       
       
       my_nchex_current = my_nchex_current + 1
       my_chex_events(my_nchex_current)%level = ilev
       my_chex_events(my_nchex_current)%i = i
       my_chex_events(my_nchex_current)%j = j
       my_chex_events(my_nchex_current)%k = k
       my_chex_events(my_nchex_current)%source_m    = 0.0_dp
       my_chex_events(my_nchex_current)%source_px   = dvx  * mass
       my_chex_events(my_nchex_current)%source_py   = dvy  * mass
       my_chex_events(my_nchex_current)%source_pz   = dvz  * mass
       my_chex_events(my_nchex_current)%source_mvsq = dvsq * mass
              
       if (my_nchex_current>=my_nchex_events) then
        call collect_sources_node        
       endif

    END DO

  END SUBROUTINE source_chex
!
!***********************************************************************
!
  SUBROUTINE source_phi(xi,yi,zi,vx,vy,vz,mass_loss)
    REAL(KIND=dp), INTENT(in)  :: xi, yi, zi, vx, vz, vy, mass_loss
    REAL(KIND=dp) :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(KIND=dp) :: dx, dy, dz, vsq
    INTEGER :: level, i, j, k, ilev
   
          
    level = -1
    
    do ilev = 0, p_nlevels - 1 
      if ((m_levels(ilev)%p_xmin<=xi).and.(xi<=m_levels(ilev)%p_xmax).and. &
          (m_levels(ilev)%p_ymin<=yi).and.(yi<=m_levels(ilev)%p_ymax).and. &
          (m_levels(ilev)%p_zmin<=zi).and.(zi<=m_levels(ilev)%p_zmax)) then
        
        level = ilev
        cycle
      else
        level = ilev-1                       
        exit 
      endif          
    enddo
    
    if (level < 0) then
      return
    endif
    
    vsq = (vx*vx + vy*vy + vz*vz)


! add source to all levels up to the one we're on

    DO ilev = 0, level

       i = floor((xi - m_levels(ilev)%p_xmin)/m_levels(ilev)%p_dx)
       j = floor((yi - m_levels(ilev)%p_ymin)/m_levels(ilev)%p_dy)
       k = floor((zi - m_levels(ilev)%p_zmin)/m_levels(ilev)%p_dz)
        
       i = min(max(0,i),m_levels(ilev)%p_nx-1)
       j = min(max(0,j),m_levels(ilev)%p_ny-1)
       k = min(max(0,k),m_levels(ilev)%p_nz-1)
       
       
       my_nchex_current = my_nchex_current + 1
       my_chex_events(my_nchex_current)%level = ilev
       my_chex_events(my_nchex_current)%i = i
       my_chex_events(my_nchex_current)%j = j
       my_chex_events(my_nchex_current)%k = k
       my_chex_events(my_nchex_current)%source_m    = mass_loss
       my_chex_events(my_nchex_current)%source_px   = vx  * mass_loss
       my_chex_events(my_nchex_current)%source_py   = vy  * mass_loss
       my_chex_events(my_nchex_current)%source_pz   = vz  * mass_loss
       my_chex_events(my_nchex_current)%source_mvsq = vsq * mass_loss
              
       if (my_nchex_current>=my_nchex_events) then
        call collect_sources_node        
       endif

    END DO
       
  END SUBROUTINE source_phi
  
  SUBROUTINE collect_sources_node
  
  INTEGER :: nn,i,j,k,ilev
  CHARACTER (LEN = 100) str_buf
    
  !write (str_buf,*) 'collect_sources_node, nch_ex : ',my_nchex_current
  !call pout(str_buf)      
  
  
  !$omp critical (update_global_sources)
   do nn = 1, my_nchex_current
    ilev =  my_chex_events(nn)%level
    i =  my_chex_events(nn)%i
    j =  my_chex_events(nn)%j
    k =  my_chex_events(nn)%k
    
    m_levels(ilev)%source_m(i,j,k)    = m_levels(ilev)%source_m(i,j,k) + my_chex_events(nn)%source_m
    m_levels(ilev)%source_px(i,j,k)   = m_levels(ilev)%source_px(i,j,k) + my_chex_events(nn)%source_px
    m_levels(ilev)%source_py(i,j,k)   = m_levels(ilev)%source_py(i,j,k) + my_chex_events(nn)%source_py
    m_levels(ilev)%source_pz(i,j,k)   = m_levels(ilev)%source_pz(i,j,k) + my_chex_events(nn)%source_pz
    m_levels(ilev)%source_mvsq(i,j,k) = m_levels(ilev)%source_mvsq(i,j,k) + my_chex_events(nn)%source_mvsq
    
    if (my_chex_events(nn)%source_m == 0.0_dp) then
    ! exclude photoionization
      m_levels(ilev)%source_nchex(i,j,k) = m_levels(ilev)%source_nchex(i,j,k) + 1
    endif
   
   enddo       
   !$omp end critical (update_global_sources)
   
   my_nchex_current = 0               
  
  END SUBROUTINE collect_sources_node
  
  

END MODULE level1_subroutines
!
!***********************************************************************
!***********************************************************************
!
MODULE level2_subroutines
!
! Second level subroutines
!
  USE global
  USE level1_subroutines
  USE random_gen_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: mc_init_call_1, mc_init_call_2, setup_level, move, split, exchange, balance_npart
  PUBLIC :: output_raw, output_raw_threads, load_raw, output_source, output_grid,collect_neutral_distr
  PUBLIC :: mc_neutrals3d
  PUBLIC :: finalize

CONTAINS

  SUBROUTINE setup_level(au_level,au_ref_ratio,au_Lo,au_Hi,au_dx,au_size, &
                         au_source_m, au_source_px, au_source_py, au_source_pz, au_source_mvsq, &
                         au_plasma, &
                         au_nchex, au_H_dens, au_H_ux, au_H_uy, au_H_uz, au_H_temp)
  
  
    INTEGER, INTENT(in)                         :: au_level, au_ref_ratio
    REAL(KIND=dp), DIMENSION(0:2), INTENT(in)   :: au_Lo
    REAL(KIND=dp), DIMENSION(0:2), INTENT(in)   :: au_Hi
    REAL(KIND=dp), INTENT(in)           :: au_dx
    INTEGER, DIMENSION(0:2), INTENT(in) :: au_size
    
    REAL(KIND=dp), DIMENSION (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1), TARGET, INTENT(in) :: au_source_m 
    REAL(KIND=dp), DIMENSION (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1), TARGET, INTENT(in) :: au_source_px 
    REAL(KIND=dp), DIMENSION (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1), TARGET, INTENT(in) :: au_source_py 
    REAL(KIND=dp), DIMENSION (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1), TARGET, INTENT(in) :: au_source_pz 
    REAL(KIND=dp), DIMENSION (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1), TARGET, INTENT(in) :: au_source_mvsq    
    
    REAL(KIND=dp), DIMENSION (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1,0:WCOMP-1), TARGET, INTENT(in) :: au_plasma    
            
    REAL(KIND=dp), DIMENSION (-1:au_size(0), -1:au_size(1), -1:au_size(2)),   TARGET, INTENT(in) :: au_nchex, au_H_dens, & 
                        au_H_ux, au_H_uy, au_H_uz, au_H_temp 
    
    m_levels(au_level)%m_ref_ratio = au_ref_ratio
    
    m_levels(au_level)%p_xmin = au_Lo(0)*au
    m_levels(au_level)%p_xmax = au_Hi(0)*au
    m_levels(au_level)%p_ymin = au_Lo(1)*au
    m_levels(au_level)%p_ymax = au_Hi(1)*au
    m_levels(au_level)%p_zmin = au_Lo(2)*au 
    m_levels(au_level)%p_zmax = au_Hi(2)*au
    
    if (au_level == 0) then
      g_xmax = m_levels(0)%p_xmax;       g_xmin = m_levels(0)%p_xmin
      g_ymax = m_levels(0)%p_ymax;       g_ymin = m_levels(0)%p_ymin
      g_zmax = m_levels(0)%p_zmax;       g_zmin = m_levels(0)%p_zmin
    endif 
    
    m_levels(au_level)%sizex = m_levels(au_level)%p_xmax - m_levels(au_level)%p_xmin
    m_levels(au_level)%sizey = m_levels(au_level)%p_ymax - m_levels(au_level)%p_ymin
    m_levels(au_level)%sizez = m_levels(au_level)%p_zmax - m_levels(au_level)%p_zmin
    
    m_levels(au_level)%p_dx = au_dx*au
    m_levels(au_level)%p_dy = au_dx*au
    m_levels(au_level)%p_dz = au_dx*au
    
    m_levels(au_level)%p_nx = au_size(0)
    m_levels(au_level)%p_ny = au_size(1) 
    m_levels(au_level)%p_nz = au_size(2)
    
    m_levels(au_level)%g_source_m    => au_source_m 
    m_levels(au_level)%g_source_px   => au_source_px 
    m_levels(au_level)%g_source_py   => au_source_py 
    m_levels(au_level)%g_source_pz   => au_source_pz 
    m_levels(au_level)%g_source_mvsq => au_source_mvsq
    
    m_levels(au_level)%g_plasma      => au_plasma
    
    allocate(m_levels(au_level)%source_m   (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1)) 
    allocate(m_levels(au_level)%source_px  (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1)) 
    allocate(m_levels(au_level)%source_py  (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1)) 
    allocate(m_levels(au_level)%source_pz  (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1)) 
    allocate(m_levels(au_level)%source_mvsq(0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1)) 
    
    allocate(m_levels(au_level)%source_nchex(-1:au_size(0), -1:au_size(1), -1:au_size(2))) 
    
    
    
    m_levels(au_level)%g_source_nchex    => au_nchex
    m_levels(au_level)%g_H_dens => au_H_dens
    m_levels(au_level)%g_H_ux   => au_H_ux
    m_levels(au_level)%g_H_uy   => au_H_uy
    m_levels(au_level)%g_H_uz   => au_H_uz
    m_levels(au_level)%g_H_temp => au_H_temp    
    
    m_levels(au_level)%num_samples = 0
    
    if (aggregate_neutral_data) then
      allocate(m_levels(au_level)%H_dens(-1:au_size(0), -1:au_size(1), -1:au_size(2)))
      allocate(m_levels(au_level)%H_ux  (-1:au_size(0), -1:au_size(1), -1:au_size(2)))
      allocate(m_levels(au_level)%H_uy  (-1:au_size(0), -1:au_size(1), -1:au_size(2)))
      allocate(m_levels(au_level)%H_uz  (-1:au_size(0), -1:au_size(1), -1:au_size(2)))
      allocate(m_levels(au_level)%H_temp(-1:au_size(0), -1:au_size(1), -1:au_size(2)))
    endif
    
        
        
    
  END SUBROUTINE setup_level
    
    
  ! Initialization mechanism sequance:
  !    1. mc_init_level
  !    2. setup_level
  !
  
  SUBROUTINE mc_init_call_1(au_nlevels)
    INTEGER, INTENT(in) :: au_nlevels
    
    p_nlevels = au_nlevels
    allocate(m_levels(0:au_nlevels-1))
    
    aggregate_neutral_data = .false.
        
  END SUBROUTINE mc_init_call_1
    
    
    
  SUBROUTINE mc_init_call_2(au_LISM_nH, au_LISM_vH, au_LISM_TH,  & 
                  au_total_neutrals, au_restart, au_load_file)
    use hdf5
    
    REAL(KIND=dp), INTENT(in) :: au_LISM_nH, au_LISM_TH
    REAL(KIND=dp), DIMENSION(0:2), INTENT(in) :: au_LISM_vH        
    REAL(KIND=dp), INTENT(in) :: au_total_neutrals
    INTEGER, INTENT(in)       :: au_restart
    CHARACTER (LEN = lenstr), INTENT(in)  :: au_load_file
 
 !!!!!!!!!     
    REAL(KIND=dp) :: achieved_density, a, b, total_neutrals
    REAL(KIND=dp) :: tmp_vx, tmp_vy, tmp_vz, volume, tot_volume, tot_volume_spliting
    REAL(KIND=dp) :: xmin, xmax, ymin, ymax, zmin, zmax        
    CHARACTER (LEN = lenstr) :: dfile
    INTEGER(8) :: npart
    
    
 !!!!!!!!
 
    INTEGER, DIMENSION(1:4) :: put_seq
    INTEGER :: n, ierr, ncores, i, kk, ilev
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
    allocate(mpi_requests(0:loc_nprocMPI-1))    
    allocate(mpi_statuses(MPI_STATUS_SIZE,0:loc_nprocMPI-1))
    
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
    
    balance_interval = 10
    
    call init_pout()  
    
    
    !$omp parallel private(put_seq,n,aux_num,i,str_buf)
#ifdef CH_OMPCPP        
      my_thread = OMP_GET_THREAD_NUM()          
      nthreads  = omp_get_num_threads()          
#endif
      
      my_rank  = loc_my_rank
      nprocMPI = loc_nprocMPI
                      
            
!      write (str_buf,*) 'Threads allocated : ',nthreads
!      call pout(str_buf)                       
                
      ! Set random number seed to be different for each processor      
      put_seq = 10*my_thread + my_rank * nprocMPI + 37 * (/ (i - 1, i = 1, 4) /)
      put_seq = my_rank*nthreads + my_thread      + 37 * (/ (i - 1, i = 1, 4) /)
      put_seq = my_rank + my_thread*loc_nprocMPI  + 37 * (/ (i - 1, i = 1, 4) /)
      CALL rand_seed(put_seq)      
      CALL rand_number(aux_num)
      my_nchex_events  = 200000 + NINT(200000*aux_num)                   
      my_nchex_current = 0

      allocate(my_chex_events(1:my_nchex_events))     
 
                                                
    !$omp end parallel
        
    ncores   = nprocMPI*nthreads
    
    
    
!    if (my_rank == 0) then
      !call init_pout()  
      !write (str_buf,*) 'Threads allocated : ',nthreads
      !call pout(str_buf)                             
!    endif
     

    kappa = 1.63_dp
! calculate gamma(kappa+1) and gamma(kappa-1/2) for heliosheath charge-exchange
    CALL initialize_gammas
    
    kappa_v_th_factor_sq = (kappa - 1.5_dp) / kappa
    kappa_v_th_factor = SQRT( kappa_v_th_factor_sq )
    IF (my_rank .EQ. 0) THEN
       PRINT *, 'kappa = ', kappa
       PRINT *, 'gamma_vrel_factor = ', gamma_vrel_factor
       PRINT *, 'kappa v_th factor = ', kappa_v_th_factor
    END IF


    call init_vr_phi_maxwel
    call init_vr_phi_kappa    


! checking get_Maxwellian_vr(vr)
  !~if (my_rank .EQ. 0) then
  !~open(unit = 668,file='get_Maxwellian_vr.txt',form='formatted',status='REPLACE' ,action='write')
  !~do i = 1, 50000
  !~  call get_Maxwellian_vr(aux_num)
  !~  write(668,*) aux_num
  !~enddo
  !~close(668)
  !~endif
       
        
! set up global parameters
    LISM_nH = 1d+6*au_LISM_nH                                           ! 0.22e6_dp  per cubic meter
    LISM_vH = 1d-2*sqrt(au_LISM_vH(0)**2+au_LISM_vH(1)**2+au_LISM_vH(2)**2) !-23.2e3_dp  meters per second
    LISM_TH = au_LISM_TH                                                ! 6200.0_dp  degrees Kelvin
    LISM_vx = 1d-2*au_LISM_vH(0)
    LISM_vy = 1d-2*au_LISM_vH(1)
    LISM_vz = 1d-2*au_LISM_vH(2)

    k_B = 1.38e-23_dp                  ! Joules per Kelvin (m^2 kg s^{-2} K^{-1})
    neutral_mass = 1.67e-27_dp         ! Kilograms
    LISM_v_thermal = SQRT(2.0_dp*k_B*LISM_TH/neutral_mass) ! 3D thermal speed

    IF (my_rank .EQ. 0) THEN
       PRINT *,'LISM velocity components: ', LISM_vx, LISM_vy, LISM_vz       
       PRINT *,'LISM thermal velocity (m/s)'
       PRINT 9,LISM_v_thermal
    END IF
    
!    ALLOCATE(cell_ratios(0:p_nlevels-1))
    cell_ratios(0) = 1
    IF (nsplit_grids >= 1) THEN
      split_xmin = (/max(-2500_dp*au,g_xmin+300*au),-350_dp*au,-140_dp*au,-65_dp*au/)
      !split_xmax = (/ 420_dp*au,                     200_dp*au, 100_dp*au, 60_dp*au/)
      split_xmax = (/ 530_dp*au,                     500_dp*au, 100_dp*au, 60_dp*au/)

      split_ymin = (/-512_dp*au,-300_dp*au,-128_dp*au,-64_dp*au/)
      split_ymax = (/ 512_dp*au, 300_dp*au, 128_dp*au, 64_dp*au/)
      
      split_zmin = (/-512_dp*au,-300_dp*au,-128_dp*au,-64_dp*au/)
      split_zmax = (/ 512_dp*au, 300_dp*au, 128_dp*au, 64_dp*au/)
            
      !DO kk = 2, nsplit_grids 
      !  split_xmin(kk) = 0.5_dp*split_xmin(kk-1)  
      !  split_ymin(kk) = 0.5_dp*split_ymin(kk-1)
      !  split_zmin(kk) = 0.5_dp*split_zmin(kk-1)
      !  split_xmax(kk) = 0.5_dp*split_xmax(kk-1)
      !  split_ymax(kk) = 0.5_dp*split_ymax(kk-1)
      !  split_zmax(kk) = 0.5_dp*split_zmax(kk-1)
      !enddo      
    endif
    DO kk = 1, nsplit_grids 
      cell_ratios(kk) = cell_ratios(kk-1)*2.0_dp
    enddo
        
    
    m_levels(p_nlevels-1)%sampling_dt =  3*ceiling(min(m_levels(p_nlevels-1)%p_dx, &
      min(m_levels(p_nlevels-1)%p_dy,m_levels(p_nlevels-1)%p_dz))/(LISM_vH*3600.0*24.0*365.0))
          
    
    do ilev = p_nlevels - 2, 0, -1
      m_levels(ilev)%sampling_dt = m_levels(ilev+1)%sampling_dt*m_levels(ilev+1)%m_ref_ratio    
    enddo
    
    IF (my_rank .eq. 0) THEN      
      do ilev = 0, p_nlevels - 1      
        PRINT *,'level ', ilev, 'sampling_dt ',m_levels(ilev)%sampling_dt
      enddo
    endif
            
        
    total_neutrals = au_total_neutrals
    
    
! set up coarse timestep which moves LISM atoms about 1.5% of a mfp
    timestep = secs_per_year             ! one year
!    timestep = 2.0_dp*secs_per_year      ! try 2 years for long tail

! define this in terms of total_neutrals first, then update below
! divide by p_nlevels so our total remains approx, even after SPLITTING
    tot_volume_spliting = (m_levels(0)%p_xmax-m_levels(0)%p_xmin)* &
                 (m_levels(0)%p_ymax-m_levels(0)%p_ymin)* &
                 (m_levels(0)%p_zmax-m_levels(0)%p_zmin)
    IF (nsplit_grids >= 1) THEN
       DO kk = 1, nsplit_grids
          xmax = split_xmax(kk);       xmin = split_xmin(kk)
          ymax = split_ymax(kk);       ymin = split_ymin(kk)
          zmax = split_zmax(kk);       zmin = split_zmin(kk)
          volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
          tot_volume_spliting = tot_volume_spliting &
                                  - (cell_ratios(kk-1))**3 * volume &
                                  + (cell_ratios(kk))**3 * volume
       END DO
    END IF
    LISM_density = total_neutrals / tot_volume_spliting
    
    
    xmax = m_levels(0)%p_xmax;       xmin = m_levels(0)%p_xmin
    ymax = m_levels(0)%p_ymax;       ymin = m_levels(0)%p_ymin
    zmax = m_levels(0)%p_zmax;       zmin = m_levels(0)%p_zmin

! choose a boundary ring thickness so that 99.9999998% of particles with
! velocity directed inward which make it into the domain come from this layer
    boundary_thickness = 6.0_dp*(LISM_v_thermal+ABS(LISM_vH))*timestep
    boundary_npartx = NINT( LISM_density * boundary_thickness &
         * (ymax-ymin+2.0*boundary_thickness) &
         * (zmax-zmin+2.0*boundary_thickness) / ncores )
    boundary_nparty = NINT( LISM_density * boundary_thickness &
         * (xmax-xmin) * (zmax-zmin) / ncores )
    boundary_npartz = NINT( LISM_density * boundary_thickness &
         * (xmax-xmin) * (ymax-ymin) / ncores )
             

    achieved_density = ncores*(boundary_npartx+boundary_nparty+boundary_npartz) &
                      / ( boundary_thickness  * &
                      ( (ymax-ymin+2.0*boundary_thickness)    &
                      * (zmax-zmin+2.0*boundary_thickness)  + &
                        (xmax-xmin) * (zmax-zmin) + &
                        (xmax-xmin) * (ymax-ymin) ) )

! NOTE: must define macro_mass in terms of EXPECTED # of particles (npart)
    total_volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
    npart = NINT( achieved_density * total_volume,8)
    mega_number = LISM_nH*total_volume/npart ! # neutrals in megaparticle
    macro_mass = neutral_mass * mega_number           ! mass of macroparticle

    IF (my_rank .eq. 0) THEN
       PRINT *,'timestep (sec)'
       PRINT 9,timestep
       PRINT *,'# of neutrals per megaparticle'
       PRINT 9,mega_number
       PRINT *,'boundary thickness'
       PRINT 9,boundary_thickness/au
       PRINT *,'LISM number density'
       PRINT 9,LISM_density
       PRINT *,'achieved number density'
       PRINT 9,achieved_density
       PRINT *,'number of particles per processor in boundary layer'
       PRINT *,2*(boundary_npartx+boundary_nparty+boundary_npartz)
       PRINT *,'expected total number of particles without splitting'
       PRINT *,npart
    END IF

    npart = NINT( achieved_density * tot_volume_spliting,8)
    my_ntotal = npart/ncores + ncores       ! allow for rounding

! allow for some extra particles due to splitting
!    my_ntotal = 1.1*my_ntotal
    
    IF (my_rank .eq. 0) THEN
      PRINT *,'expected total number of particles with splitting'
      PRINT *,npart
      PRINT *,'allocated buffer size for particles per core'
      PRINT *,my_ntotal
    endif
    
    allocate(my_particles(0:nthreads-1))     
    do i=0, nthreads-1
      my_particles(i)%my_npart = 0
      ALLOCATE(my_particles(i)%particles__(0:PNUM-1,my_ntotal))   
      my_particles(i)%particles__ = 0.0_dp
    enddo


    IF ( au_restart > 0 ) THEN
! load particle data from existing file
       CALL load_raw(au_load_file)
       CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
    END IF
    
9 FORMAT(6e13.4)

  END SUBROUTINE mc_init_call_2
  
    
  SUBROUTINE finalize
  
    USE hdf5
    implicit none
        
    integer :: i,ilev,ierr
    
    do i=0, nthreads-1      
      DEALLOCATE(my_particles(i)%particles__)         
    enddo
    deallocate(my_particles)     
    
    !$omp parallel 
      deallocate(my_chex_events)                                                      
    !$omp end parallel
    
    do ilev = 0, p_nlevels-1    
      deallocate(m_levels(ilev)%source_m ) 
      deallocate(m_levels(ilev)%source_px) 
      deallocate(m_levels(ilev)%source_py) 
      deallocate(m_levels(ilev)%source_pz) 
      deallocate(m_levels(ilev)%source_mvsq)     
      deallocate(m_levels(ilev)%source_nchex)    

      if (aggregate_neutral_data) then 
        deallocate(m_levels(ilev)%H_dens)
        deallocate(m_levels(ilev)%H_ux)
        deallocate(m_levels(ilev)%H_uy)
        deallocate(m_levels(ilev)%H_uz)
        deallocate(m_levels(ilev)%H_temp)
      endif
                  
    enddo
    
    if (pout_unit>0) then
      close(pout_unit)
    endif
    
    deallocate(m_levels)
    
    deallocate(mpi_balancing)    
    
    deallocate(mpi_requests)    
    deallocate(mpi_statuses)
    
    !print *, 'closing fortran hdf5'  
    call h5close_f(ierr)
    !print *, 'fortran hdf5 closed'  
    
  END SUBROUTINE finalize
  
  
!
!***********************************************************************
!
  SUBROUTINE move
    REAL(KIND=dp), PARAMETER :: ideal_chex_prob = 0.02_dp
    
    REAL(KIND=dp), PARAMETER :: kappa_rmax      = 500.0_dp*au
    REAL(KIND=dp), PARAMETER :: kappa_r1        = 150.0_dp*au ! from this distance probability of kappa decreases
    REAL(KIND=dp), PARAMETER :: kappa_k         = 1.0_dp/(kappa_r1-kappa_rmax)
    REAL(KIND=dp), PARAMETER :: kappa_b         = -kappa_rmax*kappa_k
    
    REAL(KIND=dp) :: massi, xi, yi, zi, vx, vz, vy, region, ri_sq, ri
    REAL(KIND=dp) :: n_p, T_p, ux, uy, uz, region_p
    REAL(KIND=dp) :: vxn, vzn, vyn
    REAL(KIND=dp) :: v_p_th, v_p_th_sq, sigma_vrelp
    REAL(KIND=dp) :: subtime, dt, ideal_dt, sigma_ds
    REAL(KIND=dp) :: beta_ph, mass_loss, dprand
    INTEGER :: my_old_npart, my_new_npart, i, icount, use_kappa, cur_my_npart
    INTEGER :: ns, nsplit, nsplit_old
    INTEGER :: count_total_dt
    LOGICAL :: nochex
    REAL(KIND=dp), DIMENSION (:,:), POINTER :: particles

    count_total_dt = 0
    my_min_dt = timestep
    my_max_dt = 0.0_dp
    my_ave_dt = 0.0_dp
    icount = 0
    my_n_chex = 0

    my_tot_mass_loss = 0.0_dp
    my_old_npart = my_particles(my_thread)%my_npart
    my_new_npart = my_particles(my_thread)%my_npart
    i = 0
    
    particles => my_particles(my_thread)%particles__
    
    DO WHILE (i .LT. my_particles(my_thread)%my_npart)
       i = i + 1

       massi  =  particles(0,i)
       xi     =  particles(1,i)
       yi     =  particles(2,i)
       zi     =  particles(3,i)
       vx     =  particles(4,i)
       vy     =  particles(5,i)
       vz     =  particles(6,i)
       region =  particles(7,i)

       ri_sq = xi*xi + yi*yi + zi*zi
       ri = SQRT(ri_sq)
       subtime = 0.0_dp
! use adaptive timestep
       DO WHILE (subtime .LT. timestep)

! obtain local plasma quantities
          CALL plasma(xi,yi,zi,n_p,T_p,ux,uy,uz,region_p)
          v_p_th_sq = 2.0_dp*k_B*T_p/neutral_mass
          v_p_th = SQRT(v_p_th_sq)

          use_kappa = 0
          IF ( (ri .LT. kappa_rmax) .AND. &  ! don't use kappa in distant tail
               (region_p .GT. 1.5_dp).AND. &  ! solar wind plasma 
               (region   .LT. 1.5_dp)       & ! interstellar particle
                ) THEN
                
            CALL rand_number(dprand)
            
            if (dprand < (ri*kappa_k + kappa_b)) then                
              use_kappa = 1
            endif
            
          END IF
          
          use_kappa = 0
          
          CALL get_sigma_vrelp(vx,vy,vz,ux,uy,uz,v_p_th_sq,sigma_vrelp,use_kappa)

          ideal_dt = ideal_chex_prob/(n_p*sigma_vrelp)
! NOTE: might need to check that dt is small enough over entire space step
! (seems ok, even at 1 AU)

          dt = MIN(ideal_dt,timestep)

          count_total_dt = count_total_dt + 1

          IF (subtime+dt .gt. timestep) THEN
             dt = timestep - subtime
          END IF

! DIAGNOSTIC:
          my_max_dt = MAX(my_max_dt,dt)
          my_min_dt = MIN(my_min_dt,dt)

          my_ave_dt = my_ave_dt + dt
          icount = icount + 1

          nochex = .TRUE.
          sigma_ds = sigma_vrelp * dt  ! use (RELATIVE) distance travelled for Prob(chex)

          nsplit_old = NINT(my_particles(my_thread)%particles__(8,i))
          DO ns = 1, nsplit_old
          
             CALL rand_number(dprand)
                                                    
             IF ( dprand .LT. n_p*sigma_ds ) THEN   ! charge exchange occurs

                CALL exchange(massi,xi,yi,zi,vx,vy,vz, &
                           ux,uy,uz,v_p_th,vxn,vyn,vzn,use_kappa)

                nsplit = NINT(particles(8,i))
                IF (nsplit .EQ. 1) THEN

                   nochex = .FALSE.
                   vx = vxn;  vy = vyn;  vz = vzn
                   region = region_p

                ELSE

! NOTE: charge-exchanged particle goes to end of particle array,
! we continue to follow parent particle
                   particles(8,i) = nsplit - 1.0

! CREATE new particle
                   IF (my_new_npart .LT. my_ntotal) THEN
                      my_new_npart = my_new_npart + 1
                      particles(0,my_new_npart) = massi
                      particles(1,my_new_npart) = xi
                      particles(2,my_new_npart) = yi
                      particles(3,my_new_npart) = zi
                      particles(4,my_new_npart) = vxn
                      particles(5,my_new_npart) = vyn
                      particles(6,my_new_npart) = vzn
                      particles(7,my_new_npart) = region_p
                      particles(8,my_new_npart) = 1.0    ! single particle
                   ELSE
                      PRINT *, my_rank, 'my_ntotal exceeded'
                      call abort
                   END IF

                END IF

             END IF

          END DO

! don't move newly charge-exchanged particles, since we don't know the correct dt
          IF ( nochex ) THEN
             ri_sq = xi*xi + yi*yi + zi*zi
! NOTE: need to photoionize each small step for accuracy close to the Sun
! photoionization; assumes beta_ph = beta_ph_E at 1au
             beta_ph = beta_ph_E * r_E_sq / ri_sq
             mass_loss = beta_ph * dt * massi
             my_tot_mass_loss = my_tot_mass_loss + MIN(massi,mass_loss)*nsplit_old

             IF ( ri_sq .LT. phi_rmax_sq ) THEN
                CALL source_phi(xi,yi,zi,vx,vy,vz, &
                                MIN(massi,mass_loss)*nsplit_old)
             END IF
             massi = massi - mass_loss

             subtime = subtime + dt
! advance particles in a straight line at current velocity (Euler's method)
             xi = xi + vx*dt
             yi = yi + vy*dt
             zi = zi + vz*dt
          END IF

       END DO

! update particle specs
       particles(0,i) = massi
       particles(1,i) = xi
       particles(2,i) = yi
       particles(3,i) = zi
       particles(4,i) = vx
       particles(5,i) = vy
       particles(6,i) = vz
       particles(7,i) = region

! new position values
       ri_sq = xi*xi + yi*yi + zi*zi

! remove particles that have moved out of the domain
! photoionization takes care of rmin bound
       IF ( (xi .GT. g_xmax) .OR. (xi .LT. g_xmin) .OR. &
            (yi .GT. g_ymax) .OR. (yi .LT. g_ymin) .OR. &
            (zi .GT. g_zmax) .OR. (zi .LT. g_zmin) .OR. &
            (massi .LT. 0.0_dp) ) THEN
          particles(:,i) = particles(:,my_particles(my_thread)%my_npart)   ! fill gap with last particle
          my_particles(my_thread)%my_npart = my_particles(my_thread)%my_npart - 1    ! stop tracking particles outside domain
          i = i - 1                  ! don't increment if we remove a particle
       END IF

    END DO

    my_ave_dt = my_ave_dt / icount

! push newly created particles over removed particles
    cur_my_npart = my_particles(my_thread)%my_npart
    particles(:,cur_my_npart+1:cur_my_npart+my_new_npart-my_old_npart) = particles(:,my_old_npart+1:my_new_npart)
! set any remaining non-particle entries to zero
    particles(:,cur_my_npart+my_new_npart-my_old_npart+1:my_new_npart) = 0.0
    my_particles(my_thread)%my_npart = my_particles(my_thread)%my_npart+my_new_npart-my_old_npart

    IF (my_particles(my_thread)%my_npart+2*(boundary_npartx+boundary_nparty+boundary_npartz) .GT. my_ntotal) THEN
       my_particles(my_thread)%my_npart = my_particles(my_thread)%my_npart - 1
       PRINT *,'warning my_ntotal exceeded by proc ',my_rank, my_ntotal
       call abort
    ELSE
       CALL boundary_particles
    END IF

  END SUBROUTINE move
!
!***********************************************************************
!
  SUBROUTINE split
    REAL(KIND=dp) :: massi, xi, yi, zi, vxi, vyi, vzi, region, nsplit
    REAL(KIND=dp) :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(KIND=dp) :: xi1, yi1, zi1, new_massi, ideal_massi, dprand
    INTEGER :: my_old_npart, part, kk, new_nsplit
    INTEGER :: level_now, level_soon, level, ilev
    REAL(KIND=dp), DIMENSION (:,:), POINTER :: particles

    my_old_npart = my_particles(my_thread)%my_npart
    part = 0
    
    particles => my_particles(my_thread)%particles__
    
    DO WHILE (part .LT. my_particles(my_thread)%my_npart)
       part = part + 1

       massi  = particles(0,part)
       xi     = particles(1,part)
       yi     = particles(2,part)
       zi     = particles(3,part)
       vxi    = particles(4,part)
       vyi    = particles(5,part)
       vzi    = particles(6,part)
       region = particles(7,part)
       nsplit = particles(8,part)
       
       if ( (xi .GT. g_xmax) .OR. (xi .LT. g_xmin) .OR. &
            (yi .GT. g_ymax) .OR. (yi .LT. g_ymin) .OR. &
            (zi .GT. g_zmax) .OR. (zi .LT. g_zmin) .OR. &
            (massi .LT. 0.0_dp) ) THEN
        cycle
       endif
            
            

! which grid are we on NOW
       level_now = 0
              
       do ilev = 1, nsplit_grids 
        if ((split_xmin(ilev)<=xi).and.(xi<=split_xmax(ilev)).and. &
            (split_ymin(ilev)<=yi).and.(yi<=split_ymax(ilev)).and. &
            (split_zmin(ilev)<=zi).and.(zi<=split_zmax(ilev))) then          
          level_now = ilev                       
          cycle
        else
          level_now = ilev-1                       
          exit 
        endif          
      enddo
                   

! which grid will get to SOON -- using 10 year ~ 50 AU at V_LISM
       xi1 = xi + 10.0_dp*secs_per_year*vxi
       yi1 = yi + 10.0_dp*secs_per_year*vyi
       zi1 = zi + 10.0_dp*secs_per_year*vzi
       level_soon = 0              
       do ilev = 1, nsplit_grids
        if ((split_xmin(ilev)<=xi1).and.(xi1<=split_xmax(ilev)).and. &
            (split_ymin(ilev)<=yi1).and.(yi1<=split_ymax(ilev)).and. &
            (split_zmin(ilev)<=zi1).and.(zi1<=split_zmax(ilev))) then
          
          level_soon = ilev                       
          cycle
        else
          level_soon = ilev-1                       
          exit 
        endif          
      enddo
              

       level = MAX(level_now,level_soon)
              
       ideal_massi = 1.0_dp / (cell_ratios(level))**3

       IF (massi .GT. 2.0_dp*ideal_massi) THEN
! we SPLIT
          particles(0,part) = massi  * 0.5
          particles(8,part) = nsplit * 2.0

       ELSE
! make sure we don't split AND recombine

! check whether to RECOMBINE
          IF ( (massi .LT. 0.5_dp*ideal_massi) &
               .AND. (region .LT. 1.5_dp) ) THEN
! don't recombine comp2 and comp3 neutrals

             new_nsplit = FLOOR(nsplit * 0.5_dp)
             IF (new_nsplit .GT. 0) THEN
                new_massi = massi * REAL(nsplit/new_nsplit,dp)
                
                particles(0,part) = new_massi
                particles(8,part) = REAL(new_nsplit)
             ELSE
! zeroth order recombination of single particles
                CALL rand_number(dprand)
                IF (dprand .LT. 0.5_dp) THEN
! double mass
                   particles(0,part) = particles(0,part) * 2.0
                ELSE
! dissintegrate
                   particles(:,part) = particles(:,my_particles(my_thread)%my_npart)
                   my_particles(my_thread)%my_npart = my_particles(my_thread)%my_npart - 1
                END IF
             END IF

          END IF

       END IF

    END DO

    particles(:,my_particles(my_thread)%my_npart+1:my_old_npart) = 0.0

  END SUBROUTINE split


!
!***********************************************************************
!
  SUBROUTINE exchange(mass,xi,yi,zi,vx,vy,vz,ux,uy,uz,v_p_th,vxn,vyn,vzn,use_kappa)
    REAL(KIND=dp) :: mass, xi, yi, zi, vx, vy, vz, ux, uy, uz
    REAL(KIND=dp) :: v_p_th, vr, phi, dprand
    REAL(KIND=dp) :: vxn, vyn, vzn
    REAL(KIND=dp) :: vv1, vv2, vv3, vv, vvsq, modvperp, modvvnew
    REAL(KIND=dp) :: r1_1, r1_2, r1_3, modr1, r1hat_1, r1hat_2, r1hat_3
    REAL(KIND=dp) :: r2_1, r2_2, r2_3, modr2, r2hat_1, r2hat_2, r2hat_3
    REAL(KIND=dp) :: vvnew1, vvnew2, vvnew3, omega
    INTEGER :: use_kappa

    my_n_chex = my_n_chex + 1

! NOTE: MASS of new particle same as OLD

! rotate vH about z-axis to bring us onto the frame of the plasma
    vv1 = vx-ux
    vv2 = vy-uy
    vv3 = vz-uz

    vvsq = vv1*vv1+vv2*vv2+vv3*vv3
    vv = SQRT(vvsq)

    IF (use_kappa .EQ. 1) THEN
       CALL get_vr_phi_kappa(vv,v_p_th,vr,phi)
    ELSE
       CALL get_vr_phi_maxwel(vv,v_p_th,vr,phi)
    END IF

    IF (ABS(vv3) .GT. ABS(vv2)) THEN
       r1_1 = vv1
       r1_2 = vv2
       r1_3 = vv3 - vvsq/vv3
    ELSEIF (ABS(vv1) .GT. ABS(vv2)) THEN
       r1_1 = vv1 - vvsq/vv1
       r1_2 = vv2
       r1_3 = vv3
    ELSE
       r1_1 = vv1
       r1_2 = vv2 - vvsq/vv2
       r1_3 = vv3
    END IF
    modr1 = SQRT(r1_1*r1_1 + r1_2*r1_2 + r1_3*r1_3)
    r1hat_1 = r1_1/modr1
    r1hat_2 = r1_2/modr1
    r1hat_3 = r1_3/modr1

    r2_1 = r1_2*vv3 - r1_3*vv2
    r2_2 = r1_3*vv1 - r1_1*vv3
    r2_3 = r1_1*vv2 - r1_2*vv1
    modr2 = SQRT(r2_1*r2_1 + r2_2*r2_2 + r2_3*r2_3)
    r2hat_1 = r2_1/modr2
    r2hat_2 = r2_2/modr2
    r2hat_3 = r2_3/modr2
    
    CALL rand_number(dprand)
    omega = 2.0_dp*pi * dprand             ! random phase for scatter by phi

    modvperp = vv*TAN(phi)
    IF (phi < 0.5_dp*pi) THEN
       vvnew1 =  vv1 + modvperp*(r1hat_1*COS(omega) + r2hat_1*SIN(omega))
       vvnew2 =  vv2 + modvperp*(r1hat_2*COS(omega) + r2hat_2*SIN(omega))
       vvnew3 =  vv3 + modvperp*(r1hat_3*COS(omega) + r2hat_3*SIN(omega))
    ELSE
       vvnew1 = -vv1 + modvperp*(r1hat_1*COS(omega) + r2hat_1*SIN(omega))
       vvnew2 = -vv2 + modvperp*(r1hat_2*COS(omega) + r2hat_2*SIN(omega))
       vvnew3 = -vv3 + modvperp*(r1hat_3*COS(omega) + r2hat_3*SIN(omega))
    END IF
    modvvnew = SQRT(vvnew1*vvnew1+vvnew2*vvnew2+vvnew3*vvnew3)

    vxn = ux + vvnew1 * vr/modvvnew
    vyn = uy + vvnew2 * vr/modvvnew
    vzn = uz + vvnew3 * vr/modvvnew

! particle position unchanged

! add ch-ex event to plasma source term
    
    CALL source_chex(mass,xi,yi,zi,vx,vy,vz,vxn,vyn,vzn)    

  END SUBROUTINE exchange
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
    integer,       intent(in) :: bnpart
    REAL(KIND=dp), intent(in) :: buffer(0:PNUM-1,1:bnpart)
    
    
    integer :: i, sP, iSenShift, avail, icomp
    
    sP        = bnpart
    iSenShift = 1
        
    do i = 0, nthreads - 1
    
      avail = my_ntotal - my_particles(i)%my_npart
      
      if (avail <= 0) then
        cycle
      endif
      
      if (sP<=avail) then
        
        my_particles(i)%particles__(:,my_particles(i)%my_npart+1:my_particles(i)%my_npart+sP) = & 
                buffer(:,iSenShift:iSenShift+sP-1)        
        
        !my_particles(i)%particles(my_particles(i)%my_npart+1:my_particles(i)%my_npart+sP,:) = & 
        !        buffer(iSenShift:iSenShift+sP-1,:)        
        
        my_particles(i)%my_npart = my_particles(i)%my_npart + sP
        exit
      else
        
        my_particles(i)%particles__(:,my_particles(i)%my_npart+1:my_ntotal) = & 
            buffer(:,iSenShift:iSenShift+avail-1)        
        
        !my_particles(i)%particles(my_particles(i)%my_npart+1:my_ntotal,:) = buffer(iSenShift:iSenShift+avail-1,:)        
                
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
    
    CHARACTER (LEN = 100) str_buf
    
    
    
    if (nthreads == 1) then
      return
    endif
  
    my_npart = 0
    do ithread = 0, nthreads - 1
      my_npart = my_npart + my_particles(ithread)%my_npart
    enddo      
                
    npart_average = my_npart/nthreads + 1    
    
    !write (str_buf,*) 'balance_node/particles on node: ',my_npart,'average: ',npart_average 
    !call pout(str_buf)
            
    do ithread = 0, nthreads - 1
      node_balancing(ithread)%procID = ithread;
      node_balancing(ithread)%part_disbalance = my_particles(ithread)%my_npart - npart_average   

      !write (str_buf,*) ' thread ',ithread, 'particles ',my_particles(ithread)%my_npart,'disbalance ', & 
      !    node_balancing(ithread)%part_disbalance
      !call pout(str_buf)
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
              
        
        if (sP < rP) then          
          
          my_particles(iR)%particles__(:,iRecShift:iRecShift+sP-1) = my_particles(iS)%particles__(:,iSenShift-sP+1:iSenShift)
          
          
          !print *,'thread',iS,'sends [',iSenShift-sP+1,iSenShift,'] to ',iR,' [',iRecShift,iRecShift+sP-1,']'
          
          my_particles(iS)%my_npart = my_particles(iS)%my_npart-sP
          my_particles(iR)%my_npart = my_particles(iR)%my_npart+sP
          
          iRecShift  = my_particles(iR)%my_npart+1
                    
          rP = rP - sP            
          
          exit
        else                  
            
          my_particles(iR)%particles__(:,iRecShift:iRecShift+rP-1) = my_particles(iS)%particles__(:,iSenShift-rP+1:iSenShift)
          
          !print *,'thread',iS,'sends [',iSenShift-rP+1,iSenShift,'] to ',iR,' [',iRecShift,iRecShift+rP-1,']'
            

          my_particles(iS)%my_npart = my_particles(iS)%my_npart-rP
          my_particles(iR)%my_npart = my_particles(iR)%my_npart+rP
          
          sP        = sP - rP
                    
                  
          iSenShift = iSenShift - rP
          jthread   = jthread  - 1
          iR         = node_balancing(jthread)%procID
          rP         = min(abs(node_balancing(jthread)%part_disbalance), my_ntotal - my_particles(iR)%my_npart)
          iRecShift  = my_particles(iR)%my_npart+1 
          
          if (sP.eq.0) then
            exit
          endif
          
        endif        
      enddo  
      
      if (node_balancing(jthread)%part_disbalance >= 0)  then
        exit
      endif
        
    enddo
    
    !~my_npart = 0
    !~do ithread = 0, nthreads - 1      
    !~  my_npart = my_npart + my_particles(ithread)%my_npart
    !~  write (str_buf,*) 'thread ',ithread, 'particles ',my_particles(ithread)%my_npart
    !~  call pout(str_buf)
    !~enddo
    !~write (str_buf,*),'particles on node: ',my_npart
    !~call pout(str_buf)
     
  end subroutine balance_node
  
 SUBROUTINE balance_npart(a_loop, a_next_balance)    
 
    integer, intent(in)  :: a_loop
    integer, intent(out) :: a_next_balance
    
    real, parameter :: max_ratio = 0.1d0   
    
    integer, parameter :: max_balance_interval = 200
    integer, parameter :: min_balance_interval = 5
    
    integer, parameter :: hist_bins = 10
    integer, dimension(0:hist_bins) :: histogram
    
    
    real(KIND=dp) :: npart_ratio,max_min_ratio,time_disbalance
    
 
    INTEGER :: tmp, min_rank, max_rank, tot_min, tot_max, ierr
    
    integer(8) :: my_npart ! number of particles on node
    integer(8) :: total_npart,npart_average, min_npart, max_npart
    integer(8), DIMENSION (0:nprocMPI-1)    :: npart
    
    type(particle_storage), pointer :: thread_particles
    
    integer :: my_disbalance, nrecv_part
    integer :: iproc,jproc,i,icomp
    integer :: iSenShift,sP,iRecShift,rP,iR,iS
    integer :: size_p
    
    REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: buffer    
    
    real(KIND=dp), DIMENSION (0:nprocMPI-1) :: timings
    integer       :: min_npart_rank, max_npart_rank,min_time_rank, max_time_rank
    real(KIND=dp) :: min_time,max_time
        
    
    
    CHARACTER (LEN = 1000) str_buf, str_buf_big
    
    INTEGER :: time1(8), time2(8), time3(8)
    real(8) :: secs1, secs2, secs3
    
    
    INTEGER :: nrequests    
    
    call date_and_time(VALUES=time1)    
    
        
    my_npart = 0
    do i = 0, nthreads - 1
      my_npart = my_npart + my_particles(i)%my_npart
    enddo      
        
    
    call MPI_allgather(my_npart, 1, MPI_LONG_LONG_INT, npart, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD, ierr)   
        
    call MPI_gather(last_time_step_sec,1,MPI_DOUBLE_PRECISION, timings, 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
               
        
    !write (str_buf,*) '             , my particles before: ', my_npart
    !call pout(str_buf) 
    
    min_npart_rank = 0; max_npart_rank = 0; min_time_rank = 0; max_time_rank = 0
    
    total_npart = 0
    min_npart   = huge(1)
    min_time    = huge(0.0_dp)
    max_npart   = 0    
    max_time    = 0.0_dp
    
    do iproc = 0, nprocMPI - 1
      total_npart = total_npart + npart(iproc)
      if (npart(iproc) < min_npart) then
        min_npart      = npart(iproc)
        min_npart_rank = iproc
      endif
      if (npart(iproc) > max_npart) then
        max_npart      = npart(iproc)
        max_npart_rank = iproc
      endif
      if (timings(iproc) < min_time) then
        min_time      = timings(iproc)
        min_time_rank = iproc
      endif
      if (timings(iproc) > max_time) then
        max_time      = timings(iproc)
        max_time_rank = iproc
      endif
      !min_npart   = min(min_npart,npart(iproc))
      !max_npart   = max(max_npart,npart(iproc))            
    enddo        
    
    npart_average = total_npart/nprocMPI + 1
    
    max_min_ratio = real(max_npart)/real(max(min_npart,1))
    
    !write (str_buf,'(a,f8.3,f8.3,f8.3,i6,i6)') 'timings (min/max, ratio, ranks)', & 
    !  min_time,max_time,max_time/min_time,min_time_rank,max_time_rank
    
!    write (str_buf_big,*) ' min/average/max number of particles ', min_npart,npart_average,max_npart, &
!       'max/min ratio =', max_min_ratio, NEW_LINE('A'), &
!       'min/max particles rank',min_npart_rank,max_npart_rank, NEW_LINE('A'), &        
!       'timings (min/max, ratio, ranks)',  min_time,max_time,max_time/min_time,min_time_rank,max_time_rank, &
!       NEW_LINE('A'),'their particles min/max ',npart(min_time_rank)/nthreads, npart(max_time_rank)/nthreads          
              
!    call pout(str_buf_big) 
    
    if (my_rank .eq. 0) then
      histogram = 0
      time_disbalance = max_time-min_time
      if (time_disbalance>1e-3) then
      do iproc = 0, nprocMPI - 1
        tmp  = nint(hist_bins*(timings(iproc)-min_time)/time_disbalance)
        
        histogram(tmp) = histogram(tmp)+1        
      enddo
      endif
      
!      write (str_buf_big,'(a,11i5)') ' histogram  ', histogram                     
!      call pout(str_buf_big)       
    endif
    
    
        
    if     (max_min_ratio < 1.05_dp) then    
      balance_interval = nint(balance_interval*1.5_dp)
    else if ((max_min_ratio > 1.1_dp).and.(max_min_ratio < 1.2_dp)) then
      balance_interval = nint(balance_interval/1.2_dp)
    else if (max_min_ratio > 1.2_dp)  then
      balance_interval = nint(balance_interval/2.0_dp)
    endif
    
    balance_interval = min(max(balance_interval,min_balance_interval),max_balance_interval)
    a_next_balance = a_loop + balance_interval
      
        
    
    total_npart_verbosity = total_npart
        
    do iproc = 0, nprocMPI - 1
      mpi_balancing(iproc)%procID = iproc;
      mpi_balancing(iproc)%part_disbalance = npart(iproc) - npart_average            
    enddo
    
    my_disbalance     = mpi_balancing(my_rank)%part_disbalance
    
!    str_buf = ' '
!    if (my_disbalance > 0) then
!      write (str_buf,*) 'I will send ', my_disbalance, ' particles '
!    elseif (my_disbalance < 0) then   
!      write (str_buf,*) 'I can accept ', abs(my_disbalance), ' particles '
!    endif
!    call pout(str_buf)
        
    
    ! allocate sender/receiver buffer. Attention, reverse order ! Need it to correctly send linearized array
    if (my_disbalance .ne.  0) then
      allocate(buffer(0:PNUM-1,1:abs(my_disbalance)))
      buffer = 0d0
    endif
    
    npart_ratio = real(abs(my_disbalance))/real(npart_average)
      
    !if ((my_disbalance > 0).and.(npart_ratio>max_ratio)) then
    if (my_disbalance > 0) then
      sP        = my_disbalance
      iSenShift = 1
      
      do i = 0, nthreads - 1
        if (sP<=my_particles(i)%my_npart) then
          
          buffer(:,iSenShift:iSenShift+sP-1) =  &
              my_particles(i)%particles__(:,my_particles(i)%my_npart-sP+1:my_particles(i)%my_npart)
          
          !buffer(iSenShift:iSenShift+sP-1,:) = my_particles(i)%particles(my_particles(i)%my_npart-sP+1:my_particles(i)%my_npart,:)
          
          my_particles(i)%my_npart = my_particles(i)%my_npart - sP
          exit
        else
          
          buffer(:,iSenShift:iSenShift+my_particles(i)%my_npart-1) = &
              my_particles(i)%particles__(:,1:my_particles(i)%my_npart)
          
          !buffer(iSenShift:iSenShift+my_particles(i)%my_npart-1,:) = my_particles(i)%particles(1:my_particles(i)%my_npart,:)
          
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
    
    !allocate(mpi_requests(0:nprocMPI-1))
    nrequests = 0
    
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
        
        if (sP < rP) then
          if (my_rank == iS) then        
!            CALL MPI_SEND(buffer(iSenShift:iSenShift+sP-1,:), &
!                     PNUM*sP,MPI_DOUBLE_PRECISION,iR, &
!                     0,MPI_COMM_WORLD,ierr)

!            write (str_buf,*) 'a) initiate transfer to ', iR, ', particles = ',sP, ',ierr = ', ierr
!            call pout(str_buf) 


            CALL MPI_ISEND(buffer(0,iSenShift), &
                     PNUM*sP,MPI_DOUBLE_PRECISION,iR, &
                     0,MPI_COMM_WORLD,mpi_requests(nrequests),ierr)
                     
            
            nrequests = nrequests + 1
                                                                        
          elseif (my_rank == iR) then  
!            CALL MPI_RECV(buffer(iRecShift:iRecShift+sP-1,:), &
!                    PNUM*sP,MPI_DOUBLE_PRECISION,iS, &
!                    0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)


!            write (str_buf,*) 'a) initiate recieve from ', iS, ', particles = ',sP, ',ierr = ', ierr
!            call pout(str_buf) 
                    
            CALL MPI_IRECV(buffer(0,iRecShift), & ! buffer(iRecShift,0)
                    PNUM*sP,MPI_DOUBLE_PRECISION,iS, &
                    0,MPI_COMM_WORLD,mpi_requests(nrequests),ierr)
                    
             
            nrequests = nrequests + 1       
                    
            nrecv_part = nrecv_part + sP 
          endif
          rP = rP - sP
          iRecShift = iRecShift + sP
          exit
        else        
          if (my_rank == iS) then        
!            CALL MPI_SEND(buffer(iSenShift:iSenShift+rP-1,:), &
!                     PNUM*rP,MPI_DOUBLE_PRECISION,iR, &
!                     0,MPI_COMM_WORLD,ierr)

 
!           write (str_buf,*) 'b) initiate transfer to ', iR, ', particles = ',rP, ',ierr = ', ierr
!           call pout(str_buf) 
            
            CALL MPI_ISEND(buffer(0,iSenShift), &
                     PNUM*rP,MPI_DOUBLE_PRECISION,iR, &
                     0,MPI_COMM_WORLD,mpi_requests(nrequests),ierr)
                     
                     
            nrequests = nrequests + 1
            
          elseif (my_rank == iR) then  
!            CALL MPI_RECV(buffer(iRecShift:iRecShift+rP-1,:), &
!                    PNUM*rP,MPI_DOUBLE_PRECISION,iS, &
!                    0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)


!            write (str_buf,*) 'b) initiate recieve from ', iS, ', particles = ',rP, ',ierr = ', ierr
!            call pout(str_buf)      
                    
            CALL MPI_IRECV(buffer(0,iRecShift), & 
                    PNUM*rP,MPI_DOUBLE_PRECISION,iS, &
                    0,MPI_COMM_WORLD,mpi_requests(nrequests),ierr)   

                    
            nrequests = nrequests + 1
            
            nrecv_part = nrecv_part + rP 
          endif
          sP        = sP - rP
          
          iSenShift = iSenShift + rP
          jproc     = jproc - 1
          iR        = mpi_balancing(jproc)%procID 
          rP        = abs(mpi_balancing(jproc)%part_disbalance)
          iRecShift = 1         
          
          if (sP.eq.0) then
            exit
          endif
        endif
        
      enddo      
          
    enddo
    
 !   write (str_buf,*) ' nrequests =', nrequests, ', disbalance = ',my_disbalance
 !   call pout(str_buf) 
    
    if (nrequests>0) then
  !    write (str_buf,*) 'call MPI_WAITALL'
  !    call pout(str_buf) 
      CALL MPI_WAITALL(nrequests,mpi_requests,mpi_statuses,ierr)   
   !   write (str_buf,*) ' MPI_WAITALL  ierr = ', ierr
   !   call pout(str_buf) 
   !   do i=0,nrequests-1
   !     write (str_buf,*) ' request ', i, 'status (source, tag, error) ', &
   !         mpi_statuses(MPI_SOURCE,i), mpi_statuses(MPI_TAG,i), mpi_statuses(MPI_ERROR,i)
   !     call pout(str_buf) 
   !   enddo
      
    endif
        
    
    if (nrecv_part > 0) then
      call add_particles_node(buffer,nrecv_part)
    endif
    
    if (my_disbalance .ne. 0) then      
      deallocate(buffer)
    endif
    
    call date_and_time(VALUES=time2)
    
    call balance_node()
    
    call date_and_time(VALUES=time3)
    
    my_npart = 0
    do i = 0, nthreads - 1
      my_npart = my_npart + my_particles(i)%my_npart
    enddo      
    
    secs1 = secs(time1)
    secs2 = secs(time2)
    secs3 = secs(time3)
    write (str_buf,'(a,i11,a,f9.3,a,f9.3,a,f9.3,a,a,i4,a)') 'balance_npart, particles : ', total_npart, &
     ', time: ',secs3-secs1, ' MPI: ', secs2-secs1, ', balance_node ', secs3-secs2, &
     NEW_LINE('A'),'             next load balance in  ', balance_interval, ' steps'       
    call pout(str_buf) 
    
 
  END SUBROUTINE balance_npart

 
!
!***********************************************************************
!
  SUBROUTINE output_raw_threads(rawfile)
    USE HDF5
    CHARACTER (LEN = lenstr)  :: rawfile
    
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: rank, ierr, i
    integer(8) :: my_npart,total_npart
    integer(8), DIMENSION (0:nprocMPI-1) :: npart    
    
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
           
    
    !CALL h5acreate_f(gr_id, 'LISM_v_thermal', H5T_NATIVE_DOUBLE, aspace_id, attr_id, ierr)
    !CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, LISM_v_thermal, adims, ierr) 
    !CALL h5aclose_f(attr_id, ierr)
    
    !CALL h5acreate_f(gr_id, 'time', H5T_NATIVE_DOUBLE, aspace_id, attr_id, ierr)
    !CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, time, adims, ierr) 
    !CALL h5aclose_f(attr_id, ierr)
    
    CALL h5sclose_f(aspace_id, ierr)
    call h5gclose_f(gr_id, ierr)
    
    
    dimsf(1) = PNUM
    dimsf(2) = total_npart
    
        
    
    ! Create the data space for the  dataset. 
    CALL h5screate_simple_f(2, dimsf, filespace, ierr)      
      ! Create the dataset with default properties.      
    CALL h5dcreate_f(file_id, 'particles', H5T_NATIVE_DOUBLE, filespace, dset_id, ierr)
    CALL h5sclose_f(filespace, ierr)
    
    ! particle components    
    my_count(1)  = dimsf(1)
    my_offset(1) = 0
    
    my_offset(2) = 0
    if (my_rank>0) then
      my_offset(2) = sum(npart(0:my_rank-1))
    endif
    
    ! Create property list for collective dataset write      
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
    
    !print *,'output_raw, total number of particles',total_npart
    
    do i = 0, nthreads - 1
      my_count(2) = my_particles(i)%my_npart
      if (my_count(2).eq.0) then
        cycle
      endif
      
      !print *,'         thread',i,'number of particles ', my_count(1)
      
      
      
      !allocate(tmp_particles(1:my_count(1),0:8))
      !tmp_particles(1:my_count(1),:) = my_particles(i)%particles(1:my_count(1),:)
      
      
      
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
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, my_particles(i)%particles__(:,1:my_count(1)), dimsf, ierr, &
                     file_space_id = filespace, &
                     mem_space_id = memspace)
                     
      !CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tmp_particles, dimsf, ierr, &
      !               file_space_id = filespace, &
      !               mem_space_id = memspace, xfer_prp = plist_id)
                            
          
      ! Close dataspaces.     
      CALL h5sclose_f(filespace, ierr)
      CALL h5sclose_f(memspace, ierr)            
      
      !deallocate(tmp_particles)
      
      my_offset(2) = my_offset(2) + my_particles(i)%my_npart
    
    enddo
        
    
    ! Close the dataset.      
    CALL h5dclose_f(dset_id, ierr)    
    
    ! Close property list.      
    CALL h5pclose_f(plist_id,  ierr)
    
    CALL h5fclose_f(file_id,  ierr)
    
    !print *, 'file ',rawfile, ' is sucusssfully closed'
    
   
  END SUBROUTINE output_raw_threads
  
  
  SUBROUTINE output_raw(rawfile)
    USE HDF5
    CHARACTER (LEN = lenstr)  :: rawfile
    
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: rank, ierr, i
    integer(8) :: my_npart,total_npart
    integer(8), DIMENSION (0:nprocMPI-1) :: npart    
    
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
    
    INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
    INTEGER(HSIZE_T), DIMENSION(2) :: my_count  
    INTEGER(HSIZE_T), DIMENSION(2) :: my_offset             
    
    REAL(KIND=dp), DIMENSION (:,:), ALLOCATABLE :: tmp_particles
    
    INTEGER :: ind
    
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
           
    
    !CALL h5acreate_f(gr_id, 'LISM_v_thermal', H5T_NATIVE_DOUBLE, aspace_id, attr_id, ierr)
    !CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, LISM_v_thermal, adims, ierr) 
    !CALL h5aclose_f(attr_id, ierr)
    
    !CALL h5acreate_f(gr_id, 'time', H5T_NATIVE_DOUBLE, aspace_id, attr_id, ierr)
    !CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, time, adims, ierr) 
    !CALL h5aclose_f(attr_id, ierr)
    
    CALL h5sclose_f(aspace_id, ierr)
    call h5gclose_f(gr_id, ierr)
    
    
    dimsf(1) = PNUM
    dimsf(2) = total_npart
    
        
    
    ! Create the data space for the  dataset. 
    CALL h5screate_simple_f(2, dimsf, filespace, ierr)      
      ! Create the dataset with default properties.      
    CALL h5dcreate_f(file_id, 'particles', H5T_NATIVE_DOUBLE, filespace, dset_id, ierr)
    CALL h5sclose_f(filespace, ierr)
        
        
    my_offset(1) = 0
    
    my_offset(2) = 0
    if (my_rank>0) then
      my_offset(2) = sum(npart(0:my_rank-1))
    endif
    
    my_count(1)  = dimsf(1)
    
    my_count(2)  = 0
    do i = 0, nthreads - 1
      my_count(2) = my_count(2) + my_particles(i)%my_npart
    enddo
    
        
    ! Create property list for collective dataset write      
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
    
    !print *,'output_raw, total number of particles',total_npart
    
    if (my_count(1)>0) then
      allocate(tmp_particles(0:dimsf(1)-1,1:my_count(2)))
      ind = 1
      do i = 0, nthreads - 1
        tmp_particles(:,ind:ind+my_particles(i)%my_npart-1) = &
          my_particles(i)%particles__(:,1:my_particles(i)%my_npart)
        ind = ind + my_particles(i)%my_npart
      enddo      
      
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
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tmp_particles, dimsf, ierr, &
                     file_space_id = filespace, &
                     mem_space_id = memspace)
                     
      !CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tmp_particles, dimsf, ierr, &
      !               file_space_id = filespace, &
      !               mem_space_id = memspace, xfer_prp = plist_id)
                            
          
      ! Close dataspaces.     
      CALL h5sclose_f(filespace, ierr)
      CALL h5sclose_f(memspace, ierr)            
      
      deallocate(tmp_particles)
                
    endif
        
    
    ! Close the dataset.      
    CALL h5dclose_f(dset_id, ierr)    
    
    ! Close property list.      
    CALL h5pclose_f(plist_id,  ierr)
    
    CALL h5fclose_f(file_id,  ierr)
    
    !print *, 'file ',rawfile, ' is sucusssfully closed'
    
   
  END SUBROUTINE output_raw
  
!
!***********************************************************************
!
  SUBROUTINE load_raw(au_load_file)
    use hdf5
    CHARACTER (LEN = lenstr)  :: au_load_file
    
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
    
    CALL h5fopen_f(au_load_file, H5F_ACC_RDONLY_F, file_id, ierr)

    call h5gopen_f(file_id, '/', gr_id, ierr)

    !CALL h5aopen_name_f(gr_id, 'time', attr_id, ierr)
    !CALL h5aread_f(attr_id, H5T_NATIVE_DOUBLE, time, adims, ierr) 
    !CALL h5aclose_f(attr_id, ierr)

    call h5gclose_f(gr_id, ierr)
    
    CALL h5dopen_f(file_id, 'particles', dset_id, ierr)            
    CALL h5dget_space_f(dset_id, filespace, ierr)
    call H5Sget_simple_extent_dims_f(filespace, dimsf, maxdimsf, ierr)
    
    load_total_npart = dimsf(2)
    new_npart = load_total_npart/(nprocMPI*nthreads)
    
    my_count(1)  = dimsf(1)
    my_offset(1) = 0
    
    
    my_offset(2) = my_rank*nthreads*new_npart
    
    ! Create property list for collective dataset write      
    !CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr) 
    !CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    
    do i = 0, nthreads - 1
      my_count(2) = new_npart
      if ((my_rank.eq.(nprocMPI-1)).and.(i.eq.(nthreads - 1))) then
        my_count(2) = load_total_npart - (nthreads*nprocMPI - 1)*new_npart 
        !my_count(1) = load_total_npart - (nthreads*(nprocMPI - 1)*new_npart + (nthreads - 1)*new_npart)
      endif
      
      !allocate(tmp_particles(1:my_count(1),0:7))
      !tmp_particles(:,:) = my_particles(i)%particles(1:my_count(1),0:8)
            
      
      !write (str_buf,*) 'read [',my_offset(1),'..',my_offset(1)+my_count(1) - 1,'] particles'
      !call pout(str_buf)      
      
      !print *, 'read [',my_offset(1),'..',my_offset(1)+my_count(1)-1,'] particles'
      
      CALL h5screate_simple_f(2, my_count, memspace, ierr) 
             
      ! Select hyperslab in the file.      
      CALL h5dget_space_f(dset_id, filespace, ierr)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,  my_offset, my_count, ierr)
                             
      
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, my_particles(i)%particles__, dimsf, ierr, &
                     file_space_id = filespace, &
                     mem_space_id = memspace)
                     
      !CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, my_particles(i)%particles(1:my_count(1),0:dimsf(2)-1), dimsf, ierr, &
      !               file_space_id = filespace, &
      !               mem_space_id = memspace, xfer_prp = plist_id)
                     
                            
          
      ! Close dataspaces.     
      CALL h5sclose_f(filespace, ierr)
      CALL h5sclose_f(memspace, ierr)            
      
      my_offset(2) = my_offset(2) + new_npart
      
      my_particles(i)%my_npart = my_count(2)
      
      !deallocate(tmp_particles)
    
    enddo
    
    ! Close property list.      
    !CALL h5pclose_f(plist_id,  ierr)
    
    ! Close the dataset.      
    CALL h5dclose_f(dset_id, ierr)    
    
    CALL h5fclose_f(file_id,  ierr)
    
    
  END SUBROUTINE load_raw
  
  
  
  
  
!
!***********************************************************************
!
  SUBROUTINE output_source(rtime)
  
  INTEGER :: i,j,k,ilev,ierr,buf_size
  REAL(KIND=dp) :: volume,rtime, factor
  
  REAL(KIND=dp), DIMENSION (:,:,:), POINTER :: source_m,source_px,source_py,source_pz,source_mvsq
  
    
  do ilev = 0, p_nlevels - 1
    volume = m_levels(ilev)%p_dx*m_levels(ilev)%p_dy*m_levels(ilev)%p_dz
        
    factor = macro_mass / volume / rtime
       
      
  
    m_levels(ilev)%source_m    = m_levels(ilev)%source_m * factor
  
    
  
    m_levels(ilev)%source_px   = m_levels(ilev)%source_px* factor
  
    
  
    m_levels(ilev)%source_py   = m_levels(ilev)%source_py* factor
  
    
  
    m_levels(ilev)%source_pz   = m_levels(ilev)%source_pz* factor
  
    
  
    m_levels(ilev)%source_mvsq = m_levels(ilev)%source_mvsq*factor*0.5_dp
  
  
  enddo
  
  !~!$OMP PARALLEL default(none) private(volume,factor,ilev) shared(m_levels,macro_mass,rtime,p_nlevels)
  !~!$OMP single
  !~do ilev = 0, p_nlevels - 1
  !~  volume = m_levels(ilev)%p_dx*m_levels(ilev)%p_dy*m_levels(ilev)%p_dz
  !~      
  !~  factor = macro_mass / volume / rtime
  !~     
  !~    
  !~  !$OMP task default(none) firstprivate(ilev,factor) shared(m_levels) 
  !~  m_levels(ilev)%source_m    = m_levels(ilev)%source_m * factor
  !~  !$omp end task
  !~  
  !~  !$OMP task default(none) firstprivate(ilev,factor) shared(m_levels) 
  !~  m_levels(ilev)%source_px   = m_levels(ilev)%source_px* factor
  !~  !$omp end task
  !~  
  !~  !$OMP task default(none) firstprivate(ilev,factor) shared(m_levels) 
  !~  m_levels(ilev)%source_py   = m_levels(ilev)%source_py* factor
  !~  !$omp end task
  !~  
  !~  !$OMP task default(none) firstprivate(ilev,factor) shared(m_levels) 
  !~  m_levels(ilev)%source_pz   = m_levels(ilev)%source_pz* factor
  !~  !$omp end task
  !~  
  !~  !$OMP task default(none) firstprivate(ilev,factor) shared(m_levels) 
  !~  m_levels(ilev)%source_mvsq = m_levels(ilev)%source_mvsq*factor*0.5_dp
  !~  !$omp end task
  !~
  !~enddo
  !~!$OMP end single
  !~!$OMP end PARALLEL
  
  
  
  
   !~source_m    => m_levels(ilev)%source_m 
    !~source_px   => m_levels(ilev)%source_px
    !~source_py   => m_levels(ilev)%source_py
    !~source_pz   => m_levels(ilev)%source_pz
    !~source_mvsq => m_levels(ilev)%source_mvsq
    
   !~!$OMP PARALLEL WORKSHARE default(none) shared(factor,source_m,source_px,source_py,source_pz,source_mvsq)
   !~ source_m    = source_m * factor
   !~ source_px   = source_px* factor
   !~ source_py   = source_py* factor
   !~ source_pz   = source_pz* factor
   !~ source_mvsq = source_mvsq*factor*0.5_dp
   !~!$OMP END PARALLEL WORKSHARE
  
  !~!$OMP PARALLEL default(none) shared(factor,source_m,source_px,source_py,source_pz,source_mvsq)
  !~  !$OMP sections
  !~  !$OMP section
  !~  source_m    = source_m * factor
  !~  !$OMP section
  !~  source_px   = source_px* factor
  !~  !$OMP section
  !~  source_py   = source_py* factor
  !~  !$OMP section
  !~  source_pz   = source_pz* factor
  !~  !$OMP section
  !~  source_mvsq = source_mvsq*factor*0.5_dp
  !~ !$OMP END sections
  !~ !$OMP END PARALLEL 
  

  do ilev = 0, p_nlevels - 1  
  
    buf_size = m_levels(ilev)%p_nx*m_levels(ilev)%p_ny*m_levels(ilev)%p_nz
    
    CALL MPI_AllReduce(m_levels(ilev)%source_m,  m_levels(ilev)%g_source_m, buf_size,&
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         
    CALL MPI_AllReduce(m_levels(ilev)%source_px,  m_levels(ilev)%g_source_px, buf_size,&
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    CALL MPI_AllReduce(m_levels(ilev)%source_py,  m_levels(ilev)%g_source_py, buf_size,&
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
                  
    CALL MPI_AllReduce(m_levels(ilev)%source_pz,  m_levels(ilev)%g_source_pz, buf_size,&
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
                  
    CALL MPI_AllReduce(m_levels(ilev)%source_mvsq,  m_levels(ilev)%g_source_mvsq, buf_size,&
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         
     
    buf_size = (m_levels(ilev)%p_nx+2)*(m_levels(ilev)%p_ny+2)*(m_levels(ilev)%p_nz+2)
              
    CALL MPI_AllReduce(m_levels(ilev)%source_nchex,  m_levels(ilev)%g_source_nchex, buf_size,&
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    
  enddo
  
  
  END SUBROUTINE output_source
  
  
  SUBROUTINE output_grid
  
    CHARACTER (LEN = 200) :: str_buf
    INTEGER :: ilev
    REAL(KIND=dp) :: coeff
    
    if (aggregate_neutral_data) then 
    
      do ilev = 0, p_nlevels - 1
      
!        write (str_buf,*) 'level ', ilev, ' number of samples', m_levels(ilev)%num_samples
!        call pout(str_buf)
        
        coeff = 1.0/m_levels(ilev)%num_samples
        
        m_levels(ilev)%g_H_dens = coeff*m_levels(ilev)%g_H_dens
        m_levels(ilev)%g_H_ux   = coeff*m_levels(ilev)%g_H_ux
        m_levels(ilev)%g_H_uy   = coeff*m_levels(ilev)%g_H_uy
        m_levels(ilev)%g_H_uz   = coeff*m_levels(ilev)%g_H_uz
        m_levels(ilev)%g_H_temp = coeff*m_levels(ilev)%g_H_temp
        
      enddo
      
    else
    
      do ilev = 0, p_nlevels - 1      
        call collect_neutral_distr(ilev)               
      enddo
    
    endif
  
  END SUBROUTINE output_grid
  
!
!***********************************************************************
!
  SUBROUTINE collect_neutral_distr(a_level)
    INTEGER, INTENT(IN) :: a_level
    REAL(KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: tmp    
    REAL(KIND=dp) :: mass, xi, yi, zi, vx, vy, vz, vrel_sq
    REAL(KIND=dp) :: dx, dy, dz, x1, x2, y1, y2, z1, z2, dxdydz
    REAL(KIND=dp) :: dx1, dx2, dy1, dy2, dz1, dz2
    REAL(KIND=dp) :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(KIND=dp) :: xmin1, ymin1, zmin1
    INTEGER :: level, i, kk, nx, ny, nz, ix, iy, iz, ierr, ithread
    REAL(KIND=dp), DIMENSION (:,:,:),  POINTER :: H_dens, H_ux, H_uy, H_uz, H_temp
    REAL(KIND=dp), DIMENSION (:,:),    POINTER :: particles
    
    CHARACTER (LEN = 200) str_buf
        

! MORE EFFICIENT TO LOOP OVER PARTICLES ONCE, BUT THIS WAY IS EASIER FOR TESTING....
    DO level = a_level, a_level   
    
       !write (str_buf,*) 'output_grid, level ', level
       !call pout(str_buf)   

       xmin = m_levels(level)%p_xmin;  xmax = m_levels(level)%p_xmax; nx = m_levels(level)%p_nx; dx = m_levels(level)%p_dx
       ymin = m_levels(level)%p_ymin;  ymax = m_levels(level)%p_ymax; ny = m_levels(level)%p_ny; dy = m_levels(level)%p_dy
       zmin = m_levels(level)%p_zmin;  zmax = m_levels(level)%p_zmax; nz = m_levels(level)%p_nz; dz = m_levels(level)%p_dz
       dxdydz = dx*dy*dz
       
       if (aggregate_neutral_data) then       
         H_dens => m_levels(level)%H_dens
         H_ux   => m_levels(level)%H_ux
         H_uy   => m_levels(level)%H_uy
         H_uz   => m_levels(level)%H_uz
         H_temp => m_levels(level)%H_temp       
       else
         H_dens => m_levels(level)%g_H_dens
         H_ux   => m_levels(level)%g_H_ux
         H_uy   => m_levels(level)%g_H_uy
         H_uz   => m_levels(level)%g_H_uz
         H_temp => m_levels(level)%g_H_temp
       endif
       
       ! because trilinear interpolation in '+' direction
       xmin1 = xmin - dx
       ymin1 = ymin - dy
       zmin1 = zmin - dz
              
       
       ALLOCATE(tmp(-1:nx,-1:ny,-1:nz))
       
       H_dens = 0.0; H_temp = 0.0;  tmp = 0.0
       H_ux = 0.0; H_uy = 0.0; H_uz = 0.0
       
       do ithread = 0, nthreads - 1          
       
       particles =>my_particles(ithread)%particles__

       DO i = 1, my_particles(ithread)%my_npart

          xi = particles(1,i)
          yi = particles(2,i)
          zi = particles(3,i)

          IF ( (xi .LT. xmax) .AND. (xi .GT. xmin1) .AND. &
               (yi .LT. ymax) .AND. (yi .GT. ymin1) .AND. &
               (zi .LT. zmax) .AND. (zi .GT. zmin1) ) THEN

             mass = particles(0,i) * particles(8,i)
             vx = particles(4,i)
             vy = particles(5,i)
             vz = particles(6,i)
             
             ix = FLOOR((xi-xmin)/dx)
             iy = FLOOR((yi-ymin)/dy)
             iz = FLOOR((zi-zmin)/dz)
             
!             ix = min(max(0,ix),m_levels(ilev)%p_nx-1)
!             iy = min(max(0,iy),m_levels(ilev)%p_ny-1)
!             iz = min(max(0,iz),m_levels(ilev)%p_nz-1)
                                      
             x1 = xmin + ix*dx; x2 = x1 + dx
             y1 = ymin + iy*dy; y2 = y1 + dy
             z1 = zmin + iz*dz; z2 = z1 + dz
          
! use trilinear interpolation
             dx1 = xi - x1;  dx2 = x2 - xi
             dy1 = yi - y1;  dy2 = y2 - yi
             dz1 = zi - z1;  dz2 = z2 - zi

! distribute number density to adjacent grid points
             H_dens(ix  ,iy  ,iz)   = H_dens(ix  ,iy  ,iz)  + dx2*dy2*dz2/dxdydz* mass
             H_dens(ix+1,iy  ,iz)   = H_dens(ix+1,iy  ,iz)  + dx1*dy2*dz2/dxdydz* mass
             H_dens(ix  ,iy  ,iz+1) = H_dens(ix  ,iy  ,iz+1)+ dx2*dy2*dz1/dxdydz* mass
             H_dens(ix+1,iy  ,iz+1) = H_dens(ix+1,iy  ,iz+1)+ dx1*dy2*dz1/dxdydz* mass
             H_dens(ix  ,iy+1,iz)   = H_dens(ix  ,iy+1,iz)  + dx2*dy1*dz2/dxdydz* mass
             H_dens(ix+1,iy+1,iz)   = H_dens(ix+1,iy+1,iz)  + dx1*dy1*dz2/dxdydz* mass
             H_dens(ix  ,iy+1,iz+1) = H_dens(ix  ,iy+1,iz+1)+ dx2*dy1*dz1/dxdydz* mass
             H_dens(ix+1,iy+1,iz+1) = H_dens(ix+1,iy+1,iz+1)+ dx1*dy1*dz1/dxdydz* mass
          
! calculate bulk velocity
             H_ux(ix  ,iy  ,iz)   = H_ux(ix  ,iy  ,iz)  + dx2*dy2*dz2/dxdydz* vx* mass
             H_ux(ix+1,iy  ,iz)   = H_ux(ix+1,iy  ,iz)  + dx1*dy2*dz2/dxdydz* vx* mass
             H_ux(ix  ,iy  ,iz+1) = H_ux(ix  ,iy  ,iz+1)+ dx2*dy2*dz1/dxdydz* vx* mass
             H_ux(ix+1,iy  ,iz+1) = H_ux(ix+1,iy  ,iz+1)+ dx1*dy2*dz1/dxdydz* vx* mass
             H_ux(ix  ,iy+1,iz)   = H_ux(ix  ,iy+1,iz)  + dx2*dy1*dz2/dxdydz* vx* mass
             H_ux(ix+1,iy+1,iz)   = H_ux(ix+1,iy+1,iz)  + dx1*dy1*dz2/dxdydz* vx* mass
             H_ux(ix  ,iy+1,iz+1) = H_ux(ix  ,iy+1,iz+1)+ dx2*dy1*dz1/dxdydz* vx* mass
             H_ux(ix+1,iy+1,iz+1) = H_ux(ix+1,iy+1,iz+1)+ dx1*dy1*dz1/dxdydz* vx* mass

             H_uy(ix  ,iy  ,iz)   = H_uy(ix  ,iy  ,iz)  + dx2*dy2*dz2/dxdydz* vy* mass
             H_uy(ix+1,iy  ,iz)   = H_uy(ix+1,iy  ,iz)  + dx1*dy2*dz2/dxdydz* vy* mass
             H_uy(ix  ,iy  ,iz+1) = H_uy(ix  ,iy  ,iz+1)+ dx2*dy2*dz1/dxdydz* vy* mass
             H_uy(ix+1,iy  ,iz+1) = H_uy(ix+1,iy  ,iz+1)+ dx1*dy2*dz1/dxdydz* vy* mass
             H_uy(ix  ,iy+1,iz)   = H_uy(ix  ,iy+1,iz)  + dx2*dy1*dz2/dxdydz* vy* mass
             H_uy(ix+1,iy+1,iz)   = H_uy(ix+1,iy+1,iz)  + dx1*dy1*dz2/dxdydz* vy* mass
             H_uy(ix  ,iy+1,iz+1) = H_uy(ix  ,iy+1,iz+1)+ dx2*dy1*dz1/dxdydz* vy* mass
             H_uy(ix+1,iy+1,iz+1) = H_uy(ix+1,iy+1,iz+1)+ dx1*dy1*dz1/dxdydz* vy* mass

             H_uz(ix  ,iy  ,iz)   = H_uz(ix  ,iy  ,iz)  + dx2*dy2*dz2/dxdydz* vz* mass
             H_uz(ix+1,iy  ,iz)   = H_uz(ix+1,iy  ,iz)  + dx1*dy2*dz2/dxdydz* vz* mass
             H_uz(ix  ,iy  ,iz+1) = H_uz(ix  ,iy  ,iz+1)+ dx2*dy2*dz1/dxdydz* vz* mass
             H_uz(ix+1,iy  ,iz+1) = H_uz(ix+1,iy  ,iz+1)+ dx1*dy2*dz1/dxdydz* vz* mass
             H_uz(ix  ,iy+1,iz)   = H_uz(ix  ,iy+1,iz)  + dx2*dy1*dz2/dxdydz* vz* mass
             H_uz(ix+1,iy+1,iz)   = H_uz(ix+1,iy+1,iz)  + dx1*dy1*dz2/dxdydz* vz* mass
             H_uz(ix  ,iy+1,iz+1) = H_uz(ix  ,iy+1,iz+1)+ dx2*dy1*dz1/dxdydz* vz* mass
             H_uz(ix+1,iy+1,iz+1) = H_uz(ix+1,iy+1,iz+1)+ dx1*dy1*dz1/dxdydz* vz* mass
 
          END IF

       END DO   ! for all particles
       
       end do   ! for all threads
       
       !write (str_buf,*) '        output_grid, density, step A '
       !call pout(str_buf)   

! need to use AllReduce so that each proc has the full normalized velocity
! and we can then correctly compute temperature below
       CALL MPI_AllReduce(H_dens, tmp, (nx+2)*(ny+2)*(nz+2),&
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       H_dens = tmp
       CALL MPI_AllReduce(H_ux, tmp, (nx+2)*(ny+2)*(nz+2),&
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       H_ux = tmp
       CALL MPI_AllReduce(H_uy, tmp, (nx+2)*(ny+2)*(nz+2),&
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       H_uy = tmp
       CALL MPI_AllReduce(H_uz, tmp, (nx+2)*(ny+2)*(nz+2),&
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       H_uz = tmp
       
       
       !write (str_buf,*) '        output_grid, density, step B '
       !call pout(str_buf)   
       
       DO kk = -1, nx
          WHERE (H_dens(kk,:,:) /= 0.0)
             H_ux(kk,:,:)   = H_ux(kk,:,:)   / H_dens(kk,:,:)
             H_uy(kk,:,:)   = H_uy(kk,:,:)   / H_dens(kk,:,:)
             H_uz(kk,:,:)   = H_uz(kk,:,:)   / H_dens(kk,:,:)
          END WHERE
       END DO

     !~ !$OMP PARALLEL default(none) shared (H_ux,H_uy,H_uz,H_dens)
     !~ !$OMP sections
     !~ 
     !~ !$OMP section
     !~ WHERE (H_dens /= 0.0)
     !~    H_ux   = H_ux   / H_dens
     !~ END WHERE
     !~ 
     !~ !$OMP section
     !~ WHERE (H_dens /= 0.0)
     !~    H_uy   = H_uy   / H_dens            
     !~ END WHERE
     !~ 
     !~ !$OMP section
     !~ WHERE (H_dens /= 0.0)
     !~    H_uz   = H_uz   / H_dens
     !~ END WHERE      
     !~ 
     !~!$OMP END sections
     !~!$OMP END PARALLEL
       

      !write (str_buf,*) '        output_grid, density, step C '
      !call pout(str_buf)   
       

       do ithread = 0, nthreads - 1     

       particles =>my_particles(ithread)%particles__     

       DO i = 1, my_particles(ithread)%my_npart

          xi = particles(1,i)
          yi = particles(2,i)
          zi = particles(3,i)

          IF ( (xi .LT. xmax) .AND. (xi .GT. xmin) .AND. &
               (yi .LT. ymax) .AND. (yi .GT. ymin) .AND. &
               (zi .LT. zmax) .AND. (zi .GT. zmin) ) THEN

             mass = particles(0,i) * particles(8,i)
             vx = particles(4,i)
             vy = particles(5,i)
             vz = particles(6,i)

             ix = FLOOR((xi-xmin)/dx)
             iy = FLOOR((yi-ymin)/dy)
             iz = FLOOR((zi-zmin)/dz)

             x1 = xmin + ix*dx; x2 = x1 + dx
             y1 = ymin + iy*dy; y2 = y1 + dy
             z1 = zmin + iz*dz; z2 = z1 + dz

! use trilinear interpolation
             dx1 = xi - x1;  dx2 = x2 - xi
             dy1 = yi - y1;  dy2 = y2 - yi
             dz1 = zi - z1;  dz2 = z2 - zi

! distribute number density to adjacent grid points
             vrel_sq = (vx - H_ux(ix  ,iy  ,iz  ))**2 &
                     + (vy - H_uy(ix  ,iy  ,iz  ))**2 &
                     + (vz - H_uz(ix  ,iy  ,iz  ))**2
             H_temp(ix  ,iy  ,iz)   = H_temp(ix  ,iy  ,iz)  + dx2*dy2*dz2/dxdydz*vrel_sq * mass
             vrel_sq = (vx - H_ux(ix+1,iy  ,iz  ))**2 &
                     + (vy - H_uy(ix+1,iy  ,iz  ))**2 &
                     + (vz - H_uz(ix+1,iy  ,iz  ))**2
             H_temp(ix+1,iy  ,iz)   = H_temp(ix+1,iy  ,iz)  + dx1*dy2*dz2/dxdydz*vrel_sq * mass
             vrel_sq = (vx - H_ux(ix  ,iy  ,iz+1))**2 &
                     + (vy - H_uy(ix  ,iy  ,iz+1))**2 &
                     + (vz - H_uz(ix  ,iy  ,iz+1))**2
             H_temp(ix  ,iy  ,iz+1) = H_temp(ix  ,iy  ,iz+1)+ dx2*dy2*dz1/dxdydz*vrel_sq * mass
             vrel_sq = (vx - H_ux(ix+1,iy  ,iz+1))**2 &
                     + (vy - H_uy(ix+1,iy  ,iz+1))**2 &
                     + (vz - H_uz(ix+1,iy  ,iz+1))**2
             H_temp(ix+1,iy  ,iz+1) = H_temp(ix+1,iy  ,iz+1)+ dx1*dy2*dz1/dxdydz*vrel_sq * mass
             vrel_sq = (vx - H_ux(ix  ,iy+1,iz  ))**2 &
                     + (vy - H_uy(ix  ,iy+1,iz  ))**2 &
                     + (vz - H_uz(ix  ,iy+1,iz  ))**2
             H_temp(ix  ,iy+1,iz)   = H_temp(ix  ,iy+1,iz)  + dx2*dy1*dz2/dxdydz*vrel_sq * mass
             vrel_sq = (vx - H_ux(ix+1,iy+1,iz  ))**2 &
                     + (vy - H_uy(ix+1,iy+1,iz  ))**2 &
                     + (vz - H_uz(ix+1,iy+1,iz  ))**2
             H_temp(ix+1,iy+1,iz)   = H_temp(ix+1,iy+1,iz)  + dx1*dy1*dz2/dxdydz*vrel_sq * mass
             vrel_sq = (vx - H_ux(ix  ,iy+1,iz+1))**2 &
                     + (vy - H_uy(ix  ,iy+1,iz+1))**2 &
                     + (vz - H_uz(ix  ,iy+1,iz+1))**2
             H_temp(ix  ,iy+1,iz+1) = H_temp(ix  ,iy+1,iz+1)+ dx2*dy1*dz1/dxdydz*vrel_sq * mass
             vrel_sq = (vx - H_ux(ix+1,iy+1,iz+1))**2 &
                     + (vy - H_uy(ix+1,iy+1,iz+1))**2 &
                     + (vz - H_uz(ix+1,iy+1,iz+1))**2
             H_temp(ix+1,iy+1,iz+1) = H_temp(ix+1,iy+1,iz+1)+ dx1*dy1*dz1/dxdydz*vrel_sq * mass
          
          END IF
       END DO
       
       end do
       
       !write (str_buf,*) '        output_grid, temperature, step A '
       !call pout(str_buf)   

       CALL MPI_AllReduce(H_temp, tmp, (nx+2)*(ny+2)*(nz+2),&
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       H_temp = tmp

       
      DO kk = -1, nx
         WHERE (H_dens(kk,:,:) /= 0.0)
            H_temp(kk,:,:) = H_temp(kk,:,:) / H_dens(kk,:,:)
         END WHERE
      END DO

! normalise grid quantities
      H_dens = H_dens * mega_number / dxdydz  ! per m^3
      H_temp = H_temp * neutral_mass / (3.0*k_B)
      
      !write (str_buf,*) '        output_grid, temperature, step B '
      !call pout(str_buf)   

!~      !$omp parallel  default(none) shared(H_dens,H_temp,neutral_mass,k_b,mega_number,dxdydz,nx) private(kk)
!~      !$omp sections
!~      
!~      !$omp section
!~      DO kk = -1, nx
!~         WHERE (H_dens(kk,:,:) /= 0.0)
!~            H_temp(kk,:,:) = H_temp(kk,:,:) / H_dens(kk,:,:)
!~         END WHERE
!~      END DO      
!~
!~      !$omp section
!~       H_dens = H_dens * mega_number / dxdydz  ! per m^3      
!~      !$omp end  sections
!~      
!~      !$omp end parallel
!~      
!~      H_temp = H_temp * neutral_mass / (3.0*k_B)
         

      DEALLOCATE(tmp)
      
      
      if (aggregate_neutral_data) then       
        m_levels(level)%g_H_dens = m_levels(level)%g_H_dens + H_dens
        m_levels(level)%g_H_ux   = m_levels(level)%g_H_ux   + H_ux
        m_levels(level)%g_H_uy   = m_levels(level)%g_H_uy   + H_uy
        m_levels(level)%g_H_uz   = m_levels(level)%g_H_uz   + H_uz
        m_levels(level)%g_H_temp = m_levels(level)%g_H_temp + H_temp   

        m_levels(level)%num_samples = m_levels(level)%num_samples + 1
      endif

    END DO  ! loop over levels
    

  END SUBROUTINE collect_neutral_distr


  function secs(sys_time)
  
  implicit none
  INTEGER,intent(in) :: sys_time(8)
  real(8) :: secs
  
  secs = sys_time(5)*3600+sys_time(6)*60+sys_time(7)+1d-3*sys_time(8)

  end function secs
!
!***********************************************************************
!***********************************************************************
!
 SUBROUTINE mc_neutrals3d(run_time_in)        
  USE global
  USE level1_subroutines
  
  IMPLICIT NONE
  
  REAL(KIND=dp), intent(in) :: run_time_in  
  
  REAL(KIND=dp) :: max_dens, max_temp  
  INTEGER :: loop, ilev, next_balance
  REAL(KIND=dp) p_mc_time,p_timestep ! private copies of coresponding variables
  REAL(KIND=dp) :: run_time
  
  INTEGER :: start1(8), finish1(8),start2(8), finish2(8), start_main(8)
  INTEGER :: start3(8), finish3(8),start4(8), finish4(8)
  real(8) :: start1_sec, finish1_sec,secs2
  
  real(8), dimension(0:nthreads-1) :: time_step_sec
  
  CHARACTER (LEN = 200) str_buf
  
  run_time = run_time_in*secs_per_year  
  
  do ilev = 0, p_nlevels - 1        
    m_levels(ilev)%source_m(:,:,:)   = 0.0_dp
    m_levels(ilev)%source_px(:,:,:)  = 0.0_dp
    m_levels(ilev)%source_py(:,:,:)  = 0.0_dp
    m_levels(ilev)%source_pz(:,:,:)  = 0.0_dp
    m_levels(ilev)%source_mvsq(:,:,:) = 0.0_dp
    m_levels(ilev)%source_nchex(:,:,:) = 0.0_dp
    m_levels(ilev)%num_samples = 0
    if (aggregate_neutral_data) then
      m_levels(ilev)%g_H_dens = 0.0_dp
      m_levels(ilev)%g_H_ux   = 0.0_dp
      m_levels(ilev)%g_H_uy   = 0.0_dp
      m_levels(ilev)%g_H_uz   = 0.0_dp
      m_levels(ilev)%g_H_temp = 0.0_dp
    endif    
  enddo
  
  !~!$omp parallel
  !call initial_particles_distribution
  !~!$omp end parallel  
  !~call output_grid
  !return
  
  !print *, 'begin parallel loop '
  
!  IF (my_rank .eq. 0)  then        
!    write (str_buf,*) 'begin parallel loop '
!    call pout(str_buf)   
!  endif
  
  call date_and_time(VALUES=start_main)
  
  next_balance = 0
  
  !$omp parallel private(p_mc_time,p_timestep,loop,start1_sec,finish1_sec, &
  !$omp&  start1, finish1, start3, finish3, start4, finish4, secs2) 
  p_mc_time  = 0.0_dp
  p_timestep = timestep
  loop = 0
  DO WHILE (p_mc_time .LT. run_time)
      
    call date_and_time(VALUES=start1)
             
    CALL move
    CALL split
          
    call date_and_time(VALUES=start3)
    secs2 = secs(start3)
    start1_sec = secs(start1)    
     
    time_step_sec(my_thread) = secs2-start1_sec
    
     
    if (loop .eq. next_balance) then
       !$omp barrier
       !$omp master
       
       last_time_step_sec = sum(time_step_sec)/nthreads
                                   
       CALL balance_npart(loop,next_balance)       
       
       !$omp end master      
       !$omp barrier
    endif   
    
    !elseif (mod(loop,6).eq.0) then
    !   !$omp barrier
    !   !$omp master
    !                 
    !   call balance_node()
    !   
    !   !$omp end master      
    !   !$omp barrier
    !   
    !endif
    
    if ((aggregate_neutral_data).and.(mod(loop,m_levels(p_nlevels-1)%sampling_dt).eq.0)) then
      !$omp barrier
      !$omp master
        
      do ilev = p_nlevels - 1, 0, -1
      if (mod(loop,m_levels(ilev)%sampling_dt).eq.0) then
        call collect_neutral_distr(ilev)
      endif
      enddo
      
      !$omp end master      
      !$omp barrier    
    endif
    
    
    

    p_mc_time = p_mc_time + p_timestep
    IF (run_time - p_mc_time .LT. p_timestep) THEN 
      p_timestep = run_time - p_mc_time  ! adjust timestep so we finish on time
    END IF

    loop = loop + 1
    
    if (my_thread.eq.0) then
      call date_and_time(VALUES=finish1)
      finish1_sec = secs(finish1)    
        
!      write (str_buf,'(a,i6,a,i10,a,f8.3,a,f8.3,a,f8.3,a)') ' loop ',loop, ', my particles ',my_particles(0)%my_npart, &
!      ' (',finish1_sec-start1_sec,' sec) [move/split ',time_step_sec(0),', balance ',finish1_sec-secs2,']'

!      call pout(str_buf)                               
    endif
     

  END DO
  
  call collect_sources_node
  !$omp end parallel
  
  call date_and_time(VALUES=finish2)
  
!  write (str_buf,'(a,f8.3,a)') 'end parallel loop (',secs(finish2)-secs(start_main),' sec)'
!  call pout(str_buf)   
  
  
  
  call date_and_time(VALUES=start1)
  CALL output_source(run_time)     
  call date_and_time(VALUES=start2)            
  CALL output_grid  
  call date_and_time(VALUES=start3)    
  
!  write (str_buf,'(a,f8.3,a)') 'output_source (',secs(start2)-secs(start1),' sec)'
!  call pout(str_buf)   
  
!  write (str_buf,'(a,f8.3,a)') 'output_grid (',secs(start3)-secs(start2),' sec)'
!  call pout(str_buf)   
  
!  write (str_buf,'(a,f8.3,a)') 'end mc_neutrals3d (total ',secs(start3)-secs(start_main),' sec)'
!  call pout(str_buf)   
  

 END SUBROUTINE mc_neutrals3d

END MODULE level2_subroutines
  
