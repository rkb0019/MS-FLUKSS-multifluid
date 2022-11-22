
  
  subroutine RUN_MC_INIT_1D(p_nz_in, p_zmin_in, p_zmax_in, restart_in, loadfile, total_neutrals, &
   lism_nh_in, lism_vh_in, lism_th_in, verbosity_in)
 
   USE neutrals_routines2_1D
    implicit none
        
   REAL*8 :: p_zmin_in, p_zmax_in
   INTEGER :: p_nz_in, total_neutrals
   byte ::  loadfile(1:60)
   integer :: restart_in, verbosity_in
   CHARACTER (LEN = 60) :: loadfilestr
   integer :: i

   REAL*8 :: lism_nh_in, lism_vh_in, lism_th_in
   REAL*8 :: lism_nh   , lism_vh   , lism_th

   loadfilestr = ' ' 

   if (restart_in == 1) then
       i = 1
	   do while (loadfile(i) .ne. 0)
		loadfilestr(i:i) = char(loadfile(i),1)		
		i = i + 1
	   enddo 

!	   do while (i<=60)
!		loadfilestr(i:i) = char(32,1)		
!		i = i + 1
!	   enddo 

!	   loadfilestr(i:) = '.'
!	   PRINT *, 'loadfile: ',loadfilestr
   endif

   lism_nh = lism_nh_in*1D+6;   ! 1/cm^3 -> 1/m^3
   lism_vh = lism_vh_in*1D-2;   ! cm/s   ->  m/s
   lism_th = lism_th_in

  
   CALL mc_init(p_nz_in, p_zmin_in, p_zmax_in, restart_in, loadfilestr, total_neutrals, &
    lism_nh, lism_vh, lism_th, verbosity_in)
  
  return 
  end

  subroutine RUN_MC_NEUTRALS_1D(run_time, min_nchex, & 
    p_dens_in, p_ux_in, p_uy_in, p_uz_in, p_temp_in, p_region_in, &
	src_momx, src_momy, src_momz, src_energy, &
	n_elements, &
	grid_dens1, grid_ux1, grid_uy1, grid_uz1, grid_temp1 )
 
  USE neutrals_routines2_1D

  implicit none
  
   integer n_elements
   REAL*8 :: p_dens_in (1:n_elements), p_temp_in(1:n_elements)
   REAL*8 :: p_ux_in(1:n_elements), p_uy_in(1:n_elements), p_uz_in(1:n_elements), p_region_in(1:n_elements)   
   REAL*8 :: src_momx(1:n_elements), src_momy(1:n_elements), src_momz(1:n_elements), src_energy(1:n_elements)
   REAL*8 :: grid_dens1(1:n_elements), grid_temp1(1:n_elements)
   REAL*8 :: grid_ux1(1:n_elements), grid_uy1(1:n_elements),  grid_uz1(1:n_elements)
   REAL*8 :: run_time, min_nchex
        
  
   CALL mc_neutrals(run_time, min_nchex, p_dens_in, p_ux_in, &
       p_uy_in, p_uz_in, p_temp_in, p_region_in, &
       src_momx, src_momy, src_momz, src_energy, &
       grid_dens1, grid_ux1, grid_uy1, grid_uz1, grid_temp1)
  
  return 
  end
    
  subroutine RUN_MC_OUTPUT_RAW_1D(chkfile)
 
   USE neutrals_routines2_1D
    implicit none
        
   
   byte ::  chkfile(1:60)   
   CHARACTER (LEN = 60) :: chkfilestr
   integer :: i

   chkfilestr = ' ' 

   
   i = 1
   do while (chkfile(i) .ne. 0)
	chkfilestr(i:i) = char(chkfile(i),1)		
	i = i + 1
   enddo 


  
   CALL output_raw(chkfilestr)
  
  return 
  end
