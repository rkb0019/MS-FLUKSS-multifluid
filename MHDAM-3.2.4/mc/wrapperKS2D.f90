
  
  subroutine RUN_MC_INIT_2D(p_xmax_in, p_zmin_in, p_zmax_in, &
                     ngrids_in, grids_nx_in, grids_nz_in, &
                     grids_boundaries_x_in, grids_boundaries_z_in, &
                     restart_in, load_file_in, total_neutrals, &
                     lism_nh_in, lism_vh_in, lism_th_in, verbosity_in, &
                     photoionize_in)
 
    USE neutrals_routines2_2D
    implicit none
       
    REAL*8 :: p_xmax_in, p_zmin_in, p_zmax_in, total_neutrals
    INTEGER :: ngrids_in   
    INTEGER, DIMENSION(1:ngrids_in) :: grids_nx_in, grids_nz_in
    REAL*8,  DIMENSION(0:ngrids_in) :: grids_boundaries_x_in, grids_boundaries_z_in
    byte ::  load_file_in(1:60)
    integer :: restart_in, verbosity_in, photoionize_in
    REAL*8 :: lism_nh_in, lism_vh_in, lism_th_in
   
    integer :: i
    CHARACTER (LEN = 60) :: loadfilestr
    REAL*8 :: lism_nh   , lism_vh   , lism_th

    loadfilestr = ' ' 

    if (restart_in == 1) then
      i = 1
      do while (load_file_in(i) .ne. 0)
        loadfilestr(i:i) = char(load_file_in(i),1)
        i = i + 1
     enddo 
    endif

    lism_nh = lism_nh_in*1D+6;   ! 1/cm^3 -> 1/m^3
    lism_vh = lism_vh_in*1D-2;   ! cm/s   ->  m/s
    lism_th = lism_th_in


    CALL mc_init(p_xmax_in, p_zmin_in, p_zmax_in, &
                     ngrids_in, grids_nx_in, grids_nz_in, &
                     grids_boundaries_x_in, grids_boundaries_z_in, &
                     restart_in, loadfilestr, total_neutrals, &
                     lism_nh, lism_vh, lism_th, verbosity_in, &
                     photoionize_in)

    return 
  end
  
  subroutine RUN_PASS_PLASMA_GRID(kk, nx, ny, grid, region)
  
    USE neutrals_routines2_2D
    implicit none
  
    REAL*8,  DIMENSION(-nx:nx-1,0:ny-1,1:4), INTENT(in) :: grid
    INTEGER, DIMENSION(    -nx:nx-1,0:ny-1), INTENT(in) :: region    
    INTEGER :: kk, nx, ny
    
    call pass_plasma_grid_mhdam(kk, ny, nx, grid, region)
    
    return 
  end

  subroutine RUN_MC_NEUTRALS_2D(run_time, min_nchex)
 
    USE neutrals_routines2_2D
    implicit none
         
    REAL*8 :: run_time, min_nchex
          
    CALL mc_neutrals(run_time, min_nchex)
  
    return 
  end
  
  subroutine RUN_RETURN_NEUTRAL_GRID(kk, nx, ny, region, grid)
    USE neutrals_routines2_2D
    implicit none
    
    REAL*8, DIMENSION(-nx:nx-1,0:ny-1,1:4), INTENT(out) :: grid
    INTEGER :: kk, nx, ny, region
    
    call return_neutral_grid_mhdam(kk, ny, nx, region, grid)
    return 
  end
  
  SUBROUTINE RUN_RETURN_SOURCE_GRID(kk, nx, ny, grid)
    USE neutrals_routines2_2D
    implicit none
    
    REAL*8, DIMENSION(-nx:nx-1,0:ny-1,1:5), INTENT(out) :: grid
    INTEGER :: kk, nx, ny
  
    call return_source_grid_mhdam(kk, ny, nx, grid)  
    return 
  end
  
    
  subroutine RUN_MC_OUTPUT_RAW_2D(chkfile)
 
    USE neutrals_routines2_2D
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
