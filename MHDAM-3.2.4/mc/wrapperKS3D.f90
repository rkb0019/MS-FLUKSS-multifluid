
  subroutine mc3d_init_call_1(au_nlevels)
 
    use global
    USE level2_subroutines
    implicit none
        
    INTEGER, INTENT(in) :: au_nlevels    
    
    call mc_init_call_1(au_nlevels)         
    return 
  end
  
  subroutine mc3d_init_call_2(au_LISM_nH, au_LISM_vH, au_LISM_TH,  & 
                  au_total_neutrals, au_restart, au_load_file)
 
    use global
    USE level2_subroutines
    implicit none
        
    REAL(KIND=dp), INTENT(in) :: au_LISM_nH, au_LISM_TH
    REAL(KIND=dp), DIMENSION(0:2), INTENT(in) :: au_LISM_vH    
    REAL(KIND=dp), INTENT(in) :: au_total_neutrals
    INTEGER, INTENT(in)       :: au_restart
    CHARACTER (LEN = lenstr), INTENT(in)  :: au_load_file
    
    
    call mc_init_call_2(au_LISM_nH, au_LISM_vH, au_LISM_TH, & 
                  au_total_neutrals, au_restart, au_load_file)         
    return 
  end
  
  subroutine setup_level_3d(au_level,au_ref_ratio,au_Lo,au_Hi,au_dx,au_size, &
                         au_source_m, au_source_px, au_source_py, au_source_pz, au_source_mvsq, &
                         au_plasma, &                         
                         au_nchex, au_H_dens, au_H_ux, au_H_uy, au_H_uz, au_H_temp)    
  
    use global
    USE level2_subroutines
    
    implicit none
  
    INTEGER, INTENT(in)                 :: au_level, au_ref_ratio
    REAL(KIND=dp), DIMENSION(0:2), INTENT(in) :: au_Lo
    REAL(KIND=dp), DIMENSION(0:2), INTENT(in) :: au_Hi
    REAL(KIND=dp), INTENT(in)           :: au_dx
    INTEGER, DIMENSION(0:2), INTENT(in) :: au_size
    
    REAL(KIND=dp), DIMENSION (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1), TARGET, INTENT(in) :: au_source_m 
    REAL(KIND=dp), DIMENSION (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1), TARGET, INTENT(in) :: au_source_px 
    REAL(KIND=dp), DIMENSION (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1), TARGET, INTENT(in) :: au_source_py 
    REAL(KIND=dp), DIMENSION (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1), TARGET, INTENT(in) :: au_source_pz 
    REAL(KIND=dp), DIMENSION (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1), TARGET, INTENT(in) :: au_source_mvsq    
    
    REAL(KIND=dp), DIMENSION (0:au_size(0)-1,0:au_size(1)-1,0:au_size(2)-1,0:WCOMP-1), TARGET, INTENT(in) :: au_plasma    
    
    REAL(KIND=dp), DIMENSION (-1:au_size(0),-1:au_size(1),-1:au_size(2)), TARGET, INTENT(in) :: au_nchex, au_H_dens, &
                              au_H_ux, au_H_uy, au_H_uz, au_H_temp 
    
    call setup_level(au_level,au_ref_ratio,au_Lo,au_Hi,au_dx,au_size, &
                         au_source_m, au_source_px, au_source_py, au_source_pz, au_source_mvsq, &
                         au_plasma, &
                         au_nchex, au_H_dens, au_H_ux, au_H_uy, au_H_uz, au_H_temp)    
    
    return 
  end

  subroutine RUN_MC_NEUTRALS_3D(run_time)
 
    use global
    USE level2_subroutines
    
    implicit none
         
    REAL(KIND=dp), INTENT(in) :: run_time
          
    CALL mc_neutrals3d(run_time)
  
    return 
  end
    
    
  subroutine RUN_MC_OUTPUT_RAW_3D(chkfile)
 
    USE level2_subroutines
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
  
  subroutine finalize_kinetic3d_f()
    
    USE level2_subroutines 
    call finalize
                
    return 
  end
