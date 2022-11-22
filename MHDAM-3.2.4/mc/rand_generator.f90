
! Random generator used in fortran. Original file is taken from www.fortran.com
! Rewritten by Sergey Borovikov to be thread friendly 

module random_gen_mod
implicit none
private

public :: rand_seed,rand_number
INTEGER(4),private :: x=123456789, y=362436069, z=521288629, w=916191069

!$OMP threadprivate (x,y,z,w)

contains

   FUNCTION kiss ()
      integer(4) :: kiss
! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123;  Default seeds x,y,z,w.
!  Set your own seeds with statement i=kisset(ix,iy,iz,iw).
!
      x = 69069 * x + 1327217885
      y = m (m (m (y, 13), - 17), 5)
      z = 18000 * iand (z, 65535) + ishft (z, - 16)
      w = 30903 * iand (w, 65535) + ishft (w, - 16)
      kiss = x + y + ishft (z, 16) + w
   contains
      function m(k, n)
         integer(4) :: m, k, n
         m = ieor (k, ishft (k, n) )
      end function m
   END FUNCTION kiss
   
   subroutine rand_seed(seeds)
      integer(4), dimension(1:4) :: seeds
      x = seeds(1)
      y = seeds(2)
      z = seeds(3)
      w = seeds(4)
   end subroutine rand_seed
   
   subroutine rand_number_orig(n)
      integer(8) :: v2_32 = 4294967296_8
      integer(8) :: v2_31 = 2147483648_8
      real(8) :: n,tmp
      integer(4) :: kiss_value
      integer(8) :: kiss_value8
                  
      kiss_value  = kiss()
      kiss_value8 = kiss_value 
      kiss_value8 = kiss_value8 + v2_31
      
      tmp = real(kiss_value8,8)/real(v2_32,8)
      
      n = tmp
      
      !PRINT *,kiss_value,kiss_value8,n
   end subroutine rand_number_orig
   
   subroutine rand_number(n)      
      integer(8), parameter :: v2_31 = 2147483648_8
      real(8), intent(out) :: n
      real(8) :: tmp
      integer(4) :: kiss_value,shift1,shift2
      integer(8) :: kiss_value8
            
      
            
      x = 69069 * x + 1327217885
      shift1 = ieor (y,      ishft (y,       13) )
      shift2 = ieor (shift1, ishft (shift1, -17) )
      y      = ieor (shift2, ishft (shift2,   5) )                  
      z = 18000 * iand (z, 65535) + ishft (z, - 16)
      w = 30903 * iand (w, 65535) + ishft (w, - 16)
      kiss_value = x + y + ishft (z, 16) + w
                                            
      kiss_value8 = kiss_value 
      kiss_value8 = kiss_value8 + v2_31
      
      tmp = kiss_value8*2.3283064365386962890625D-10
      
      n = tmp
      
      !PRINT *,kiss_value,kiss_value8,n
   end subroutine rand_number
      

   
end module random_gen_mod

