!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Esta subrutina se encarga de poner a todas los segmentos dentro del slab

subroutine pxssv(il)

use system
use MPI
use const
use transform
use solventchains
implicit none
    
integer j, ii, jj,i, il
real*8 xx(3)
real*8 x(3)
real*8 v(3)
integer testsystem
integer testsystemr
integer testsystemc
integer testsystem_cylinder
integer testsystem_superellipse
integer testsystem_cube
integer testsystem_cuboctahedron
integer flag
integer aa
real*8 com(3)

integer, external :: PBCREFI, PBCSYMI

com = 0.0
do j = 1, longsv
 com(:) = com(:) + in1sv(j,:)
enddo

com = com/float(longsv)

ntranssv(il) = ingsv

do j = 1, longsv

  x(:) = in1sv(j,:) - com(:) ! substract center of coordinates
  v = MATMUL(MAT,x) ! transformed space

  aa = nint(v(1)/delta)
  pxsv(il,j) = aa

  aa = nint(v(2)/delta)
  pysv(il,j) = aa

  aa = nint(v(3)/delta)
  pzsv(il,j) = aa

enddo ! j
return
end
      
