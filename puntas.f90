subroutine puntas(counter)

use cube
use system
use mmask
use kaist
implicit none
integer i,j,k
integer nx, ny, nz
real*8 vertices(1,3) !Número de centros de las gaussianas, 3
real*8 Xarray(3)
real*8, dimension(:,:), allocatable :: vert_test
real*8 points_test(3)
real*8 sigma
integer counter
character*5 title
real*8, parameter :: pi  = 4 * atan(1.0_8) 

! Vértices del cubo

vertices(1,:) = (/c_cube(1) + l_cube/2.0,c_cube(2) + l_cube/2.0,c_cube(3) + l_cube/2.0/) !Solo este vertice para el octavo del cubo
!vertices(2,:) = (/c_cube(1) - l_cube/2.0,c_cube(2) + l_cube/2.0,c_cube(3) + l_cube/2.0/)
!vertices(3,:) = (/c_cube(1) + l_cube/2.0,c_cube(2) - l_cube/2.0,c_cube(3) + l_cube/2.0/)
!vertices(4,:) = (/c_cube(1) + l_cube/2.0,c_cube(2) + l_cube/2.0,c_cube(3) - l_cube/2.0/)
!vertices(5,:) = (/c_cube(1) - l_cube/2.0,c_cube(2) - l_cube/2.0,c_cube(3) + l_cube/2.0/)
!vertices(6,:) = (/c_cube(1) + l_cube/2.0,c_cube(2) - l_cube/2.0,c_cube(3) - l_cube/2.0/)
!vertices(7,:) = (/c_cube(1) - l_cube/2.0,c_cube(2) + l_cube/2.0,c_cube(3) - l_cube/2.0/)
!vertices(8,:) = (/c_cube(1) - l_cube/2.0,c_cube(2) - l_cube/2.0,c_cube(3) - l_cube/2.0/)

!Caras del cubo (octavo)

!vertices(1,:) = (/c_cube(1) - l_cube/2.0,c_cube(2) - l_cube/2.0,c_cube(3) + l_cube/2.0/)
!vertices(2,:) = (/c_cube(1) + l_cube/2.0,c_cube(2) - l_cube/2.0,c_cube(3) - l_cube/2.0/)
!vertices(3,:) = (/c_cube(1) - l_cube/2.0,c_cube(2) + l_cube/2.0,c_cube(3) - l_cube/2.0/)

!Aristas del cubo (octavo)

!vertices(1,:) = (/c_cube(1) - l_cube/2.0,c_cube(2) + l_cube/2.0,c_cube(3) + l_cube/2.0/)
!vertices(2,:) = (/c_cube(1) + l_cube/2.0,c_cube(2) - l_cube/2.0,c_cube(3) + l_cube/2.0/)
!vertices(3,:) = (/c_cube(1) + l_cube/2.0,c_cube(2) + l_cube/2.0,c_cube(3) - l_cube/2.0/)

! testeo con un punto arbitrario

allocate(vert_test(1,3)) !Número de centros de las gaussianas, 3
vert_test = vertices
sigma = kp

! Matriz con las posiciones de cada celdilla del array y la gaussiana

do i = 1,dimx
 do j = 1,dimy
  do k = 1,dimz
   Xarray = (/delta/2.0 + delta*(i-1),delta/2.0 + delta*(j-1),delta/2.0 + delta*(k-1)/)
   mask(i,j,k) = gauss_interact(Xarray,vert_test,sigma)
  enddo
 enddo
enddo

title = 'punta'
call savetodisk(mask, title, counter)

contains

double precision function  gauss_interact(X,points,sigma)
real*8 X(3)
real*8, dimension(:,:), allocatable :: points
integer ii
real*8 sigma, xx, yy, zz

gauss_interact = 0.0


do ii = 1,size(points,1)


 xx = (X(1) - points(ii,1))
 yy = (X(2) - points(ii,2))
 zz = (X(3) - points(ii,3))
 gauss_interact = gauss_interact + dexp(-(1.0/(2.0*sigma**2))*(xx**2 + yy**2 + zz**2))
end do

end function

end subroutine

