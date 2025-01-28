integer function testsystem_cuboctahedron(x)
use system
use transform
use ellipsoid
use COrotMod

implicit none
real*8 x(3), xx(3), v(3), maxx(3)
integer j, i
real*8 vect
real*8 mmmult
real, external :: PBCSYMR, PBCREFR
real*8 dims(3)
logical test1, test2, test3, test4
logical test5, test6, test7

dims(1) = delta*dimx
dims(2) = delta*dimy
dims(3) = delta*dimz

maxx(1) = float(dimx)*delta
maxx(2) = float(dimy)*delta
maxx(3) = float(dimz)*delta

! collision with walls and out of system

testsystem_cuboctahedron = 0

v = MATMUL(MAT,x) ! to transformed space

! out-of-boundaries check is performed in transformed space

if (v(1).le.0.0) then
 if (PBC(1).eq.2)testsystem_cuboctahedron = -1
 if (PBC(1).eq.0)testsystem_cuboctahedron = -2
endif

if (v(1).gt.(float(dimx)*delta)) then
 if (PBC(2).eq.2)testsystem_cuboctahedron = -1
 if (PBC(2).eq.0)testsystem_cuboctahedron = -2
endif

if (v(2).le.0.0) then
 if (PBC(3).eq.2)testsystem_cuboctahedron = -1
 if (PBC(3).eq.0)testsystem_cuboctahedron = -2
endif

if (v(2).gt.(float(dimy)*delta)) then
 if (PBC(4).eq.2)testsystem_cuboctahedron = -1
 if (PBC(4).eq.0)testsystem_cuboctahedron = -2
endif

if (v(3).le.0.0) then
 if (PBC(5).eq.2)testsystem_cuboctahedron = -1
 if (PBC(5).eq.0)testsystem_cuboctahedron = -2
endif

if (v(3).gt.(float(dimz)*delta)) then
 if (PBC(6).eq.2)testsystem_cuboctahedron = -1
 if (PBC(6).eq.0)testsystem_cuboctahedron = -2
endif



if (testsystem_cuboctahedron.eq.0) then ! saves some time

do i = 1, 3 ! put into cell
if(v(i).lt.0.0) then
   if(PBC(2*i-1).eq.1)v(i) = PBCSYMR(v(i),maxx(i))
   if(PBC(2*i-1).eq.3)v(i) = PBCREFR(v(i),maxx(i))
endif
if(v(i).gt.dims(i)) then
   if(PBC(2*i).eq.1)v(i) = PBCSYMR(v(i),maxx(i))
   if(PBC(2*i).eq.3)v(i) = PBCREFR(v(i),maxx(i))
endif
enddo


! collision with particles


do j = 1, NNN

xx(:) = x(:) - Rell(:,j) ! distance to center of CO

v = xx
xx = MATMUL(MAT,v) ! to transformed space

! Looks for near neighbor, only important for PBC
if(PBC(1).eq.1)xx(1) = xx(1) - nint(xx(1) / (float(dimx)*delta)) * float(dimx)*delta
if(PBC(3).eq.1)xx(2) = xx(2) - nint(xx(2) / (float(dimy)*delta)) * float(dimy)*delta
if(PBC(5).eq.1)xx(3) = xx(3) - nint(xx(3) / (float(dimz)*delta)) * float(dimz)*delta
xx(:) = MATMUL(IMAT,xx) ! to real space


! Check in real space for collision with CO

call BetweenPlanes(plane1(:,j),klocta1(j),klocta1b(j),xx,test1)
call BetweenPlanes(plane2(:,j),klocta2(j),klocta2b(j),xx,test2)
call BetweenPlanes(plane3(:,j),klocta3(j),klocta3b(j),xx,test3)
call BetweenPlanes(plane4(:,j),klocta4(j),klocta4b(j),xx,test4)
call BetweenPlanes(plane5(:,j),klcubex1(j),klcubex2(j),xx,test5)
call BetweenPlanes(plane6(:,j),klcubey1(j),klcubey2(j),xx,test6)
call BetweenPlanes(plane7(:,j),klcubez1(j),klcubez2(j),xx,test7)
if(test1)then
  if(test2)then
    if(test3)then
      if(test4)then
        if(test5)then
          if(test6)then
            if(test7)then
               testsystem_cuboctahedron = -1
             return
           endif
         endif
       endif
     endif
   endif
 endif
endif

enddo

endif ! testsystem = 0

return

end function

