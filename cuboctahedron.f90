subroutine update_matrix_cuboctahedron(flag)
use system
use ellipsoid !cuidado agregar channel quizas
use ematrix
use MPI
use const
use chainsdat
use molecules
!use channel
use transform, only : MAT, IMAT
use rotchain

implicit none

integer npoints ! points per cell for numerical integration 
integer counter
character*5 title
logical flag
integer j,ix,iy,iz
real*8 lcubeL, lcubeS, loctaL, loctaS
real*8 center(3)
real pnumber
real*8 area
real*8 sumpolseg 
real*8 sstemp,vvtemp, maxss
real*8 cutarea
real*8 temp
real*8 temp2
real*8 sumvoleps1, sumvolprot1, sumvolq1, sumvolx1
integer ncha1
real*8 volx1(maxvolx)
real*8 com1(maxvolx,3)
integer p1(maxvolx,3)
integer i
real*8 volxx1(dimx,dimy,dimz)
real*8 volxx(dimx,dimy,dimz)
real*8 x(3), v(3), hcyl
integer nbands
real*8 COvol

sumpolseg = 0.0
cutarea = 0.0 ! throw away cells that have less area than cutarea x area of the cell with largest area  


! clear all (variables for the sum from all CO)
voleps = 0.0
volprot = 0.0
volq = 0.0
volx = 0.0
volxx = 0.0
com = 0.0
ncha = 0

do j=1,NNN ! loop ovr the particles

lcubeL = Lcubell(j) + delta 
lcubeS = Lcubell(j) - delta
loctaL = Loctall(j) + delta
loctaS = Loctall(j) - delta

 flag = .false.

 !change center 

 center(1) = Rell(1,j) 
 center(2) = Rell(2,j)
 center(3) = Rell(3,j)
 
 npoints = 50
 call integrate_cuboctahedron(lcubeL,loctaL,center,rotmatCO(:,:,j),j,npoints,voleps1,sumvoleps1,flag)

 flag = .false. ! not a problem if eps lays outside boundaries
 
 npoints = 50
 call integrate_cuboctahedron(Lcubell(j),Loctall(j),center,rotmatCO(:,:,j),j,npoints,volprot1,sumvolprot1,flag)

 npoints = 50
 call integrate_cuboctahedron(lcubeS,loctaS,center,rotmatCO(:,:,j),j,npoints,volq1,sumvolq1,flag)

 npoints = 200
 call newintegrateg_cuboctahedron(Lcubell(j),Loctall(j),center,rotmatCO(:,:,j),j,npoints,volx1,sumvolx1,com1,p1,ncha1,volxx1)

!! eps
 voleps1 = voleps1-volprot1
 voleps1 = voleps1*eeps(j)

!! charge
 volq1 = volprot1-volq1
 temp = sum(volq1)

 if(temp.ne.0.0)volq1 = volq1/temp*echarge(j)/(delta**3) ! sum(volq) is echarge

 area =  3.0**(1.0/2.0)*Loctall(j)**2 + 3.0*(1.0-3.0**(1.0/2.0))*(Loctall(j) - Lcubell(j))**2.0 ! MARIO

!! Normalize volx1 and volxx1 so that the integral is equal to the total number of ligands on the CO j

 volx1 = volx1/sumvolx1*area*sigma(j)
 volxx1 = volxx1/sumvolx1*area*sigma(j)

!! Sum of number of polymer segments

 sumpolseg = sumpolseg + area*sigma(j)*long

!! volume
 volprot1 = volprot1 * 0.99
 volprot = volprot+volprot1

! CHECK COLLISION HERE...
 if(maxval(volprot).gt.1.0) then ! collision
   flag=.true. 
 endif
 
 voleps = voleps + voleps1
 volq = volq + volq1 

! add com1 and volx to list

 volxx = volxx + volxx1 !actualizo volxx, verificar

do i = 1, ncha1
ncha = ncha + 1
volx(ncha)=volx1(i)
com(ncha,:)=com1(i,:)
p0(ncha,:)=p1(i,:)
enddo

enddo ! NNN

! save total CO density
title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

!stop

! print information summary 
if (verbose.ge.2) then

!temp = volume of cuboctahedron
COvol = 0.0
do j = 1, NNN
   COvol = COvol + (1.0/6.0)*Loctall(j)**3 - 0.5*(Loctall(j) -Lcubell(j))**3
enddo

if (rank.eq.0) then
write(stdout,*) 'cuboctahedron:', 'Total nanocuboct volumen real space= ', COvol
write(stdout,*) 'cuboctahedron:', 'Total discretized volumen =', (sum(volprot))*delta**3
write(stdout,*) 'cuboctahedron:', 'total number of segments in system =', sumpolseg
endif
endif

do j = 1, NNN
 area = 3.0**(1.0/2.0)*Loctall(j)**2 + 3.0*(1.0-3.0**(1.0/2.0))*(Loctall(j) - Lcubell(j))**2.0 ! MARIO
! area = 3.0**(1.0/2.0)*Loctall(j)**2 + 6.0*(1.0-2.0**(1.0/2.0))*(Loctall(j) - Lcubell(j))**2.0 ! LEO
 if (rank.eq.0) write(stdout,*) 'cuboctahedron:', ' Total nanocuboct #',j,' area', area
enddo
!

title = 'aveps'
counter = 1
call savetodisk(voleps, title, counter)

title = 'avcha'
counter = 1
call savetodisk(volq, title, counter)

title = 'avgrf'
counter = 1
call savetodisk(volxx, title, counter)

end subroutine


double precision function intcell_cuboctahedron(lcube,locta,center,rotmatrix,npart,ix,iy,iz,n)
use system
use transform
use COrotMod

implicit none
real*8 lcube,locta
real*8 center(3)
real*8 rotmatrix(3,3)
integer npart
integer ix,iy,iz,ax,ay,az
integer cc
real*8 vect
integer n
real*8 mmmult
real*8 dr(3), dxr(3)
logical test1, test2, test3, test4
logical test5, test6, test7

cc = 0
do ax = 1, n
do ay = 1, n
do az = 1, n

dr(1) = ix*delta-(ax)*delta/float(n) 
dr(2) = iy*delta-(ay)*delta/float(n) 
dr(3) = iz*delta-(az)*delta/float(n) 

! dr in transformed space
dxr = MATMUL(IMAT, dr)

dxr(1) = dxr(1) - center(1)
dxr(2) = dxr(2) - center(2)
dxr(3) = dxr(3) - center(3)

call BetweenPlanes(plane1(:,npart),klocta1(npart),klocta1b(npart),dxr,test1)
call BetweenPlanes(plane2(:,npart),klocta2(npart),klocta2b(npart),dxr,test2)
call BetweenPlanes(plane3(:,npart),klocta3(npart),klocta3b(npart),dxr,test3)
call BetweenPlanes(plane4(:,npart),klocta4(npart),klocta4b(npart),dxr,test4)
call BetweenPlanes(plane5(:,npart),klcubex1(npart),klcubex2(npart),dxr,test5)
call BetweenPlanes(plane6(:,npart),klcubey1(npart),klcubey2(npart),dxr,test6)
call BetweenPlanes(plane7(:,npart),klcubez1(npart),klcubez2(npart),dxr,test7)
if(test1)then
 if(test2)then
  if(test3)then
   if(test4)then
    if(test5)then
     if(test6)then
      if(test7)then
       cc = cc + 1
      endif
     endif
    endif
   endif
  endif
 endif
endif

enddo
enddo
enddo

intcell_cuboctahedron  = float(cc)/(float(n)**3)
end function

subroutine newintegrateg_cuboctahedron(lcube,locta,center,rotmatrix,npart,npoints,volx1,sumvolx1,com1,p1,ncha1,volxx1)
use system
use transform
use chainsdat
use ematrix
use const
use COrotMod

implicit none
real*8 center(3)
real*8 rotmatrix(3,3)
integer npart
real*8 sumvolx1
integer npoints
!real*8 AAA(3,3), AAAX(3,3)
integer indexvolx(dimx,dimy,dimz)
integer listvolx(ncha,3)
real*8 Rell(3), Aell(3)
real*8 radio
real*8 phi, dphi, tetha,dtetha, as, ds
integer mphi, mtetha
integer ix,iy,iz,jx,jy,jz
real*8 x(3), v(3)
integer i,j
integer ncount
real*8 comshift ! how far from the surface of the sphere the grafting point is
integer ncha1 ! count for current sphere
real*8 volx1(maxvolx)
real*8 com1(maxvolx,3)
integer p1(maxvolx,3)
real*8 volxx1(dimx,dimy,dimz)
integer flagin
integer dims(3), is(3), js(3)
real*8 lcuber, pasoc, pasoo
real*8 xx, yy, zz
real*8 vector(3)
real*8 lcube, locta
real*8 xtotest(3)
logical test1, test2, test3, test4
logical test5, test6, test7

indexvolx = 0
ncha1 = 0
volx1 = 0.0
sumvolx1 = 0.0 ! total volumen, including that outside system
com1 = 0.0
p1 = 0
volxx1 = 0.0

! This routine determines the surface coverage and grafting positions only for cuboctahedron!
lcuber = lcube/locta

pasoc = delta/float(npoints)/locta

pasoo = pasoc*(1.0/3.0)**(1.0/4.0) ! % different integration steps are needed to have the same area element
                          ! % Lets R(u,v) = x(u,v)i + y(u,v)j + z(u,v)k
                          ! % then dA = dR/du x dR/dv * du*dv
                          ! % For octahedro x = u, y = v, z = u+v-1
                          ! % For cube x = v-1, y = 0, z = 0  

!%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%% OCTAHEDRO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx = 0.0

do while (xx < 1.0)
    yy = 0.0
    do while (yy < (1.0-xx))
    zz = -xx-yy+1 

    xtotest(1) = xx*locta/2.0
    xtotest(2) = yy*locta/2.0
    xtotest(3) = zz*locta/2.0
    xtotest = MATMUL(rotmatrix,xtotest)
   
    call BetweenPlanes(plane5(:,npart),klcubex1(npart),klcubex2(npart),xtotest,test5)
    call BetweenPlanes(plane6(:,npart),klcubey1(npart),klcubey2(npart),xtotest,test6)
    call BetweenPlanes(plane7(:,npart),klcubez1(npart),klcubez2(npart),xtotest,test7)
   
    if((test5).and.(test6).and.(test7))then
       x(1) = xx
       x(2) = yy
       x(3) = -xx-yy+1
       x = MATMUL(rotmatrix,x)
       call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

       x(1) = -xx
       x(2) = -yy
       x(3) = -xx-yy+1
       x = MATMUL(rotmatrix,x)
       call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
        
       x(1) = -xx
       x(2) = yy
       x(3) = -xx-yy+1
       x = MATMUL(rotmatrix,x)
       call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
        
       x(1) = xx
       x(2) = -yy
       x(3) = -xx-yy+1
       x = MATMUL(rotmatrix,x)
       call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

       x(1) = xx
       x(2) = yy
       x(3) = xx+yy-1
       x = MATMUL(rotmatrix,x)
       call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
        
       x(1) = -xx
       x(2) = -yy
       x(3) = +xx+yy-1
       x = MATMUL(rotmatrix,x)
       call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
        
       x(1) = -xx
       x(2) = yy
       x(3) = +xx+yy-1
       x = MATMUL(rotmatrix,x)
       call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
        
       x(1) = xx
       x(2) = -yy
       x(3) = +xx+yy-1 
       x = MATMUL(rotmatrix,x)
       call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      
     
     endif  
        yy = yy + pasoo
    enddo    
    xx = xx + pasoo
enddo


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%% CUBO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xx = -lcuber
do while (xx < lcuber)
  yy = -lcuber
  do while (yy < lcuber)

  zz = lcuber    
  
  xtotest(1) = xx*locta/2.0
  xtotest(2) = yy*locta/2.0
  xtotest(3) = zz*locta/2.0

  xtotest = MATMUL(rotmatrix,xtotest)
    
  call BetweenPlanes(plane1(:,npart),klocta1(npart),klocta1b(npart),xtotest,test1)
  call BetweenPlanes(plane2(:,npart),klocta2(npart),klocta2b(npart),xtotest,test2)
  call BetweenPlanes(plane3(:,npart),klocta3(npart),klocta3b(npart),xtotest,test3)
  call BetweenPlanes(plane4(:,npart),klocta4(npart),klocta4b(npart),xtotest,test4)

  if((test1).and.(test2).and.(test3).and.(test4))then 
        x(1) = -lcuber
        x(2) = xx
        x(3) = yy
        x = MATMUL(rotmatrix,x)
        call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

        x(1) = lcuber
        x(2) = xx
        x(3) = yy
        x = MATMUL(rotmatrix,x)
        call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

        x(1) = xx
        x(2) = -lcuber
        x(3) = yy
        x = MATMUL(rotmatrix,x)
        call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

        x(1) = xx
        x(2) = lcuber
        x(3) = yy
        x = MATMUL(rotmatrix,x)
        call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

        x(1) = xx
        x(2) = yy
        x(3) = -lcuber
        x = MATMUL(rotmatrix,x)
        call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)      

        x(1) = xx
        x(2) = yy
        x(3) = lcuber
        x = MATMUL(rotmatrix,x)
        call integrar_matrices(x,center,locta,indexvolx,ncha1,p1,volxx1,volx1,com1,sumvolx1)     
    endif
        
  yy = yy + pasoc
  enddo  
  xx = xx + pasoc
enddo 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aca termina de tirar puntos en la superficie del cubooctahedro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do i = 1, ncha1
com1(i,:) = com1(i,:)/volx1(i)
! Moves the position of the first segment lseg/2 away from the surface to prevent collision due to round errors.
vector(:) = com1(i,:)-center(:)
vector(:) = vector(:)/norm2(vector)
com1(i,:) = com1(i,:) + 1.5*lseg*vector(:)
enddo

end


subroutine integrar_matrices(x, center, locta, indexvolx, ncha1, p1, volxx1, volx1, com1, sumvolx1)
use system
use const
use ematrix
use transform, only : MAT, IMAT
implicit none
real*8 x(3), center(3)
integer flagin
integer j
integer is(3), js(3), dims(3)
integer, external :: PBCSYMI
integer jx,jy,jz
integer ncha1
integer indexvolx(dimx,dimy,dimz)
integer p1(maxvolx,3)
real*8 volxx1(dimx,dimy,dimz)
real*8 com1(maxvolx,3)
real*8 volx1(maxvolx)
real*8 v(3)
real*8 locta, lcube
real*8 sumvolx1

x(:) = x(:)*locta/2.0 + center(:)

dims(1) = dimx
dims(2) = dimy
dims(3) = dimz

! x in real space, v in transformed space
    v = MATMUL(MAT,x)

! PBC

flagin = 1

do j = 1,3

    is(j) = floor(v(j)/delta)+1
    js(j) = is(j)

select case (PBC((j-1)*2+1))
  case (0 , 2)
    if(is(j).lt.1) then
    write(stdout,*) 'Error in newintegrateg: out of boundary'
    endif
  case (1)
    js(j)=PBCSYMI(is(j), dims(j)) 
  case (3)
    if(v(j).lt.0.0)flagin=0
endselect

select case (PBC((j-1)*2+2))
  case (0 , 2)
    if(is(j).gt.dims(j)) then
    write(stdout,*) 'Error in newintegrateg: out of boundary'
    endif
  case (1)
    js(j)=PBCSYMI(is(j), dims(j)) 
  case (3)
    if(v(j).gt.float(dims(j))*delta)flagin=0
endselect
enddo

jx = js(1)
jy = js(2)
jz = js(3)

if(flagin.eq.1) then

! increase counter
if(indexvolx(jx,jy,jz).eq.0) then

 if(ncha1.eq.maxvolx) then
   write(stdout,*) 'ellipsoid: increase maxvolx'
   stop
 endif

 ncha1 = ncha1 + 1
 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz
endif

! agrega el punto
volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
com1(indexvolx(jx,jy,jz),:) = com1(indexvolx(jx,jy,jz),:) + x(:)
endif

sumvolx1 = sumvolx1 + 1.0

end



subroutine integrate_cuboctahedron(lcube,locta,center,rotmatrix,npart,npoints,volprot,sumvolprot,flag)
use system
use transform
use const
use COrotMod
implicit none
real*8 sumvolprot
integer npoints
integer npart
real*8 lcube, locta
real*8 center(3)
real*8 rotmatrix(3,3)
real*8 volprot(dimx,dimy,dimz)
real*8 dr(3), dxr(3)
integer ix,iy,iz,ax,ay,az
real*8 vect
logical flagin, flagout
real*8 intcell_cuboctahedron
real*8 mmmult
integer jx,jy, jz
logical flag
integer RdimZ
real*8 box(4)
real*8 x(3), v(3)
integer xmin,xmax,ymin,ymax,zmin,zmax
integer i,j
real*8 maxAell
logical flagsym
real*8 voltemp
real*8 Rpos(3)
logical test1, test2, test3, test4
logical test5, test6, test7

volprot = 0.0
sumvolprot = 0.0 ! total volumen, including that outside the system
maxAell = locta/2.0 ! maximum size CO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create a box in transformed coordinate enclosing the CO
Rpos(1) = center(1)
Rpos(2) = center(2)
Rpos(3) = center(3)

!!!!!!!!!!!!!!!! xmin !!!!!!!!!!!!!!!
x(1) = Rpos(1)-maxAell
do j = 1, 2
do i = 1, 2
x(2) = Rpos(2)+maxAell*(-1.0)**j
x(3) = Rpos(3)+maxAell*(-1.0)**i
v = MATMUL(MAT,x)
box(i+2*(j-1)) = v(1)
enddo
enddo
xmin = int(minval(box)/delta)-2
!!!!!!!!!!!!!! xmax !!!!!!!!!!!!!!!!!!!!!1
x(1) = Rpos(1)+maxAell
do j = 1, 2
do i = 1, 2
x(2) = Rpos(2)+maxAell*(-1.0)**j
x(3) = Rpos(3)+maxAell*(-1.0)**i
v = MATMUL(MAT,x)
box(i+2*(j-1)) = v(1)
enddo
enddo
xmax = int(maxval(box)/delta)+2
!!!!!!!!!!!!!!!! ymin !!!!!!!!!!!!!!!
x(2) = Rpos(2)-maxAell
do j = 1, 2
do i = 1, 2
x(1) = Rpos(1)+maxAell*(-1.0)**j
x(3) = Rpos(3)+maxAell*(-1.0)**i
v = MATMUL(MAT,x)
box(i+2*(j-1)) = v(2)
enddo
enddo
ymin = int(minval(box)/delta)-2
!!!!!!!!!!!!!!!! ymax !!!!!!!!!!!!!!!
x(2) = Rpos(2)+maxAell
do j = 1, 2
do i = 1, 2
x(1) = Rpos(1)+maxAell*(-1.0)**j
x(3) = Rpos(3)+maxAell*(-1.0)**i
v = MATMUL(MAT,x)
box(i+2*(j-1)) = v(2)
enddo
enddo
ymax = int(maxval(box)/delta)+2
!!!!!!!!!!!!!!!! zmin !!!!!!!!!!!!!!!
x(3) = Rpos(3)-maxAell
do j = 1, 2
do i = 1, 2
x(1) = Rpos(1)+maxAell*(-1.0)**j
x(2) = Rpos(2)+maxAell*(-1.0)**i
v = MATMUL(MAT,x)
box(i+2*(j-1)) = v(3)
enddo
enddo
zmin = int(minval(box)/delta)-2
!!!!!!!!!!!!!!!! zmax !!!!!!!!!!!!!!!
x(3) = Rpos(3)+maxAell
do j = 1, 2
do i = 1, 2
x(1) = Rpos(1)+maxAell*(-1.0)**j
x(2) = Rpos(2)+maxAell*(-1.0)**i
v = MATMUL(MAT,x)
box(i+2*(j-1)) = v(3)
enddo
enddo
zmax = int(maxval(box)/delta)+2

! Make a list of the cells that have no CO, those that have part CO and those that have full CO
! Consider boundary conditions

do ix = xmin, xmax
do iy = ymin, ymax
do iz = zmin, zmax

jx=ix
jy=iy
jz=iz

if(PBC(1).eq.1) then
 jx=mod(ix+dimx-1,dimx)+1
endif
if(PBC(3).eq.1) then
 jy=mod(iy+dimy-1,dimy)+1
endif
if(PBC(5).eq.1) then
 jz=mod(iz+dimz-1,dimz)+1
endif

flagin=.false.
flagout=.false.

do ax = 0,1
do ay = 0,1
do az = 0,1

! v in transformed space
v(1) = float(ax+ix-1)*delta
v(2) = float(ay+iy-1)*delta
v(3) = float(az+iz-1)*delta

! x in real space, v in transformed space
    x = MATMUL(IMAT,v)

x(1) = x(1) - center(1)
x(2) = x(2) - center(2)
x(3) = x(3) - center(3)

call BetweenPlanes(plane1(:,npart),klocta1(npart),klocta1b(npart),x,test1)
call BetweenPlanes(plane2(:,npart),klocta2(npart),klocta2b(npart),x,test2)
call BetweenPlanes(plane3(:,npart),klocta3(npart),klocta3b(npart),x,test3)
call BetweenPlanes(plane4(:,npart),klocta4(npart),klocta4b(npart),x,test4)
call BetweenPlanes(plane5(:,npart),klcubex1(npart),klcubex2(npart),x,test5)
call BetweenPlanes(plane6(:,npart),klcubey1(npart),klcubey2(npart),x,test6)
call BetweenPlanes(plane7(:,npart),klcubez1(npart),klcubez2(npart),x,test7)
if(test1)then
 if(test2)then
  if(test3)then
   if(test4)then
    if(test5)then
     if(test6)then
      if(test7)then
        flagin=.true.
      else
       flagout=.true.
      endif
     else
      flagout=.true.
     endif
    else
     flagout=.true.
    endif
   else
    flagout=.true.
   endif
  else
   flagout=.true.
  endif
 else
  flagout=.true.
 endif
else 
 flagout=.true.
endif

enddo
enddo
enddo

!! Check particle out of system


    flagsym = .false.
    if (jx.lt.1) then
       if(PBC(1).ne.3) then
         write(stdout,*) 'cuboctahedron:','update_matrix: ix', ix
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jy.lt.1) then
       if(PBC(3).ne.3) then
         write(stdout,*) 'cuboctahedron:','update_matrix: iy', iy
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jz.lt.1) then
       if(PBC(5).ne.3) then
         write(stdout,*) 'cuboctahedron:','update_matrix: iz', iz
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jx.gt.dimx) then
       if(PBC(2).ne.3) then
         write(stdout,*) 'cuboctahedron:','update_matrix: ix', ix
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jy.gt.dimy) then
       if(PBC(4).ne.3) then
         write(stdout,*) 'cuboctahedron:','update_matrix: iy', iy
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jz.gt.dimz) then
       if(PBC(6).ne.3) then
         write(stdout,*) 'cuboctahedron:','update_matrix: iz', iz
         stop
       else
         flagsym = .true.
       endif
    endif

if(flagsym.eqv..false.) then  ! cell is not out of system due to reflection symmetry
if((flagin.eqv..true.).and.(flagout.eqv..false.)) then ! cell all inside channel
    voltemp = 1.0
endif
if((flagin.eqv..false.).and.(flagout.eqv..true.)) then ! cell all outside channel
    voltemp = 0.0
endif
if((flagin.eqv..true.).and.(flagout.eqv..true.)) then ! cell part inside annd outside channel
    voltemp = intcell_cuboctahedron(lcube,locta,center,rotmatrix,npart,ix,iy,iz,npoints)
endif

sumvolprot = sumvolprot + voltemp
volprot(jx,jy,jz) = voltemp
endif ! flagsym


enddo ! ix
enddo ! iy
enddo ! iz

end subroutine


subroutine BetweenPlanes(plane,s1,s2,point,test)

implicit none

logical test
real*8 point(3)
real*8 plane(3)
real*8 s1, s2, h, d

h = s2-s1
d = DOT_PRODUCT(point,plane) - s1

if((abs(d).lt.abs(h)).and.((d*h).gt.0.0)) then
   test = .true.
else
   test = .false.
endif

end subroutine BetweenPlanes
