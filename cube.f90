subroutine update_matrix_cube(flag)
use system
use cube !cuidado agregar channel quizas
use ematrix
use MPI
use const
use chainsdat
use molecules
use channel
use transform, only : MAT, IMAT
use rotchain

implicit none

integer npoints ! points per cell for numerical integration 
integer counter
character*5 title
logical flag
integer j,ix,iy,iz
real*8 l_cubeL, l_cubeS
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


cutarea = 0.0 ! throw away cells that have less area than cutarea x area of the cell with largest area  
sumpolseg = 0.0
l_cubeL = l_cube + 2*delta
l_cubeS = l_cube - 2*delta

! clear all
voleps = 0.0
volprot = 0.0
volq = 0.0
volx = 0.0
volxx = 0.0
com = 0.0
ncha = 0

! channel center in x, y plane

 originc(1) = float(dimx)*delta/2.0 
 originc(2) = float(dimy)*delta/2.0 

 npoints = 50

 flag = .false.

 call integrate_cube(l_cubeL, c_cube,npoints, voleps1, sumvoleps1, flag)

 flag = .false. ! not a problem if eps lays outside boundaries

 
 call integrate_cube(l_cube,c_cube,npoints, volprot1, sumvolprot1, flag)

 call integrate_cube(l_cubeS,c_cube,npoints, volq1, sumvolq1, flag)

 call newintegrateg_cube(l_cube,c_cube,cubeR,l_pol,npoints,volx1,sumvolx1, com1, p1, ncha1, volxx1)

!! eps
 voleps1 = voleps1-volprot1
 voleps1 = voleps1*eepsc

!! charge
 volq1 = volprot1-volq1
 temp = sum(volq1)
 volq1 = volq1/temp*echargec/(delta**3) ! sum(volq) is echarge

area = 6.0*l_cube**2

!! volume
 volprot1 = volprot1 * 0.9999
 volprot = volprot+volprot1

! CHECK COLLISION HERE...
 if(maxval(volprot).gt.1.0) then ! collision
   flag=.true. 
 endif
 
 voleps = voleps + voleps1
 volq = volq + volq1 

! add com1 and volx to list

 volxx = volxx1

ncha = ncha1
do i = 1, ncha
volx(i)=volx1(i)
com(i,:)=com1(i,:)
p0(i,:)=p1(i,:)
rotangle(i) = 0.0
enddo

title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

sumpolseg = ncha

if (verbose.ge.2) then
temp = l_cube**3
!do j = 1, NNN
!temp = temp + 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)
!enddo
if (rank.eq.0) then
write(stdout,*) 'channel:', 'update_matrix: Total nanocube volumen real space= ', temp
write(stdout,*) 'channel:', 'update_matrix: Total discretized volumen =', (sum(volprot))*delta**3
write(stdout,*) 'channel:', 'number of polymers in system =', sumpolseg
write(stdout,*) 'channel:', 'surface area =', area
write(stdout,*) 'channel:', 'surface density =', sumpolseg/area
endif
endif

title = 'aveps'
counter = 1
!call savetodisk(voleps, title, counter)

title = 'avcha'
counter = 1
!call savetodisk(volq, title, counter)

title = 'avgrf'
counter = 1
call savetodisk(volxx, title, counter)

end subroutine

subroutine integrate_cube(l_cube,c_cube, npoints,volprot,sumvolprot,flag)
use system
use transform

implicit none
real*8 sumvolprot
integer npoints
real*8 l_cube
real*8 c_cube(3)
real*8 volprot(dimx,dimy,dimz)
real*8 dr(3), dxr(3)
integer ix,iy,iz,ax,ay,az
real*8 vect
logical flagin, flagout
real*8 intcell_cube
real*8 mmmult
integer jx,jy, jz
logical flag
integer RdimZ
real*8 box(4)
real*8 x(3), v(3)
integer xmin,xmax,ymin,ymax,zmin,zmax
integer i,j

logical flagsym
real*8 voltemp

volprot = 0.0
sumvolprot = 0.0 ! total volumen, including that outside the system

! scan over all cells

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz

flagin = .false.
flagout = .false.

do ax = 0,1
do ay = 0,1
do az = 0,1

! v in transformed space
v(1) = float(ax+ix-1)*delta
v(2) = float(ay+iy-1)*delta
v(3) = float(az+iz-1)*delta

! x in real space, v in transformed space
    x = MATMUL(IMAT,v)

x(1) = x(1) - c_cube(1)
x(2) = x(2) - c_cube(2)
x(3) = x(3) - c_cube(3)

if (((abs(x(1)).gt.(l_cube/2)).or.(abs(x(2)).gt.(l_cube/2))).or.(abs(x(3)).gt.(l_cube/2)))flagout=.true.
if (((abs(x(1)).lt.(l_cube/2))).and.(abs(x(2)).lt.(l_cube/2)).and.(abs(x(3)).lt.(l_cube/2)))flagin=.true.

enddo
enddo
enddo

if((flagin.eqv..true.).and.(flagout.eqv..false.)) then ! cell all inside channel
    voltemp = 1.0
endif
if((flagin.eqv..false.).and.(flagout.eqv..true.)) then ! cell all outside channel
    voltemp = 0.0
endif
if((flagin.eqv..true.).and.(flagout.eqv..true.)) then ! cell part inside annd outside channel
    voltemp = intcell_cube(l_cube, c_cube,ix,iy,iz,npoints)
endif

sumvolprot = sumvolprot + voltemp
volprot(ix,iy,iz) = voltemp

enddo ! ix
enddo ! iy
enddo ! iz


end subroutine

double precision function intcell_cube(l_cube,c_cube,ix,iy,iz,n)
use system
use transform

implicit none
real*8 l_cube
real*8 c_cube(3)
integer ix,iy,iz,ax,ay,az
integer cc
real*8 vect
integer n
real*8 mmmult
real*8 dr(3), dxr(3)

cc = 0
do ax = 1, n
do ay = 1, n
do az = 1, n

dr(1) = ix*delta-(ax)*delta/float(n) 
dr(2) = iy*delta-(ay)*delta/float(n) 
dr(3) = iz*delta-(az)*delta/float(n) 

! dr in transformed space
dxr = MATMUL(IMAT, dr)

dxr(1) = dxr(1) - c_cube(1)
dxr(2) = dxr(2) - c_cube(2)
dxr(3) = dxr(3) - c_cube(3)

if (((abs(dxr(1)).lt.(l_cube/2)).and.(abs(dxr(2)).lt.(l_cube/2))).and.(abs(dxr(3)).lt.(l_cube/2)))cc=cc+1 ! integra dentro del cubo

enddo
enddo
enddo

intcell_cube  = float(cc)/(float(n)**3)
end function

subroutine newintegrateg_cube(l_cube,c_cube,cubeR,l_pol,npoints,volx1,sumvolx1,com1,p1,ncha1,volxx1)
use system
use transform
use chainsdat
use ematrix
use const
implicit none
real*8 sumvolx1
integer l_pol, npoints
integer indexvolx(dimx,dimy,dimz)
integer listvolx(ncha,3)
real*8 sep !separacion entre polimeros en una cara
real*8 l_cube, c_cube(3)
integer cubeR
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
integer jjjz, jjjt, npointz, npointt
integer RdimZ

pi=acos(-1.0)

dims(1) = dimx
dims(2) = dimy
dims(3) = dimz

indexvolx = 0
ncha1 = 0
volx1 = 0.0
sumvolx1 = 0.0 ! total volumen, including that outside system
com1 = 0.0
p1 = 0
volxx1 = 0.0

! This routine determines the surface coverage and grafting positions only for cube
!

if (cubeR.eq.0) then

sep = l_cube/float(l_pol)

do ix =1,l_pol
do iy = 1,l_pol
        
! Cara de ABAJO

x(1) = c_cube(1) - l_cube/2.0 + float(ix)*sep - sep/2.0
x(2)= c_cube(2) - l_cube/2.0 + float(iy)*sep - sep/2.0
x(3) = c_cube(3) - l_cube/2.0

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

! increase counter
 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
sumvolx1 = sumvolx1 + 1.0
com1(ncha1,:) = x(:)
com1(ncha1,3) = com1(ncha1,3) - lseg/2.0

! Cara de ARRIBA

x(1) = c_cube(1) - l_cube/2.0 + float(ix)*sep -  sep/2.0
x(2)= c_cube(2) - l_cube/2.0 + float(iy)*sep -  sep/2.0
x(3) = c_cube(3) + l_cube/2.0

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

! increase counter
 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
sumvolx1 = sumvolx1 + 1.0
com1(ncha1,:) = x(:)
com1(ncha1,3) = com1(ncha1,3) + lseg/2.0

! Cara de IZQUIERDA

x(1) = c_cube(1) - l_cube/2.0 + float(ix)*sep - sep/2.0
x(2)= c_cube(2) - l_cube/2.0
x(3) = c_cube(3) - l_cube/2.0 + float(iy)*sep -  sep/2.0

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

! increase counter
 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
sumvolx1 = sumvolx1 + 1.0
com1(ncha1,:) = x(:)
com1(ncha1,2) = com1(ncha1,2) - lseg/2.0

! Cara de DERECHA zx

x(1) = c_cube(1) - l_cube/2.0 + float(ix)*sep - sep/2.0
x(2)= c_cube(2) + l_cube/2.0
x(3) = c_cube(3) - l_cube/2.0 + float(iy)*sep - sep/2.0

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

!increase counter
 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
sumvolx1 = sumvolx1 + 1.0
com1(ncha1,:) = x(:)
com1(ncha1,2) = com1(ncha1,2) + lseg/2.0

! Cara de ATRAS

x(1) = c_cube(1) - l_cube/2.0
x(2)= c_cube(2) - l_cube/2.0 + float(ix)*sep -  sep/2.0
x(3) = c_cube(3) - l_cube/2.0 + float(iy)*sep - sep/2.0

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

 !increase counter
 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
sumvolx1 = sumvolx1 + 1.0
com1(ncha1,:) = x(:)
com1(ncha1,1) = com1(ncha1,1) - lseg/2.0

! Cara de ADELANTE

x(1) = c_cube(1) + l_cube/2.0
x(2)= c_cube(2) - l_cube/2.0 + float(ix)*sep - sep/2.0
x(3) = c_cube(3) - l_cube/2.0 + float(iy)*sep - sep/2.0

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

! increase counter
 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
sumvolx1 = sumvolx1 + 1.0
com1(ncha1,:) = x(:)
com1(ncha1,1) = com1(ncha1,1) + lseg/2.0

enddo
enddo

endif

if (cubeR.eq.1) then

sep = l_cube/float(l_pol)

do ix =1,l_pol
do iy = 1,l_pol

! Cara de ARRIBA

x(1) = c_cube(1) - l_cube/2.0 + float(ix)*sep -  sep/2.0
x(2)= c_cube(2) - l_cube/2.0 + float(iy)*sep -  sep/2.0
x(3) = c_cube(3) + l_cube/2.0

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

! increase counter
 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
sumvolx1 = sumvolx1 + 1.0
com1(ncha1,:) = x(:)
com1(ncha1,3) = com1(ncha1,3) + lseg/2.0

! Cara de DERECHA zx

x(1) = c_cube(1) - l_cube/2.0 + float(ix)*sep - sep/2.0
x(2)= c_cube(2) + l_cube/2.0
x(3) = c_cube(3) - l_cube/2.0 + float(iy)*sep - sep/2.0

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

!increase counter
 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
sumvolx1 = sumvolx1 + 1.0
com1(ncha1,:) = x(:)
com1(ncha1,2) = com1(ncha1,2) + lseg/2.0

! Cara de ADELANTE

x(1) = c_cube(1) + l_cube/2.0
x(2)= c_cube(2) - l_cube/2.0 + float(ix)*sep - sep/2.0
x(3) = c_cube(3) - l_cube/2.0 + float(iy)*sep - sep/2.0

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

! increase counter
 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
sumvolx1 = sumvolx1 + 1.0
com1(ncha1,:) = x(:)
com1(ncha1,1) = com1(ncha1,1) + lseg/2.0

enddo
enddo

endif

if (cubeR.eq.2) then

sep = 2.0*(l_cube/((l_pol*2.0)-1.0))

! Cara de ARRIBA

do ix =1,l_pol
do iy = 1,l_pol

x(1) = c_cube(1) - l_cube/2.0 + float(ix)*sep - sep/2.0
x(2)= c_cube(2) - l_cube/2.0 + float(iy)*sep - sep/2.0
x(3) = c_cube(3) + l_cube/2.0

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

! increase counter
 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
sumvolx1 = sumvolx1 + 1.0
com1(ncha1,:) = x(:)
com1(ncha1,3) = com1(ncha1,3) + lseg/2.0

enddo
enddo

!Cara de DERECHA zx

do ix =1,l_pol
do iy = 1,l_pol - 1

x(1) = c_cube(1) - l_cube/2.0 + float(ix)*sep - sep/2.0
x(2)= c_cube(2) + l_cube/2.0
x(3) = c_cube(3) - l_cube/2.0 + float(iy)*sep - sep/2.0

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

!increase counter
 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
sumvolx1 = sumvolx1 + 1.0
com1(ncha1,:) = x(:)
com1(ncha1,2) = com1(ncha1,2) + lseg/2.0

enddo
enddo

! Cara de ADELANTE

do ix =1,l_pol - 1
do iy = 1,l_pol - 1

x(1) = c_cube(1) + l_cube/2.0
x(2)= c_cube(2) - l_cube/2.0 + float(ix)*sep - sep/2.0
x(3) = c_cube(3) - l_cube/2.0 + float(iy)*sep - sep/2.0

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

!increase counter
 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
sumvolx1 = sumvolx1 + 1.0
com1(ncha1,:) = x(:)
com1(ncha1,1) = com1(ncha1,1) + lseg/2.0

enddo
enddo

endif

end subroutine
