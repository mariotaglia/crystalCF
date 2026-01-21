
subroutine update_matrix_geom(flag)
use system
use ellipsoid
use geometries
use ematrix
use MPI
use const
use chainsdat
use molecules
use rotchain
use ellipsoid_create
use cuboctahedron_create
implicit none
integer npoints ! points per cell for numerical integration 
integer counter
character*5 title
logical flag
integer j
real pnumber
real*8 area_ell
real*8 sumpolseg_ell
real*8 maxss
real*8 cutarea_ell
real*8 temp
real*8 temp2
real*8 sumvoleps1, sumvolprot1, sumvolx1
integer ncha1
real*8 volx1(maxvolx)
real*8 com1(maxvolx,3)
integer p1(maxvolx,3)
integer i
real*8 volxx1(dimx,dimy,dimz)
real*8 volxx(dimx,dimy,dimz)

real*8 lcubeL, lcubeS, loctaL, loctaS
real*8 center_co(3)
real*8 area_co
real*8 sumpolseg_co
real*8 cutarea_co
real*8 COvol

call make_ellipsoid(NNNell) ! update matrixes for all particles

!cutarea = 0.01 ! throw away cells that have less area than cutarea x area of the cell with largest area  
cutarea_ell = 0.0 ! throw away cells that have less area than cutarea x area of the cell with largest area  
sumpolseg_ell = 0.0

! clear all
voleps = 0.0
volprot = 0.0
volx = 0.0
volxx = 0.0
com = 0.0
ncha = 0

do j = 1, NNNell
   seed = seed_np
  ! rotate ellipsoid matrixes according to current rotation matrix

   call rotv(AAA(:,:,j), rotmatrix(:,:,ids_ell(j)))
   call rotv(AAAS(:,:,j), rotmatrix(:,:,ids_ell(j)))
   call rotv(AAAL(:,:,j), rotmatrix(:,:,ids_ell(j)))
   call rotv(AAAX(:,:,j), rotmatrix(:,:,ids_ell(j)))

   call rotvo(orient(:,j), rotmatrix(:,:,ids_ell(j)))

   npoints = 50

   flag = .false.

   call integrate_ellipsoid(AAAL(:,:,j),AellL(:,j), Rell(:,ids_ell(j)),npoints, voleps1 , sumvoleps1, flag)
   flag = .false. ! not a problem if eps lays outside boundaries
   call integrate_ellipsoid(AAA(:,:,j),Aell(:,j), Rell(:,ids_ell(j)),npoints, volprot1, sumvolprot1, flag)

   npoints = 100000000
   call newintegrateg_ellipsoid(Aell(:,j),Rell(:,ids_ell(j)),npoints,volx1,sumvolx1, com1, p1, ncha1, volxx1)

  !! volume
   temp = 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)/(sumvolprot1*delta**3) ! rescales volume
   volprot1 = volprot1*temp                                                 ! WARNING: transformation should mantain cell volumen
   sumvolprot1 = sumvolprot1*temp

  !! eps
   voleps1 = voleps1-volprot1
   voleps1 = voleps1*eeps(ids_ell(j))

  !! grafting
   pnumber = 1.6075

   area_ell = (Aell(1,j)*Aell(2,j))**pnumber 
   area_ell = area_ell+(Aell(1,j)*Aell(3,j))**pnumber 
   area_ell = area_ell+(Aell(2,j)*Aell(3,j))**pnumber 
   area_ell = 4.0*pi*(area_ell/3.0)**(1.0/pnumber) ! approximate (< 1% error) area of elipsoid, see wikipedia

   temp2 = maxval(volx1)

   where(volx1<temp2*cutarea_ell) ! remove cells with very little area
   volx1 = 0.0
   end where 

   volx1 = volx1/sumvolx1*area_ell*sigma(ids_ell(j))
   volxx1 = volxx1/sumvolx1*area_ell*sigma(ids_ell(j))

   maxss = 1.0d100

   sumpolseg_ell = sumpolseg_ell + area_ell*sigma(ids_ell(j))*longp(ids_ell(j))

  !! volume  
   volprot1 = volprot1 * 0.99
   volprot = volprot+volprot1

  ! CHECK COLLISION HERE...
   if(maxval(volprot).gt.1.0) then ! collision
     flag=.true. 
     exit
   endif
   
   voleps = voleps + voleps1

  ! add com1 and volx to list

   volxx = volxx + volxx1

do i = 1, ncha1
  ncha = ncha+1
  volx(ncha)=volx1(i)
  com(ncha,:)=com1(i,:)
  p0(ncha,:)=p1(i,:)
  longc(ncha) = longp(ids_ell(j))
enddo

end do 

!===========================================

do j=1,NNNco ! loop ovr the cuboctahedron particles
  seed = seed_np
  lcubeL = Lcubell(j) + delta 
  lcubeS = Lcubell(j) - delta
  loctaL = Loctall(j) + delta
  loctaS = Loctall(j) - delta

   flag = .false.

   !change center 

   center_co(1) = Rell(1,ids_co(j)) 
   center_co(2) = Rell(2,ids_co(j))
   center_co(3) = Rell(3,ids_co(j))

   npoints = 50
   call integrate_cuboctahedron(lcubeL,loctaL,center_co,rotmatrix(:,:,ids_co(j)), j, &
    npoints,voleps1,sumvoleps1,flag)

   flag = .false. ! not a problem if eps lays outside boundaries
   
   npoints = 50
   call integrate_cuboctahedron(Lcubell(j),Loctall(j),center_co,rotmatrix(:,:,ids_co(j)), &
    j,npoints,volprot1,sumvolprot1,flag)

   npoints = 200
   call newintegrateg_cuboctahedron(Lcubell(j),Loctall(j),center_co,rotmatrix(:,:,ids_co(j)), j, &
    npoints,volx1,sumvolx1,com1,p1,ncha1,volxx1)

  !! eps
   voleps1 = voleps1-volprot1
   voleps1 = voleps1*eeps(ids_co(j))

   !area =  3.0**(1.0/2.0)*Loctall(j)**2 + 3.0*(1.0-3.0**(1.0/2.0))*(Loctall(j) - Lcubell(j))**2.0 ! MARIO

   if (Loctall(j).le.Lcubell(j)*2) then
    if (Loctall(j).le.Lcubell(j)) then
      area_co = 3.0**(1.0/2.0) *Loctall(j)**3 ! Octahedra
    else
      area_co = 3.0**(1.0/2.0)*Loctall(j)**2 + 3.0*(1.0 - 3.0**(1.0/2.0))*(Loctall(j) - Lcubell(j))**2.0
    endif

   elseif (Loctall(j).gt.Lcubell(j)*2) then
    if (Loctall(j).ge.Lcubell(j)*3) then
      area_co = 6*Lcubell(j)**2 ! Cube
    else
      if (Lcubell(j).gt.0) then
        area_co = 6*Lcubell(j)**2 + (1.0/32.0)*(2*3**(1.0/2.0) - 6)*Loctall(j)**2 *(3-Loctall(j)/Lcubell(j))**2
      else
        area_co = 0
      endif
    endif

   endif

  !! Normalize volx1 and volxx1 so that the integral is equal to the total number of ligands on the CO j

   volx1 = volx1/sumvolx1*area_co*sigma(ids_co(j))
   volxx1 = volxx1/sumvolx1*area_co*sigma(ids_co(j))

  !! Sum of number of polymer segments
   sumpolseg_co = sumpolseg_co + area_co*sigma(ids_co(j))*longp(ids_co(j))

  !! volume
   volprot1 = volprot1 * 0.99
   volprot = volprot+volprot1

  ! CHECK COLLISION HERE...
   if(maxval(volprot).gt.1.0) then ! collision
     flag=.true. 
   endif
   
   voleps = voleps + voleps1

  ! add com1 and volx to list

   volxx = volxx + volxx1 !actualizo volxx, verificar

  do i = 1, ncha1
  ncha = ncha+1
  volx(ncha) = volx1(i)
  com(ncha,:) = com1(i,:)
  p0(ncha,:) = p1(i,:)
  longc(ncha) = longp(ids_co(j))
  enddo

enddo


if(ncha.ge.maxvolx) then
write(stdout,*) 'geometries:', 'increase maxvolx'
call endall
endif

title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

temp = 0
do j = 1, NNNell
temp = temp + 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)
enddo
COvol = 0.0
do j = 1, NNNco
 if (Loctall(j).le.Lcubell(j)*2) then
  if (Loctall(j).le.Lcubell(j)) then
   COvol = COvol + (1.0/6.0)*Loctall(j)**3
  else
   COvol = COvol + (1.0/6.0)*Loctall(j)**3 - 0.5*(Loctall(j) -Lcubell(j))**3
  endif
 elseif (Loctall(j).gt.Lcubell(j)*2) then
   if (Loctall(j).ge.Lcubell(j)*3) then
    COvol = COvol + Lcubell(j)**3 - 4*((3.0/2.0)*Lcubell(j) - Loctall(j)/2)**3
   else
    COvol = COvol + Lcubell(j)**3
  endif
 endif
  
 if(COvol.le.0.0)COvol=Lcubell(j)**3

enddo

if (rank.eq.0) then
write(stdout,*) 'geometries:', 'update_matrix: Total volumen real space= ', temp+COvol
write(stdout,*) 'geometries:', 'update_matrix: Total discretized volumen =', sum(volprot)*delta**3
write(stdout,*) 'geometries:', 'number of monomers in system =', sumpolseg_ell + sumpolseg_co
endif

title = 'aveps'
counter = 1
call savetodisk(voleps, title, counter)

title = 'avgrf'
counter = 1
call savetodisk(volxx, title, counter)

end subroutine