subroutine fkfun(x,f,ier2)

use system
use chainsdat
use molecules
use const
use results
use kai
use MPI
use fields_fkfun
use kinsol
use conformations
use ematrix
use ellipsoid
use transform
use kaist
use mparameters_monomer
use mmask
use solventchains
implicit none
real*8 intq, intxh
real*8 eta
integer*4 ier2
integer ncells
real*8 x(*),f(*)
real*8 protemp
integer i,j, ix, iy, iz, ii, ax, ay, az
integer im, ip
integer jx, jy, jz, jj
real*8 xpot(dimx, dimy, dimz, 0:N_monomer) ! 0 is solvent
real*8 xh_tosend(dimx,dimy,dimz)
real*8 qsv_tosend(dimx,dimy,dimz)
integer iii
integer, external :: PBCSYMI, PBCREFI

! poor solvent 
real*8 sttemp
! MPI
integer tag
parameter(tag = 0)
integer err
real*8 avpol_tosend(dimx,dimy,dimz,N_monomer)
real*8 avpol_temp(dimx,dimy,dimz,N_monomer)
real*8 q_tosend, sumtrans_tosend
real*8 fv, fv2

!-----------------------------------------------------
! Common variables

shift = 1.0d100

ncells = dimx*dimy*dimz ! numero de celdas

! Jefe

if(rank.eq.0) then ! llama a subordinados y pasa vector x
   flagsolver = 1
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
   CALL MPI_BCAST(x, eqs*ncells , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
endif

!------------------------------------------------------
! DEBUG
!      if(iter.gt.2000) then
!      do i = 1, n
!      write(stdout,*)i, x(i)
!      enddo
!      endif


! Recupera xh y psi desde x()
xtotalsum = 0.0


! ELECTRO
! psi = 0.0
do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
     xtotalsum(ix,iy,iz)= 1.0-exp(-x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)))  ! xtotalsum is the sum of polymers (ip >= 1) and solvent (ip = 0)

     do ip = 1, N_poorsol
      xtotal(ix,iy,iz,ip) = x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ ip*ncells) ! input, xtotal for polymers
     enddo

! ELECTRO
! if(electroflag.eq.1)psi(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells)   

  enddo
 enddo
enddo


! ELECTRO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! Boundary conditions electrostatic potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reflection or PBC, (PBC = 1 or 3)
! 
!do jx = 0, dimx+1
!do jy = 0, dimy+1
!do jz = 0, dimz+1
!
!ix=jx
!iy=jy
!iz=jz ! these lines are necessary for PBC = 0 or 2
!
!if (PBC(1).eq.1)ix = PBCSYMI(jx,dimx)
!if (PBC(3).eq.1)iy = PBCSYMI(jy,dimy)
!if (PBC(5).eq.1)iz = PBCSYMI(jz,dimz)
!
!if (PBC(1).eq.3)ix = PBCREFI(jx,dimx)
!if (PBC(3).eq.3)iy = PBCREFI(jy,dimy)
!if (PBC(5).eq.3)iz = PBCREFI(jz,dimz)
!
!   psi(jx, jy, jz) = psi(ix, iy, iz)
!enddo
!enddo
!enddo
!
!! Bulk or Wall, PBC = 0 or 2
!
!select case (PBC(1)) ! x = 0
!case(0) ! set bulk 
!   psi(0,:,:) = 0.0 
!case(2)
!   psi(0,:,:) = psi(1,:,:) ! zero charge
!endselect
!
!select case (PBC(2)) ! x = dimx
!case(0) ! set bulk 
!   psi(dimx+1,:,:) = 0.0  
!case(2)
!   psi(dimx+1,:,:) = psi(dimx,:,:) ! zero charge
!endselect
!
!select case (PBC(3)) ! y = 0
!case(0) ! set bulk 
!   psi(:,0,:) = 0.0  
!case(2)
!   psi(:,0,:) = psi(:,1,:) ! zero charge
!endselect
!
!select case (PBC(4)) ! y = dimy
!case(0) ! set bulk 
!   psi(:,dimy+1,:) = 0.0
!case(2)
!   psi(:,dimy+1,:) = psi(:,dimy,:) ! zero charge
!endselect
!
!select case (PBC(5)) ! z = 0
!case(0) ! set bulk 
!   psi(:,:,0) = 0.0  
!case(2)
!   psi(:,:,0) = psi(:,:,1) ! zero charge
!endselect
!
!select case (PBC(6)) ! z = dimz
!case(0) ! set bulk 
!   psi(:,:,dimz+1) = 0.0
!case(2)
!   psi(:,:,dimz+1) = psi(:,:,dimz) ! zero charge
!endselect
!
!! volume fraction and frdir



! ELECTRO
!
!fdis = 0.0
!
!
!do ix=1,dimx
! do iy=1,dimy
!  do iz=1,dimz
!    xpos(ix, iy, iz) = 0.0!expmupos*(xtotal(ix, iy, iz,0)**vsalt)*dexp(-psi(ix, iy, iz)*zpos) ! ion plus volume fraction 
!    xneg(ix, iy, iz) = 0.0!expmuneg*(xtotal(ix, iy, iz, 0)**vsalt)*dexp(-psi(ix, iy, iz)*zneg) ! ion neg volume fraction
!    xHplus(ix, iy, iz) = 0.0!expmuHplus*(xtotal(ix, iy, iz, 0))*dexp(-psi(ix, iy, iz))           ! H+ volume fraction
!    xOHmin(ix, iy,iz) = 0.0!expmuOHmin*(xtotal(ix,iy,iz, 0))*dexp(+psi(ix,iy,iz))           ! OH-  volume fraction
!
!     do im =1,N_monomer
!        if (zpol(im).eq.1) then !BASE
!          fdis(ix,iy,iz,im) = 0.0! 1.0 /(1.0 + xOHmin(ix,iy,iz)/(K0(im)*xtotal(ix,iy,iz,0)))
!        else if (zpol(im).eq.-1) then !ACID
!          fdis(ix,iy,iz,im) = 0.0! 1.0 /(1.0 + xHplus(ix,iy,iz)/(K0(im)*xtotal(ix,iy,iz,0)))
!        endif
!     enddo
!
!   enddo
! enddo  
!enddo
!


! ELECTRO
!xtotal(:,:,:,0)=xtotalsum(:,:,:) -xpos(:,:,:)-xneg(:,:,:)-xHplus(:,:,:)-xOHmin(:,:,:)
!xtotal(:,:,:,0)=1.0-xh(:,:,:)-xpos(:,:,:)-xneg(:,:,:)-xHplus(:,:,:)-xOHmin(:,:,:) 


! solvent from difference
xtotal(:,:,:,0)=xtotalsum(:,:,:)
do ip = 1, N_poorsol
  xtotal(:,:,:,0) = xtotal(:,:,:,0)-xtotal(:,:,:,ip) ! get solvent from difference
enddo

! Compute dielectric permitivity
! ELECTRO
!xtotalsum = 0.0 ! sum of all polymers
!do ip = 0, N_poorsol
!xtotalsum(:,:,:) = xtotalsum(:,:,:) + xtotal(:,:,:,ip)
!enddo
 
!call dielectfcn(xtotalsum,volprot,epsfcn,Depsfcn)

!------------------------------------------------------------------------
! PDFs polimero
!------------------------------------------------------------------------

! Calcula xpot

sttemp = st/vsol

do im = 0, N_monomer ! loop over different monomer types

do ix=1,dimx
 do iy=1,dimy
   do iz=1,dimz
     fv = (1.0 - volprot(ix,iy,iz))

! PACKING
!     xpot(ix, iy, iz, im) = xh(ix,iy,iz)**vpol


! LOCAL HS
     eta = xtotalsum(ix,iy,iz)
              xpot(ix, iy, iz, im) = (-(8.0*eta-(9.0*(eta**2))+(3.0*(eta**3))) &
              /((1.0-eta)**3))

! ELECTRO
!     if(zpol(im).ne.0.0) then
!         xpot(ix,iy,iz,im) =  xpot(ix,iy,iz,im)/fdis(ix,iy,iz,im)*dexp(-psi(ix,iy,iz)*zpol(im))
!     endif
 
! Dielectrics
!     gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz))**2 
!     gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))**2 
!     xpot(ix, iy, iz) = xpot(ix,iy,iz)*exp(-Depsfcn(ix,iy,iz)*(gradpsi2)*constqE)
!     xpot(ix,iy,iz,im) = xpot(ix,iy,iz,im)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0*vpol/fv)


! Poor solvent

     protemp=0.0

     do ax = -Xulimit,Xulimit 
      do ay = -Xulimit,Xulimit
       do az = -Xulimit,Xulimit

            jx = ix+ax
            jy = iy+ay
            jz = iz+az

            if(jx.lt.1) then
            if(PBC(1).eq.1)jx = PBCSYMI(jx,dimx)
            if(PBC(1).eq.3)jx = PBCREFI(jx,dimx)
            endif

            if(jx.gt.dimx) then
            if(PBC(2).eq.1)jx = PBCSYMI(jx,dimx)
            if(PBC(2).eq.3)jx = PBCREFI(jx,dimx)
            endif

            if(jy.lt.1) then
            if(PBC(3).eq.1)jy = PBCSYMI(jy,dimy)
            if(PBC(3).eq.3)jy = PBCREFI(jy,dimy)
            endif

            if(jy.gt.dimy) then
            if(PBC(4).eq.1)jy = PBCSYMI(jy,dimy)
            if(PBC(4).eq.3)jy = PBCREFI(jy,dimy)
            endif


            if(jz.lt.1) then
            if(PBC(5).eq.1)jz = PBCSYMI(jz,dimz)
            if(PBC(5).eq.3)jz = PBCREFI(jz,dimz)
            endif

            if(jz.gt.dimz) then
            if(PBC(6).eq.1)jz = PBCSYMI(jz,dimz)
            if(PBC(6).eq.3)jz = PBCREFI(jz,dimz)
            endif


            if((jx.ge.1).and.(jx.le.dimx)) then
            if((jy.ge.1).and.(jy.le.dimy)) then
            if((jz.ge.1).and.(jz.le.dimz)) then
                fv = (1.0-volprot(jx,jy,jz))

               do ip = 0, N_poorsol
               protemp=protemp + Xu(ax,ay,az)*st_matrix(hydroph(im),ip)*sttemp*xtotal(jx,jy,jz,ip)*fv
               enddo ! ip

            endif
            endif
            endif

       enddo
      enddo
     enddo

     xpot(ix,iy,iz,im) =xpot(ix,iy,iz,im) + protemp

   enddo ! ix
  enddo ! iy
enddo !iz
enddo ! N_monomer


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE SOLVENT VOLUME FRACTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


xh = 0.0
xh_tosend = 0.0
qsv = 0.0
qsv_tosend = 0.0
sumprolnpro = 0.0
rhosv = 0.0
sumprotrans = 0.0


do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz ! loop over position COM of solvent molecule

iii = ix+dimx*(iy-1)+dimx*dimy*(iz-1)     ! number of cell

if (mod(iii-1,size).eq.rank) then ! each processor runs on different cells
! rank 0 takes cell 1,1,1


do i = 1, cuantassv ! loop over sv conformations
prosv = -benergy*ntranssv(i) ! energy of trans bonds

do j = 1, longsv ! loop over segment

            jx = ix+pxsv(i,j)
            jy = iy+pysv(i,j)
            jz = iz+pzsv(i,j)

! CHECK PBC

            if(jx.lt.1) then
            if(PBC(1).eq.1)jx = PBCSYMI(jx,dimx)
            if(PBC(1).eq.3)jx = PBCREFI(jx,dimx)
            endif

            if(jx.gt.dimx) then
            if(PBC(2).eq.1)jx = PBCSYMI(jx,dimx)
            if(PBC(2).eq.3)jx = PBCREFI(jx,dimx)
            endif

            if(jy.lt.1) then
            if(PBC(3).eq.1)jy = PBCSYMI(jy,dimy)
            if(PBC(3).eq.3)jy = PBCREFI(jy,dimy)
            endif

            if(jy.gt.dimy) then
            if(PBC(4).eq.1)jy = PBCSYMI(jy,dimy)
            if(PBC(4).eq.3)jy = PBCREFI(jy,dimy)
            endif

            if(jz.lt.1) then
            if(PBC(5).eq.1)jz = PBCSYMI(jz,dimz)
            if(PBC(5).eq.3)jz = PBCREFI(jz,dimz)
            endif

            if(jz.gt.dimz) then
            if(PBC(6).eq.1)jz = PBCSYMI(jz,dimz)
            if(PBC(6).eq.3)jz = PBCREFI(jz,dimz)
            endif

            if((jx.ge.1).and.(jx.le.dimx)) then
            if((jy.ge.1).and.(jy.le.dimy)) then
            if((jz.ge.1).and.(jz.le.dimz)) then
 
            prosv = prosv+xpot(jx, jy, jz, 0)

            endif     
            endif     
            endif     
            
enddo ! j

   prosv = dexp(prosv)
   qsv_tosend(ix,iy,iz) = qsv_tosend(ix,iy,iz) + prosv
   sumprolnpro(ix,iy,iz) = sumprolnpro(ix,iy,iz) + prosv*dlog(prosv)
   sumprotrans(ix,iy,iz) = sumprotrans(ix,iy,iz) + prosv*ntranssv(i)

   fv = (1.0-volprot(ix,iy,iz))

do j=1,longsv ! calculate xhtemp

            jx = ix+pxsv(i,j)
            jy = iy+pysv(i,j)
            jz = iz+pzsv(i,j)

! CHECK PBC

            if(jx.lt.1) then
            if(PBC(1).eq.1)jx = PBCSYMI(jx,dimx)
            if(PBC(1).eq.3)jx = PBCREFI(jx,dimx)
            endif

            if(jx.gt.dimx) then
            if(PBC(2).eq.1)jx = PBCSYMI(jx,dimx)
            if(PBC(2).eq.3)jx = PBCREFI(jx,dimx)
            endif

            if(jy.lt.1) then
            if(PBC(3).eq.1)jy = PBCSYMI(jy,dimy)
            if(PBC(3).eq.3)jy = PBCREFI(jy,dimy)
            endif

            if(jy.gt.dimy) then
            if(PBC(4).eq.1)jy = PBCSYMI(jy,dimy)
            if(PBC(4).eq.3)jy = PBCREFI(jy,dimy)
            endif

            if(jz.lt.1) then
            if(PBC(5).eq.1)jz = PBCSYMI(jz,dimz)
            if(PBC(5).eq.3)jz = PBCREFI(jz,dimz)
            endif

            if(jz.gt.dimz) then
            if(PBC(6).eq.1)jz = PBCSYMI(jz,dimz)
            if(PBC(6).eq.3)jz = PBCREFI(jz,dimz)
            endif

            if((jx.ge.1).and.(jx.le.dimx)) then
            if((jy.ge.1).and.(jy.le.dimy)) then
            if((jz.ge.1).and.(jz.le.dimz)) then

              fv = (1.0-volprot(jx,jy,jz))
              fv2 = (1.0-volprot(ix,iy,iz))

              xh_tosend(jx,jy,jz) = xh_tosend(jx,jy,jz) + prosv*fv2/fv*vsol
            endif     
            endif     
            endif     
            
enddo ! j

enddo  ! i

endif ! processor

enddo ! ix
enddo ! iy 
enddo ! iz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE POLYMER VOLUME FRACTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

avpol = 0.0
avpol_tosend = 0.0
q = 0.0
sumtrans = 0.0

do jj = 1, cpp(rank+1)
   ii = cppini(rank+1)+jj

   q_tosend=0.0
   sumtrans_tosend = 0.0
   avpol_temp = 0.0

 do i=1,newcuantas(ii)
   pro(i, jj)=dlog(shift)
   do j=1,long
    ax = px(i, j, jj) ! cada uno para su cadena...
    ay = py(i, j, jj)
    az = pz(i, j, jj)         
    pro(i, jj) = pro(i, jj) + xpot(ax, ay, az, segtype(j))

   enddo
    pro(i,jj) = pro(i,jj) -benergy*ntrans(i,ii) ! energy of trans bonds
    pro(i,jj)=dexp(pro(i,jj))

   do j=1,long
   fv = (1.0-volprot(px(i,j, jj),py(i,j, jj),pz(i,j, jj)))
   im = segtype(j)
    avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj), im)= &
    avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj), im)+pro(i, jj)*vsol/(delta**3)/fv* &
    ngpol(ii)*sc ! ngpol(ii) has the number of chains grafted to the point ii
   enddo

   q_tosend=q_tosend+pro(i, jj)
   sumtrans_tosend = sumtrans_tosend+ntrans(i, ii)*pro(i,jj)

 enddo ! i
! norma 
    
avpol_tosend=avpol_tosend + avpol_temp/q_tosend

q(ii) = q_tosend ! no la envia ahora
sumtrans(ii) = sumtrans_tosend/q_tosend

!write(stdout,*) rank+1,jj,ii,q(ii)
enddo ! jj

!------------------ MPI ----------------------------------------------

call MPI_Barrier(MPI_COMM_WORLD, err)
call MPI_REDUCE(avpol_tosend, avpol, ncells*N_monomer, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
call MPI_REDUCE(xh_tosend, xh, ncells, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
call MPI_REDUCE(qsv_tosend, qsv, ncells, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

! Subordinados
if(rank.ne.0) then
!!!!!!!!!!! IMPORTANTE, LOS SUBORDINADOS TERMINAN ACA... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  goto 3333
endif


!!!!!!!!!!!!!!!!!!!!!!! Normalize solvent !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        intxh = 0.0
        do ix = 1, dimx
        do iy = 1, dimy
        do iz = 1, dimz
                fv = (1.0-volprot(ix,iy,iz))
                intxh  = intxh + xh(ix,iy,iz)*fv*delta**3
        enddo
        enddo
        enddo
! intx is the integral of xh before normalization

        intq = 0.0
        do ix = 1, dimx
        do iy = 1, dimy
        do iz = 1, dimz
                fv = (1.0-volprot(ix,iy,iz))
                intq = intq + qsv(ix,iy,iz)*fv*delta**3
        enddo
        enddo
        enddo
! intq = integral of qsv



if (flagmu.eq.0) then ! calculate using constant phi


!!!! Normalize xh
        xh = xh/intxh*kp*float(dimx*dimy*dimz)*delta**3
! musolv, see notes
        musolv = dlog(kp*float(dimx*dimy*dimz)*(delta**3)/float(longsv)/intq) 

else if (flagmu.eq.1) then  ! calculate using constant expmu

        musolv = kp
! xh, see notes
        xh = xh*exp(musolv)*intq*float(longsv)/intxh

else if (flagmu.eq.2) then ! calculate using constant Nsolv

        musolv = dlog(kp*vsol/intq) 

        !!!! Normalize xh
        xh = xh*exp(musolv)*intq*float(longsv)/intxh

endif

!!!! phisolv

phisolv = 0.0
do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz
fv = (1.0-volprot(ix,iy,iz))
phisolv = phisolv + xh(ix,iy,iz)*fv
enddo
enddo
enddo
phisolv = phisolv/float(dimx*dimy*dimz) 

!! CHECK MUSOLV

!do ix = 1, dimx
!do iy = 1, dimy
!do iz = 1, dimz

!print*, musolv, dlog(rhosv(ix,iy,iz)*vsol/qsv(ix,iy,iz))

!enddo
!enddo
!enddo
!stop

!! CHECK AVERAGE SOLV DENSITY

!if(rank.eq.0)write(stdout,*)'Target kp', kp

! FROM XH
!temp = 0.0
!do ix = 1,dimx
!do iy = 1,dimy
!do iz = 1,dimz
!  fv = (1.0-volprot(ix,iy,iz))
!  temp = temp + fv*xh(ix,iy,iz)
!enddo
!enddo
!enddo
!temp = temp/float(dimx*dimy*dimz) ! average xh 

!if(rank.eq.0)write(stdout,*)'kp from xh', temp

!! FROM RHOSV
!temp = 0.0
!do ix = 1,dimx
!do iy = 1,dimy
!do iz = 1,dimz
!  fv = (1.0-volprot(ix,iy,iz))
!  temp = temp + fv*rhosv(ix,iy,iz)
!enddo
!enddo
!enddo
!temp = temp*float(longsv)*vsol/float(dimx*dimy*dimz) ! average xh 

!if(rank.eq.0)write(stdout,*)'kp from rhosv', temp

!stop


!!!!!!!!!!!!!!!!!!!!!!! FIN MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------------------------
!   Construye Ecuaciones a resolver 
!----------------------------------------------------------------------------------------------

! Qtot



! ELECTRO
!qtot = 0.0
!do ix=1,dimx
!do iy=1,dimy
!do iz=1,dimz
!  
! fv = (1.0-volprot(ix,iy,iz))
!
! qtot(ix, iy, iz) =  (zpos*xpos(ix, iy, iz)+zneg*xneg(ix, iy, iz))/vsalt + xHplus(ix, iy, iz) - xOHmin(ix, iy, iz)
!
! do im = 1, N_monomer
!     qtot(ix, iy, iz) =  qtot(ix,iy,iz) + avpol(ix,iy,iz,im)*zpol(im)/vpol*fdis(ix,iy,iz,im)
! enddo
!
! qtot(ix, iy,iz) = qtot(ix,iy,iz)*fv + volq(ix,iy,iz)*vsol    ! OJO
!
!enddo
!enddo
!enddo
!
! Volume fraction

do ix=1,dimx
do iy=1,dimy
do iz=1,dimz

! ELECTRO
!f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))= xh(ix,iy,iz) + &
!      xneg(ix, iy, iz) + xpos(ix, iy, iz) + xHplus(ix, iy, iz) + &
!      xOHmin(ix, iy, iz) -1.000000d0

f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))= -xtotalsum(ix,iy,iz)+xh(ix,iy,iz) ! xtotalsum = solvent + polymers
do im = 1, N_monomer
    f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))= f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) + avpol(ix,iy,iz,im)
enddo ! im


enddo
enddo
enddo

! Poor solvent

do ix=1,dimx
do iy=1,dimy
do iz=1,dimz

do ip = 1, N_poorsol
  f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) = xtotal(ix,iy,iz,ip)

  do im = 1, N_monomer
   if(hydroph(im).eq.ip) then 
    f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) = f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) - avpol(ix,iy,iz,im)
   endif
  enddo ! im

enddo ! ip
enddo ! ix
enddo ! iy
enddo ! iz


! ELECTRO
!
!if(electroflag.eq.1) then
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Poisson equatio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!
!! Some auxialiary variables, see Notes Poisson eq. non-cubic grid
!!
!
!MV(1) = MAT(1,1)
!MV(2) = MAT(1,2)  
!MV(3) = MAT(1,3)
!
!MU(1) = MAT(2,1)
!MU(2) = MAT(2,2)  
!MU(3) = MAT(2,3)
!
!MW(1) = MAT(3,1)
!MW(2) = MAT(3,2)  
!MW(3) = MAT(3,3)
!
!MVV = DOT_PRODUCT(MV,MV)
!MUU = DOT_PRODUCT(MU,MU)
!MWW = DOT_PRODUCT(MW,MW)
!
!MVU = DOT_PRODUCT(MV,MU)
!MVW = DOT_PRODUCT(MV,MW)
!MUW = DOT_PRODUCT(MU,MW)
!
!do ix=1,dimx
!do iy=1,dimy
!do iz=1,dimz
!
!psivv = psi(ix+1,iy,iz)-2*psi(ix,iy,iz)+psi(ix-1,iy,iz)
!psiuu = psi(ix,iy+1,iz)-2*psi(ix,iy,iz)+psi(ix,iy-1,iz)
!psiww = psi(ix,iy,iz+1)-2*psi(ix,iy,iz)+psi(ix,iy,iz-1)
!
!psivu = (psi(ix+1,iy+1,iz)+psi(ix-1,iy-1,iz)-psi(ix+1,iy-1,iz)-psi(ix-1,iy+1,iz))/4.0
!psivw = (psi(ix+1,iy,iz+1)+psi(ix-1,iy,iz-1)-psi(ix+1,iy,iz-1)-psi(ix-1,iy,iz+1))/4.0
!psiuw = (psi(ix,iy+1,iz+1)+psi(ix,iy-1,iz-1)-psi(ix,iy+1,iz-1)-psi(ix,iy-1,iz+1))/4.0
!
!psiv(1) = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))/2.0
!psiv(2) = (psi(ix,iy+1,iz)-psi(ix,iy-1,iz))/2.0
!psiv(3) = (psi(ix,iy,iz+1)-psi(ix,iy,iz-1))/2.0
!
!epsv(1) = (epsfcn(ix+1,iy,iz)-epsfcn(ix-1,iy,iz))/2.0
!epsv(2) = (epsfcn(ix,iy+1,iz)-epsfcn(ix,iy-1,iz))/2.0
!epsv(3) = (epsfcn(ix,iy,iz+1)-epsfcn(ix,iy,iz-1))/2.0
!
!psitemp = epsfcn(ix,iy,iz)*(MVV*psivv+MUU*psiuu+MWW*psiww+2.0*MVU*psivu+2.0*MVW*psivw+2.0*MUW*psiuw)
!psitemp = psitemp + DOT_PRODUCT(MATMUL(TMAT,epsv),MATMUL(TMAT,psiv))
!
!! OJO CHECK!!!!
!
!f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells)=(psitemp + qtot(ix, iy, iz)*constq)/(-2.0)
!
!
!enddo
!enddo
!enddo
!
!endif ! electroflag
! 
norma = 0.0

do i = 1, eqs*ncells
 norma = norma + (f(i))**2
enddo

iter = iter + 1
if(verbose.ge.3) then
if(rank.eq.0)write(stdout,*)'fkfun:', iter, norma, q(1)
endif

if(isnan(norma)) then
    if(rank.eq.0)write(stdout,*)'Norma is NaN, stop'
    f(1:eqs*ncells) = 0.0
endif
if(iter.gt.maxiters) then
    if(rank.eq.0)write(stdout,*)'Iter > Maxiters, stop'
    f(1:eqs*ncells) = 0.0
endif


       
     

3333 continue
ier2 = 0.0 

return
end
