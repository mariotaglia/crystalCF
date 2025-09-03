subroutine fkfun(x,f,ier2)

use chainsdat, only : longcha
use molecules, only : benergy, vsol
use const, only : stdout
use results, only : xtotalsum, avpol
use kai, only : Xu, Xulimit
use MPI
use fields_fkfun, only : xtotal, sumprolnpro, sumprotrans, long, phisolv, musolv, &
    pro, prosv, newcuantas, ngpol, cpp, cppini, segtype, xh, shift, sumtrans, &
    q, qsv, rhosv
use kinsol, only : maxiters, iter, norma
use conformations, only : px,py,pz, ntrans
use ematrix, only : dimx, dimy, dimz, eqs, volprot, pbc, delta, flagmu
use kaist, only : kp, st
use mparameters_monomer, only : N_monomer, N_poorsol, hydroph, st_matrix
use solventchains, only : pxsv, pysv, pzsv, ntranssv, longsv, cuantassv
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


do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
     xtotalsum(ix,iy,iz)= 1.0-exp(-x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)))  ! xtotalsum is the sum of polymers (ip >= 1) and solvent (ip = 0)

     do ip = 1, N_poorsol
      xtotal(ix,iy,iz,ip) = x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ ip*ncells) ! input, xtotal for polymers
     enddo

  enddo
 enddo
enddo

! solvent from difference
xtotal(:,:,:,0)=xtotalsum(:,:,:)
do ip = 1, N_poorsol
  xtotal(:,:,:,0) = xtotal(:,:,:,0)-xtotal(:,:,:,ip) ! get solvent from difference
enddo

! Calcula xpot

sttemp = st/vsol

do im = 0, N_monomer ! loop over different monomer types

do ix=1,dimx
 do iy=1,dimy
   do iz=1,dimz
     fv = (1.0 - volprot(ix,iy,iz))

! LOCAL HS
     eta = xtotalsum(ix,iy,iz)
              xpot(ix, iy, iz, im) = (-(8.0*eta-(9.0*(eta**2))+(3.0*(eta**3))) &
              /((1.0-eta)**3))

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
   do j=1,longcha(ii)
    ax = px(i, j, jj) ! cada uno para su cadena...
    ay = py(i, j, jj)
    az = pz(i, j, jj)         
    pro(i, jj) = pro(i, jj) + xpot(ax, ay, az, segtype(j))

   enddo
    pro(i,jj) = pro(i,jj) -benergy*ntrans(i,ii) ! energy of trans bonds
    pro(i,jj)=dexp(pro(i,jj))

   do j=1,longcha(ii)
   fv = (1.0-volprot(px(i,j, jj),py(i,j, jj),pz(i,j, jj)))
   im = segtype(j)
    avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj), im)= &
    avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj), im)+pro(i, jj)*vsol/(delta**3)/fv* &
    ngpol(ii) ! ngpol(ii) has the number of chains grafted to the point ii
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

! Volume fraction

do ix=1,dimx
do iy=1,dimy
do iz=1,dimz

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

norma = 0.0

do i = 1, eqs*ncells
 norma = norma + (f(i))**2
enddo

iter = iter + 1
if(rank.eq.0)write(stdout,*)'fkfun:', iter, norma, q(1)

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
