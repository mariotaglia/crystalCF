!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!    Free Energy Calculation...
!
!
!
subroutine Free_Energy_Calc(looped)

use system
use const
use fields_fkfun
use MPI
use molecules
use kai
use results
use ematrix
use montecarlo
use ellipsoid
use transform
use kaist
use conformations
use mparameters_monomer
use solventchains
implicit none

integer looped
real*8  q_tosend(ncha), sumtrans_tosend(ncha)
real*8  q0(ncha), sumtrans0(ncha)
integer newcuantas0(ncha)
real*8 F_Mix_s
real*8 Free_energy2, sumrho, suma, mupol, sumHS
real*8 temp
real*8 F_trans, F_Conf, F_vdW, F_eps, F_HS
real*8 F_conf_sv, F_trans_sv
real*8 Free_Energy_plusSv
real*8 pro0(cuantas, maxcpp)
real*8 entropy(dimx,dimy,dimz)
character*5  title
real*8 eta

! MPI
integer stat(MPI_STATUS_SIZE) 
integer source
integer dest
integer tag
parameter(tag = 0)
integer err

! Dummies
integer ix, iy, iz, i, ii, ax, ay, az, jj
integer jx, jy, jz,iii
integer im, ip, ipp
real*8 fv, fv2

integer, external :: PBCSYMI, PBCREFI

! Solvent data for F_conf
real*8 sumprolnpro0(dimx,dimy,dimz)
real*8 sumprotrans0(dimx,dimy,dimz)


!!!!!!!!!!!!!!!!!!!!!!!!! MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Polymer
!

entropy = 0.0
q0 = 0.0
q_tosend = 0.0
sumtrans_tosend = 0.0

if(rank.ne.0) then
       dest = 0
! Envia q

       do jj = 1, cpp(rank+1)
       iii = cppini(rank+1)+jj
       q_tosend(iii) = q(iii)
       enddo

        call MPI_REDUCE(q_tosend, q0, ncha, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

! newcuantas
        call MPI_REDUCE(newcuantas, newcuantas0, ncha, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)

! Envia pro
        CALL MPI_SEND(pro, cuantas*cpp(rank+1) , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD,err)

! sum trans
       do jj = 1, cpp(rank+1)
       iii = cppini(rank+1)+jj
       sumtrans_tosend(iii) = sumtrans(iii)
       enddo
        call MPI_REDUCE(sumtrans_tosend, sumtrans0, ncha, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)


!!!!! Solvent

      call MPI_REDUCE(sumprolnpro, sumprolnpro0, dimx*dimy*dimz, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
      call MPI_REDUCE(sumprotrans, sumprotrans0, dimx*dimy*dimz, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

      goto 888 !!!! PROCESSORS NOT EQ 1 GO TO THE END 

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Free_Energy = 0.0
Free_Energy2 = 0.0




!! ELECTRO      
!! 2. Mezcla ion positivo
!
!      F_Mix_pos = 0.0 
!
!      do ix = 1, dimx
!      do iy = 1, dimy
!      do iz = 1, dimz
!      
!      fv=(1.0-volprot(ix,iy,iz))
!
!      F_Mix_pos = F_Mix_pos + xpos(ix, iy,iz) &
!      *(dlog(xpos(ix, iy, iz)/vsalt)-1.0-dlog(expmupos) + dlog(vsalt))*fv
!
!      F_Mix_pos = F_Mix_pos - xposbulk &
!      *(dlog(xposbulk/vsalt)-1.0-dlog(expmupos) + dlog(vsalt))*fv
!
!      enddo
!      enddo
!      enddo
!      F_Mix_pos = F_Mix_pos * delta**3/vsol/vsalt
!      Free_Energy = Free_Energy + F_Mix_pos
!
!! 3. Mezcla ion negativo
!
!      F_Mix_neg = 0.0
!
!      do ix = 1, dimx
!      do iy = 1, dimy
!      do iz = 1, dimz
!
!      fv=(1.0-volprot(ix,iy,iz))
!
!      F_Mix_neg = F_Mix_neg + xneg(ix, iy,iz) &
!      *(dlog(xneg(ix, iy, iz)/vsalt)-1.0- dlog(expmuneg) + dlog(vsalt))*fv
!
!      F_Mix_neg = F_Mix_neg - xnegbulk &
!      *(dlog(xnegbulk/vsalt)-1.0- dlog(expmuneg) + dlog(vsalt))*fv
!
!      enddo 
!      enddo 
!      enddo 
!      F_Mix_neg = F_Mix_neg * delta**3/vsol/vsalt
!      Free_Energy = Free_Energy + F_Mix_neg
!
!! 4. Mezcla protones
!
!      F_Mix_Hplus = 0.0
!
!      do ix = 1, dimx
!      do iy = 1, dimy
!      do iz = 1, dimz
!
!      fv=(1.0-volprot(ix,iy,iz))
!
!      F_Mix_Hplus = F_Mix_Hplus &
!     +xHplus(ix, iy, iz)*(dlog(xHplus(ix,iy,iz))-1.0 -dlog(expmuHplus))*fv
!
!      F_Mix_Hplus = F_Mix_Hplus &
!     -xHplusbulk*(dlog(xHplusbulk)-1.0 -dlog(expmuHplus))*fv
!
!      enddo
!      enddo
!      enddo
!      F_Mix_Hplus = F_Mix_Hplus * delta**3/vsol
!      Free_Energy = Free_Energy + F_Mix_Hplus
!
!! 5. Mezcla hidroxilos
!
!      F_Mix_OHmin = 0.0
!
!      do ix = 1, dimx
!      do iy = 1, dimy
!      do iz = 1, dimz
!
!      fv=(1.0-volprot(ix,iy,iz))
!
!      F_Mix_OHmin = F_Mix_OHmin + xOHmin(ix, iy,iz)*(dlog(xOHmin(ix, iy, iz))-1.0-dlog(expmuOHmin))*fv
!
!      F_Mix_OHmin = F_Mix_OHmin - xOHminbulk*(dlog(xOHminbulk)-1.0-dlog(expmuOHmin))*fv
!
!      enddo
!      enddo
!      enddo
!      F_Mix_OHmin = F_Mix_OHmin * delta**3/vsol
!      Free_Energy = Free_Energy + F_Mix_OHmin
!
! 6. Entropia interna polimero

      F_Conf = 0.0

! Jefe

       if (rank.eq.0) then ! Igual tiene que serlo, ver arriba

       do jj = 1, cpp(rank+1)
       iii = jj
       q_tosend(iii) = q(iii)
       enddo

        call MPI_REDUCE(q_tosend, q0, ncha, &
        MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

        call MPI_REDUCE(newcuantas, newcuantas0, ncha, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)

       do jj = 1, cpp(rank+1)
       do i = 1, newcuantas0(jj)
       iii = jj
      
         F_Conf = F_Conf + (pro(i, jj)/q0(iii)) &
      *dlog((pro(i, jj))/q0(iii))*ngpol(iii)

       entropy(p0(iii,1),p0(iii,2),p0(iii,3)) =  - dlog(q0(iii)/shift) 
       enddo
       enddo 

         do ii = 2, size ! loop sobre los procesadores restantes

        source = ii-1

        call MPI_RECV(pro0, cuantas*cpp(ii), &
        MPI_DOUBLE_PRECISION, source, tag, MPI_COMM_WORLD,stat, err)


       do jj = 1, cpp(ii)
!       write(stdout,*) ii, jj, pro0(10,jj)
       iii = cppini(ii)+jj
       do i = 1, newcuantas0(iii)

         F_Conf = F_Conf + (pro0(i, jj)/q0(iii))*dlog((pro0(i, jj))/q0(iii))*ngpol(iii)

       entropy(p0(iii,1),p0(iii,2),p0(iii,3)) =  - dlog(q0(iii)/shift) 

       enddo
       enddo

         enddo ! ii

       endif ! rank

      Free_Energy = Free_Energy + F_Conf

if(rank.eq.0) then

!      title = 'entpy'
!      call savetodisk(entropy, title, looped)
 
      open (unit=8, file='entropy.out', form='unformatted')
      write(8)dimx,dimy,dimz
      write(8)entropy
      close(8)
endif


! 6.5 Energy of trans bonds

      F_trans = 0.0

! Jefe

if (rank.eq.0) then ! Igual tiene que serlo, ver arriba

       do jj = 1, cpp(rank+1) ! sumtrans in rank 0
       iii = jj
       sumtrans_tosend(iii) = sumtrans(iii)
       enddo

        call MPI_REDUCE(sumtrans_tosend, sumtrans0, ncha, &
        MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

       do ii = 1, ncha
       F_trans = F_trans + sumtrans0(ii)*ngpol(ii)*benergy
       enddo  

       endif ! rank

      Free_Energy = Free_Energy + F_trans

if(rank.eq.0) then
      title = 'entpy'
      call savetodisk(entropy, title, looped)
 
      open (unit=8, file='entropy.out', form='unformatted')
      write(8)dimx,dimy,dimz
      write(8)entropy
      close(8)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solvent (processor #0)


sumprolnpro0 = 0.0
sumtrans0 = 0.0

if (rank.eq.0) then 
call MPI_REDUCE(sumprolnpro, sumprolnpro0, dimx*dimy*dimz, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
call MPI_REDUCE(sumprotrans, sumprotrans0, dimx*dimy*dimz, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif

rhosv = qsv*dexp(musolv)/vsol 

!!!!!!!!!!!!!!!!!!! ONLY CALCULATE SOLVENT CONTRIBUTIONS IF KP != 0 !!!!!!!!!!!!!!!!
if((kp.ne.0.0).or.(flagmu.eq.1)) then 

! 1. Mezcla solvente

F_Mix_s = 0.0 

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz
   fv=(1.0-volprot(ix,iy,iz))
   F_Mix_s = F_Mix_s + rhosv(ix,iy,iz)*(dlog(rhosv(ix, iy, iz)*vsol)-1.0-musolv)*fv
!      F_Mix_s = F_Mix_s - xsolbulk*(dlog(xsolbulk)-1.0)*fv
enddo      
enddo      
enddo      

F_Mix_s = F_Mix_s * delta**3
Free_Energy = Free_Energy + F_Mix_s
      

! 6-bis

! Conformational entropy of the solvent

F_conf_sv = 0.0

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz

fv=(1.0-volprot(ix,iy,iz))
F_conf_sv = F_conf_sv + (sumprolnpro0(ix,iy,iz)/qsv(ix,iy,iz) - dlog(qsv(ix,iy,iz)))*rhosv(ix,iy,iz)*fv


enddo
enddo
enddo

F_conf_sv = F_conf_sv*(delta**3)
Free_Energy = Free_Energy + F_conf_sv 

! 6.5 bis

! Gauche energy solvent

F_trans_sv = 0.0

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz

fv=(1.0-volprot(ix,iy,iz))
F_trans_sv = F_trans_sv + sumprotrans0(ix,iy,iz)/qsv(ix,iy,iz)*rhosv(ix,iy,iz)*fv*benergy

enddo
enddo
enddo


F_trans_sv = F_trans_sv*(delta**3)
Free_Energy = Free_Energy + F_trans_sv 

endif ! solvent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!      
!! ELECTRO
!! 7. Chemical Equilibrium
!      F_Eq = 0.0 
!            
!
!      do ix  = 1, dimx
!      do iy  = 1, dimy
!      do iz  = 1, dimz
!
!      do im = 1, N_monomer
!      
!      fv=(1.0-volprot(ix,iy,iz))
!
!      if(zpol(im).ne.0) then
!
!      F_Eq = F_Eq + fdis(ix,iy,iz,im)*dlog(fdis(ix,iy,iz,im)) &
!      *avpol(ix,iy,iz,im)/vpol*fv
!
!      F_Eq = F_Eq + (1.0-fdis(ix,iy,iz,im)) &
!      *dlog(1.0-fdis(ix,iy,iz,im))*avpol(ix,iy,iz,im)/vpol*fv
!
!      F_Eq = F_Eq + (1.0-fdis(ix,iy,iz,im))*dlog(K0(im))*avpol(ix,iy,iz,im)/vpol*fv
!
!      select case (zpol(im))
!      case (-1) ! acid
!       F_Eq = F_Eq + (1.0-fdis(ix,iy,iz,im))*(-dlog(expmuHplus))*avpol(ix,iy,iz,im)/vpol*fv
!      case (1) ! base
!       F_Eq = F_Eq + (1.0-fdis(ix,iy,iz,im))*(-dlog(expmuOHmin))*avpol(ix,iy,iz,im)/vpol*fv
!      end select
!
!      endif ! zpol
!
!      enddo ! im
!   
!      enddo
!      enddo
!      enddo
!
!      F_eq = F_eq *delta**3/vsol
!
!      Free_Energy = Free_Energy + F_Eq


! 8.vdW ! Ojo, los kai son negativos => atraccion

       F_vdW = 0.0

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz

      fv=(1.0-volprot(ix,iy,iz))

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
                fv2 = (1.0-volprot(jx,jy,jz)) 

            do ip = 0, N_poorsol
            do ipp = 0, N_poorsol
 
                F_vdW = F_vdW - 0.5000*delta**3*xtotal(ix,iy,iz,ip) &
        *xtotal(jx,jy,jz,ipp)*Xu(ax, ay, az)*st*st_matrix(ip,ipp)*fv*fv2/(vsol*vsol)
 
            enddo ! ip
            enddo ! ipp

            endif
            endif
            endif

            enddo
            enddo
            enddo

      enddo
      enddo
      enddo

      Free_Energy = Free_Energy + F_vdW

!! ELECTRO      
!! 9. Electrostatic ! OJO
!
!      F_electro = 0.0    
!
!      do ix  = 1, dimx
!      do iy  = 1, dimy
!      do iz  = 1, dimz
!
!      F_electro = F_electro &
!       + delta**3*psi(ix, iy, iz)*qtot(ix, iy, iz)/2.0/vsol
!
!      enddo
!      enddo
!      enddo
!  
!      print*, F_electro
!
!      Free_Energy = Free_Energy + F_electro
!
! 10. Pol-prot

      F_eps = 0.0 

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz
      fv=(1.0-volprot(ix,iy,iz))
      do im = 1, N_monomer
      F_eps = F_eps - avpol(ix,iy,iz,im)*voleps(ix,iy,iz)*(delta**3)/vsol*fv
      enddo
      enddo
      enddo
      enddo

      Free_Energy = Free_Energy + F_eps

! 10 HS contribution

      F_HS = 0.0
      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz
         fv=(1.0-volprot(ix,iy,iz))
         eta = xtotalsum(ix,iy,iz)
         F_HS = F_HS + eta*(4.-3.*eta)/((1.-eta)**2)*eta*(delta**3)/vsol*fv
      enddo
      enddo
      enddo

      Free_Energy = Free_Energy + F_HS

      write(stdout,*) 'Free_Energy_Calc: Free energy(1) = ', Free_energy

! Calculation of Omega + muSolv*NSolv

Free_Energy_plusSV = Free_Energy

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz

fv=(1.0-volprot(ix,iy,iz))
Free_Energy_plusSV = Free_Energy_plusSV + musolv*rhosv(ix,iy,iz)*fv*delta**3

enddo
enddo
enddo

! minimal F

      Free_Energy2 = 0.0

!      sumpi = 0.0
      sumrho=0.0
      sumHS = 0.0 

! ELECTRO
!        sumel=0.0
!        sumdiel = 0.0

        do ix=1,dimx
        do iy=1,dimy
        do iz=1,dimz

        fv=(1.0-volprot(ix,iy,iz))
        sumrho = sumrho - rhosv(ix, iy, iz)*fv

!           sumpi = sumpi+dlog(xh(ix, iy, iz))*fv     
!           sumpi = sumpi-dlog(xsolbulk)*fv

! ELECTRO 
!           sumrho = sumrho + ( - xh(ix, iy, iz) -xHplus(ix, iy, iz) &
!        - xOHmin(ix, iy, iz) - (xpos(ix, iy, iz)+xneg(ix, iy, iz))/vsalt)*fv! sum over  rho_i i=+,-,s

!           sumrho = sumrho - ( - xsolbulk -xHplusbulk &
!       -xOHminbulk - (xposbulk+xnegbulk)/vsalt)*fv ! sum over  rho_i i=+,-,s
     
!          sumrho = sumrho + ( - xh(ix, iy, iz))*fv! sum over  rho_i i=+,-,s
!           sumrho = sumrho - ( - xsolbulk)*fv ! sum over  rho_i i=+,-,s

! ELECTRO           
!         sumel = sumel - qtot(ix, iy, iz)*psi(ix, iy, iz)/2.0     
!         sumel = sumel + volq(ix,iy,iz)*psi(ix,iy,iz)*vsol                   


!         psiv(1) = psi(ix+1,iy,iz)-psi(ix,iy,iz)
!         psiv(2) = psi(ix,iy+1,iz)-psi(ix,iy,iz)
!         psiv(3) = psi(ix,iy,iz+1)-psi(ix,iy,iz)
!
!         gradpsi2 = DOT_PRODUCT(MATMUL(TMAT, psiv), MATMUL(TMAT, psiv))
!         sumdiel = sumdiel + 0.5/constq*xtotalsum(ix,iy,iz)*gradpsi2*Depsfcn(ix,iy,iz)
!
         eta = xtotalsum(ix,iy,iz)

         sumHS = sumHS + (-4.*eta + 2.*eta**2)/((1.-eta)**3) * fv * eta 

         enddo
         enddo
         enddo
         
 !        sumpi = (delta**3/vsol)*sumpi
         sumrho = (delta**3)*sumrho

!         print*, 'sumrho', sumrho


         sumHS = (delta**3/vsol)*sumHS

! ELECTRO
!         sumel = (delta**3/vsol)*sumel
!         sumdiel = (delta**3/vsol)*sumdiel
!         suma = sumpi + sumrho + sumel + sumdiel


!         suma = sumpi + sumrho
         suma = sumHS + sumrho

         do ii = 1, ncha
         Free_Energy2 = Free_Energy2-dlog(q0(ii)/shift)*ngpol(ii) 
         enddo

         Free_Energy2 = Free_Energy2 + suma - F_vdW

      write(stdout,*) 'Free_Energy_Calc: Free energy(2) = ', Free_energy2

! Guarda energia libre


        mupol = 0.0
        do ii = 1, ncha
        mupol = mupol - dlog(q0(ii)/shift)*ngpol(ii)
        enddo

        temp = sum(ngpol)
        mupol = mupol/temp

          if(rank.eq.0) then

         write(301,*)looped, Free_energy
          flush(301)
         write(302,*)looped, F_Mix_s 

! ELECTRO
!         write(303,*)looped, F_Mix_pos
!         write(304,*)looped, F_Mix_neg
!         write(305,*)looped, F_Mix_Hplus
!         write(306,*)looped, F_Mix_OHmin
!         write(308,*)looped, F_Eq
!         write(311,*)looped, F_electro


         write(3071,*)looped, F_trans
         write(315,*)looped, F_trans_sv
         write(316,*)looped, F_conf_sv


         write(307,*)looped, F_Conf
         write(309,*)looped, F_vdW
         write(311,*)looped, F_HS
         write(410,*)looped, F_eps

         write(420,*)looped, Free_Energy_plusSV
         flush(420)
 
         write(312,*)looped, Free_energy2

         write(313,*)looped, mupol

         endif
 
 888     call MPI_BCAST(free_energy, 1, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, err)

         return

         end




