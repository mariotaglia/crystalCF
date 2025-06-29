
subroutine initmpi
use MPI
use chainsdat
implicit none

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

end subroutine

subroutine initconst
use const
use molecules
use ellipsoid
use mparameters_monomer
implicit none
pi = acos(-1.0)

vsol = vsol0
error = 1e-4 ! para comparar con la norma...
errel=1d-6
itmax=200

eqs=(1+N_poorsol)
end subroutine

subroutine initall
use molecules
use const
use MPI
use ellipsoid
use chainsdat

use mparameters_monomer
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open common files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

if(rank.eq.0) then
       open(unit=301, file='F_tot_gcanon.dat', access='APPEND')
       open(unit=302, file='F_mixs.dat',  access='APPEND')
       open(unit=307, file='F_conf.dat',  access='APPEND')
       open(unit=3071, file='F_trans.dat',  access='APPEND')
       open(unit=309, file='F_vdW.dat',  access='APPEND')
       open(unit=311, file='F_HS.dat',  access='APPEND')
       open(unit=410, file='F_eps.dat',  access='APPEND')
       open(unit=312, file='F_tot2.dat',  access='APPEND')
       open(unit=314, file='F_mixpos2.dat',  access='APPEND')
       open(unit=315, file='F_trans_sv.dat',  access='APPEND')
       open(unit=316, file='F_conf_sv.dat',  access='APPEND')
       open(unit=420, file='F_tot_canon.dat',  access='APPEND')
endif
end subroutine

subroutine endall
use MPI
implicit none

!!!!!!!!!!!!!!!!!!!!!!
! Close common files
!!!!!!!!!!!!!!!!!!!!!!

close(301)
close(302)
close(307)
close(3071)
close(309)
close(310)
close(311)
close(312)
close(313)
close(315)
close(316)

call MPI_FINALIZE(ierr) ! finaliza MPI    
stop

end subroutine



subroutine savedata(cccc)
use system
use results
use const
use molecules
use chainsdat
use kai
use ematrix
use fields_fkfun
use MPI
use kinsol
use kaist
use mparameters_monomer
use channel
use solventchains
implicit none
integer cccc
character*20 filename
character*5  title
real*8 temp(dimx,dimy,dimz)
real*8 sumpol
integer ix,iy,iz, im
real*8 minpol
integer minpolpos(3)
!----------------------------------------------------------
!  OUTPUT
!----------------------------------------------------------

if(rank.eq.0) then ! solo el jefe escribe a disco....
  ! Guarda infile
!  write(filename,'(A4, I3.3, A4)')'out.', cccc, '.dat'
!  open(unit=45, file=filename)
!   do i = 1, 2*n
!    write(45, *)x1(i)
!   enddo
!  close(45)

!!!!!!!!!!!!!!!!!!! Guarda archivos !!!!!!!!!!!!!!!!!!!!!!!!!!!!


! minpol

  temp = 1.d50
  do ix = 1, dimx
  do iy = 1, dimy
  do iz = 1, dimz
  if(volprot(ix,iy,iz).eq.0.0) then
          
  temp(ix,iy,iz) = 0.0

  do im = 1, N_monomer
     temp(:,:,:) =  temp(:,:,:) + avpol(:,:,:, im)
  enddo

  endif

  enddo
  enddo
  enddo

  minpolpos = minloc(temp)
  minpol = minval(temp)

! Polimero, todo

  temp = 0.0
  do im = 1, N_monomer
     temp(:,:,:) =  temp(:,:,:) + avpol(:,:,:, im)*(1.0 - volprot(:,:,:))
  enddo

 title = 'avpol'
  call savetodisk(temp, title, cccc)

! Polimero, por tipo
  
  do im = 1, N_monomer
    temp(:,:,:) = avpol(:,:,:,im)*(1.0 - volprot(:,:,:))
    write(title,'(A3, I2.2)')'avp',im
    call savetodisk(temp, title, cccc)
  enddo

! Solvente
  temp(:,:,:) = xh(:,:,:)*(1.0 - volprot(:,:,:))

  title = 'avsol'
  call savetodisk(temp, title, cccc)

! save volprot for supercell
if(rank.eq.0) then
open (unit=8, file='out.par', form='unformatted')
do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
  xpar(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) = volprot(ix,iy,iz)
  enddo
 enddo
enddo
write(8)xpar
close(8)
endif

! system

  write(filename,'(A7, I4.4, A4)')'system.', cccc, '.dat'
  open (unit=310, file=filename)
  write(310,*)'st          = ',st
  write(310,*)'kp          = ',kp
  write(310,*)'eepsc       = ',eepsc ! Surface-polymer interaction
  write(310,*)'fnorm       = ',norma ! residual size of iteration vector
  write(310,*)'length seg  = ', lseg ! value see subroutine cadenas
  write(310,*)'length seg  = ', lseg ! value see subroutine cadenas
  write(310,*)'delta       = ',delta
  write(310,*)'vsol        = ',vsol
  write(310,*)'long        = ',long
  write(310,*)'iterations  = ',iter
  write(310,*)'sigma cad/nm2 = ',ncha/(dimx*dimy*delta*delta)
  write(310,*)'musolv =    ', musolv
  write(310,*)'phisolv =    ', phisolv
  write(310,*)'Nsolv =       ', phisolv*(dimx*dimy*dimz)*(delta**3)/vsol/longsv
  write(310,*)'GIT version = ', _VERSION
  write(310,*)'phisolv Oh =', xh(dimx/2+1,1,1)
  write(310,*)'phipol Oh =', avpol(dimx/2+1,1,1,1)
  write(310,*)'phisolv Th =', xh(dimx/2+1, 1, dimz/4+1)
  write(310,*)'phipol Th =', avpol(dimx/2+1, 1, dimz/4+1,1)
  write(310,*)'minpol =', minpol
  write(310,*)'minpolpos =', minpolpos
 

  sumpol = 0.0
  do ix = 1, dimx
  do iy = 1, dimy
  do iz = 1, dimz
  do im = 1, N_monomer
  sumpol = sumpol + avpol(ix,iy,iz,im)*(delta**3)*(1.0-volprot(ix,iy,iz))/vsol
  enddo
  enddo
  enddo
  enddo

  write(310,*)'Number of segments =          ', sumpol
  close(310)

endif

end subroutine

subroutine store2disk(counter) ! saves state to disk
use ellipsoid
use kinsol
use montecarlo
use ematrix
use results
use MPI
use const
implicit none
integer counter
character*20 filename

if(rank.eq.0) then
open (unit=8, file='out.out', form='unformatted')
write(8)counter
write(8)free_energy
write(8)xflag
close(8)
endif

if(rank.eq.0) then
write(filename,'(A4, I4.4, A4)')'out.', counter, '.dat'
open(unit=8, file=filename, form='unformatted')
write(8)counter
write(8)free_energy
write(8)xflag
close(8)
endif
end subroutine

subroutine retrivefromdisk(counter) ! saves state to disk
use ellipsoid
use kinsol
use montecarlo
use ematrix
use results
use const
implicit none
integer counter

open (unit=8, file='in.in', form='unformatted')
read(8)counter
read(8)free_energy
read(8)xflag
counter = 0
close(8)
end subroutine

subroutine mirror
use const
use kinsol
use mparameters_monomer
implicit none
real*8 xh(dimx,dimy,dimz)
real*8 xtotal(dimx,dimy,dimz,N_poorsol)
integer ip
integer ix,iy,iz
real*8 temp
integer ncells


ncells = dimx*dimy*dimz

do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
     xh(ix,iy,iz)=xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1))

     do ip = 1, N_poorsol
     xtotal(ix,iy,iz,ip) = xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells)
     enddo

  enddo
 enddo
enddo 
   
do ix=1,int(dimx/2)
 do iy=1,dimy
  do iz=1,dimz
      
     temp = xh(ix,iy,iz)
     xh(ix,iy,iz) = xh(dimx-ix,iy,iz)
     xh(dimx-ix,iy,iz) = temp

     do ip = 1, N_poorsol
       temp = xtotal(ix,iy,iz,ip)
       xtotal(ix,iy,iz,ip) = xtotal(dimx-ix,iy,iz,ip)
       xtotal(dimx-ix,iy,iz,ip) = temp
     enddo

  enddo
 enddo
enddo
    
  
do ix=1,dimx
   do iy=1,dimy
      do iz=1,dimz
        xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1))= xh(ix,iy,iz)

        do ip = 1, N_poorsol
         xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) =  xtotal(ix,iy,iz,ip)
        enddo

      enddo
   enddo
enddo


endsubroutine





