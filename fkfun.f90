subroutine fkfun(x,f,ier2)
use system
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

type :: IndexMap
        integer :: jj       ! 1ra prioridad: Punto de injerto
        integer :: ii       ! 2da prioridad: Cadena
        integer :: i_conf   ! 3ra prioridad: Configuración
        integer(kind=8) :: pos ! Posición física en el archivo 90[cite: 4, 6]
end type IndexMap

type(IndexMap), allocatable, save :: sorted_idx(:)
integer, save :: total_records = 0
integer :: k
integer(kind=8) :: file_size_bytes


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
integer :: id_cha, l_cha, ntrans_val
integer :: jj_read, i_read
real*8  :: pro_val
integer(kind=8) :: pos_read

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

if (flag_write_pxyz == 1) then 
    ! Solo entramos aquí si el arreglo NO ha sido asignado aún
    if (.not. allocated(sorted_idx)) then
        rewind(90)
        rewind(91)

        inquire(unit=91, size=file_size_bytes)
        total_records = int(file_size_bytes / 20) 
        
        allocate(sorted_idx(total_records))

        do k = 1, total_records
            read(91) sorted_idx(k)%jj, sorted_idx(k)%ii, sorted_idx(k)%i_conf, sorted_idx(k)%pos
        end do

        call quicksort_idx(sorted_idx, 1, total_records)
    endif
endif

do jj = 1, cpp(rank+1)
    ii = cppini(rank+1) + jj

    q_tosend = 0.0
    sumtrans_tosend = 0.0
    avpol_temp = 0.0

    do i = 1, newcuantas(ii)
        pro(i, jj) = dlog(shift)
        
        if (flag_write_pxyz.eq.1) then
            pos_read = find_pos_in_index(jj, ii, i)
            read(90, pos=pos_read) jj_read, id_cha, i_read, ntrans_val, l_cha, &
            px(1,1:l_cha,jj_read), & 
            py(1,1:l_cha,jj_read), & 
            pz(1,1:l_cha,jj_read)   
            if (id_cha /= ii) then 
                write(stdout,*)'MISMATCH id_cha', jj_read, id_cha, i_read, ii, jj, ii
                stop
            endif  
            do j = 1, l_cha
                ax = px(1, j, jj)
                ay = py(1, j, jj)
                az = pz(1, j, jj)         
                pro(i, jj) = pro(i, jj) + xpot(ax, ay, az, segtype(j))
            enddo        
            pro(i, jj) = pro(i, jj) - benergy*ntrans(i,ii)
            pro(i, jj) = dexp(pro(i, jj))

            do j = 1, l_cha
                ax = px(1, j, jj)
                ay = py(1, j, jj)
                az = pz(1, j, jj)
                
                fv = (1.0 - volprot(ax, ay, az))
                im = segtype(j)
                
                avpol_temp(ax, ay, az, im) = avpol_temp(ax, ay, az, im) + &
                     pro(i, jj) * vsol / (delta**3) / fv * ngpol(ii)
            enddo

            q_tosend = q_tosend + pro(i, jj)
            sumtrans_tosend = sumtrans_tosend + ntrans_val*pro(i, jj)

        else
           
           do j=1,longcha(ii)
            ax = px(i, j, jj) ! cada uno para su cadena...
            ay = py(i, j, jj)
            az = pz(i, j, jj)         
            pro(i, jj) = pro(i, jj) + xpot(ax, ay, az, segtype(j))
           enddo
            
           pro(i,jj) = pro(i,jj) -benergy*ntrans(i,ii) ! energy of trans bonds
           pro(i,jj) = dexp(pro(i,jj))

           do j=1,longcha(ii)
               fv = (1.0-volprot(px(i,j, jj),py(i,j, jj),pz(i,j, jj)))
               im = segtype(j)
               avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj), im)= &
               avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj), im)+pro(i, jj)*vsol/(delta**3)/fv* &
               ngpol(ii) ! ngpol(ii) has the number of chains grafted to the point ii
           enddo
            q_tosend=q_tosend+pro(i, jj)
            sumtrans_tosend = sumtrans_tosend+ntrans(i, ii)*pro(i,jj)
        endif      
        
    enddo ! Fin bucle i (configuraciones de la cadena ii)

    ! --- 4. NORMALIZACIÓN POR CADENA ---
    if (q_tosend > 0.0d0) then
        avpol_tosend = avpol_tosend + avpol_temp / q_tosend
        q(ii) = q_tosend 
        sumtrans(ii) = sumtrans_tosend / q_tosend
    endif

enddo ! Fin bucle jj (total de puntos de injerto del rank)

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

contains

integer function compare_jj_ii_i(m1, m2)
    ! NO repitas 'implicit none' aquí si ya está arriba
    type(IndexMap), intent(in) :: m1, m2
    
    if (m1%jj /= m2%jj) then
        compare_jj_ii_i = m1%jj - m2%jj
    else if (m1%ii /= m2%ii) then
        compare_jj_ii_i = m1%ii - m2%ii
    else
        compare_jj_ii_i = m1%i_conf - m2%i_conf
    end if
end function compare_jj_ii_i

recursive subroutine quicksort_idx(a, first, last)
    type(IndexMap), intent(in out) :: a(:)
    integer, intent(in) :: first, last
    
    ! Renombramos para evitar conflicto con la 'x' global de fkfun
    integer :: i_ptr, j_ptr 
    type(IndexMap) :: pivot, temp

    i_ptr = first
    j_ptr = last
    pivot = a((first + last) / 2)

    do
        do while (compare_jj_ii_i(a(i_ptr), pivot) < 0)
            i_ptr = i_ptr + 1
        end do
        do while (compare_jj_ii_i(a(j_ptr), pivot) > 0)
            j_ptr = j_ptr - 1
        end do
        
        if (i_ptr <= j_ptr) then
            temp = a(i_ptr)
            a(i_ptr) = a(j_ptr)
            a(j_ptr) = temp
            i_ptr = i_ptr + 1
            j_ptr = j_ptr - 1
        end if
        
        if (i_ptr > j_ptr) exit
    end do

    if (first < j_ptr) call quicksort_idx(a, first, j_ptr)
    if (i_ptr < last) call quicksort_idx(a, i_ptr, last)
end subroutine quicksort_idx


function find_pos_in_index(target_jj, target_ii, target_i) result(found_pos)
    integer, intent(in) :: target_jj, target_ii, target_i
    integer(kind=8) :: found_pos
    integer :: low, high, mid
    integer :: cmp

    found_pos = -1 
    low = 1
    high = total_records

    do while (low <= high)
        mid = (low + high) / 2
        
        ! Lógica de comparación jerárquica jj > ii > i
        if (sorted_idx(mid)%jj /= target_jj) then
            cmp = target_jj - sorted_idx(mid)%jj
        else if (sorted_idx(mid)%ii /= target_ii) then
            cmp = target_ii - sorted_idx(mid)%ii
        else
            cmp = target_i - sorted_idx(mid)%i_conf
        end if

        if (cmp == 0) then
            found_pos = sorted_idx(mid)%pos
            return
        else if (cmp > 0) then
            low = mid + 1
        else
            high = mid - 1
        end if
    end do
end function find_pos_in_index


end