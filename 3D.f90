subroutine solve(flagcrash)

! Subroutine solve
!  
! This routine sets the initial guess and calls to kinsol solver
! 
! MPI implementation 
! The master process calls kinsol, which calls fkfun
! Slaves process call fkfun directly
! If the solver does not converge in that iteration, repeat
! If a solution is found, stop and store solution
!


use system, only : dimx, dimy, dimz, eqs
use const, only : error, stdout, infile
use kinsol, only : ier, xflag, norma, iter
use MPI
use mparameters_monomer, only : N_poorsol

implicit none
external fcn
integer i, ix, iy, iz, ip
integer flagcrash

!-----  kinsol solution variables -----------

real*8 x1(eqs*dimx*dimy*dimz),xg1(eqs*dimx*dimy*dimz)
real*8 f(eqs*dimx*dimy*dimz)
       
integer ncells

! Volumen fraction
real*8 xh(dimx, dimy, dimz), xtotal(dimx,dimy,dimz,N_poorsol)

real*8 temp

! MPI
integer tag, source
parameter(tag = 0)
integer err
integer ier_tosend
double  precision norma_tosend

! number of equations
ncells = dimx*dimy*dimz

! Initial guess
if((infile.eq.2).or.(infile.eq.-1).or.(infile.eq.3)) then
  do i = 1, eqs*ncells  
      xg1(i) = xflag(i)     
      x1(i) = xflag(i)
  enddo
endif

if(infile.eq.0) then
  do i=1,ncells
    xg1(i)=0.1
    x1(i)=0.1
  enddo

  do i = ncells+1,(N_poorsol+1)*ncells
    xg1(i)=0.1
    x1(i)=0.1
  enddo
endif ! infile

!--------------------------------------------------------------
! Solve               
!--------------------------------------------------------------

! master process 
if(rank.eq.0) then ! only master process calls kinsol solver
   iter = 0
   write(stdout,*) 'solve: Enter solver ', eqs*ncells, ' eqs'

   if(infile.ge.0) then
    call call_kinsol(x1, xg1, ier)
   endif
   if(infile.eq.-1) then
    call fkfun(x1, f, ier)
   endif
   flagsolver = 0
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
endif
  
! Slave processes
if(rank.ne.0) then
  do
     flagsolver = 0
     source = 0
     CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
     if(flagsolver.eq.1) then
        call call_fkfun(x1) ! the solver hasn't converged yet => slaves call fkfun 
     endif ! flagsolver
     if(flagsolver.eq.0) exit ! the solver has converged => slaves exit the loop
   enddo
endif

! Recover ier and norm, so the slaves know if the solver converged

! Master
if (rank.eq.0) then
   norma_tosend = norma
   ier_tosend = ier
   CALL MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
   CALL MPI_BCAST(ier_tosend,1, MPI_INTEGER,0,MPI_COMM_WORLD,err)
endif

! Slaves
if (rank.ne.0) then
   CALL MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
   CALL MPI_BCAST(ier_tosend, 1, MPI_INTEGER,0,MPI_COMM_WORLD,err)
   norma = norma_tosend
   ier = ier_tosend
endif

! Recover xh and xtotal from x1 vector in last iteration of the solver (these variables are not global)
do ix=1,dimx
   do iy=1,dimy
      do iz=1,dimz
       xh(ix,iy,iz)=x1(ix+dimx*(iy-1)+dimx*dimy*(iz-1))

       do ip=1, N_poorsol
          xtotal(ix,iy,iz,ip)=x1(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells)
       enddo
      enddo
   enddo  
enddo

! Check if the solution converged, if not, try to recover
if(infile.ne.-1) then
  if((ier.lt.0).or.(.not.((norma.gt.0).or.(norma.lt.0))).or.(norma.gt.error)) then ! not converged
    if(rank.eq.0)write(stdout,*) 'solve: Error in solver: ', ier
    if(rank.eq.0)write(stdout,*) 'solve: norm ', norma
    flagcrash = 1
    return
  endif
endif    

! Converged, store xflag to be used in next call to the solver
do i = 1, eqs*ncells
  xflag(i) = x1(i) 
enddo
infile = 2 ! use xflag in next call to the solver
flagcrash = 0
return

end subroutine


