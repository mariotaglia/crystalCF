subroutine allocation

use system
use fields_fkfun
use conformations
use chainsdat
use kinsol
use results
use ematrix
use mkinsol
use ellipsoid
use MPI
use kai
use mparameters_monomer
use mmask
use solventchains
implicit none


allocate(mask(dimx,dimy,dimz))

! fields_fkfun
ALLOCATE(xtotal(dimx, dimy, dimz, 0:N_poorsol))
ALLOCATE(xh(dimx, dimy, dimz))
ALLOCATE(rhosv(dimx, dimy, dimz))
ALLOCATE(qsv(dimx, dimy, dimz))
ALLOCATE(sumprolnpro(dimx, dimy, dimz))
ALLOCATE(sumprotrans(dimx, dimy, dimz))

! kinsol
ALLOCATE (xflag(eqs*dimx*dimy*dimz))
ALLOCATE (xpar(dimx*dimy*dimz))

! results
ALLOCATE (avpol(dimx, dimy, dimz, N_monomer))
ALLOCATE (xtotalsum(dimx, dimy, dimz))

! ematrix
ALLOCATE (volprot(dimx,dimy,dimz))
ALLOCATE (volprot1(dimx,dimy,dimz))
ALLOCATE (voleps(dimx,dimy,dimz))
ALLOCATE (voleps1(dimx,dimy,dimz))

! mkinsol
ALLOCATE (pp(eqs*dimx*dimy*dimz))

! chainsdat
allocate(in1(long,3))
allocate(in1sv(longsv,3))
allocate(cpp(size))
allocate(cppini(size))
end subroutine
