subroutine allocatecpp

!
! Subroutine allocatecpp
!
! This routine allocates variables depedent on maxcpp 
! maxcpp is the maximum grafting points per process
!
        
        
use fields_fkfun
use conformations
use chainsdat
use solventchains
implicit none
ALLOCATE(px(cuantas, long, maxcpp))
ALLOCATE(py(cuantas, long, maxcpp))
ALLOCATE(pz(cuantas, long, maxcpp))
ALLOCATE(pro(cuantas, maxcpp))

ALLOCATE(pxsv(cuantassv, longsv))
ALLOCATE(pysv(cuantassv, longsv))
ALLOCATE(pzsv(cuantassv, longsv))


end subroutine
