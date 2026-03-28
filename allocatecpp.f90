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
ALLOCATE(px(1, long, maxcpp))
ALLOCATE(py(1, long, maxcpp))
ALLOCATE(pz(1, long, maxcpp))
ALLOCATE(pro(cuantas, maxcpp))

ALLOCATE(pxsv(cuantassv, longsv))
ALLOCATE(pysv(cuantassv, longsv))
ALLOCATE(pzsv(cuantassv, longsv))


end subroutine
