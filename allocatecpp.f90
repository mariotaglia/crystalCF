subroutine allocatecpp

!
! Subroutine allocatecpp
!
! This routine allocates variables depedent on maxcpp 
! maxcpp is the maximum grafting points per process
!
        
use system        
use fields_fkfun
use conformations
use chainsdat
use solventchains
implicit none

if (flag_write_pxyz.eq.1) then 
ALLOCATE(px(1, long, maxcpp))
ALLOCATE(py(1, long, maxcpp))
ALLOCATE(pz(1, long, maxcpp))

else
ALLOCATE(px(cuantas, long, maxcpp))
ALLOCATE(py(cuantas, long, maxcpp))
ALLOCATE(pz(cuantas, long, maxcpp))
endif

ALLOCATE(pro(cuantas, maxcpp))

ALLOCATE(pxsv(cuantassv, longsv))
ALLOCATE(pysv(cuantassv, longsv))
ALLOCATE(pzsv(cuantassv, longsv))


end subroutine
