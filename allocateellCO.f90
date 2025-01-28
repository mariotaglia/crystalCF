subroutine allocateellCO

use system 
use ellipsoid
use COrotmod

! cuboctahedron
ALLOCATE (Rellf(3,NNN))
ALLOCATE (Rell(3,NNN))
ALLOCATE (Loctall(NNN))
ALLOCATE (Lcubell(NNN))
ALLOCATE (rotmatCO(3,3,NNN))
ALLOCATE (echarge(NNN))
ALLOCATE (sigma(NNN))
ALLOCATE (eeps(NNN))
ALLOCATE (klocta1(NNN))
ALLOCATE (klocta2(NNN))
ALLOCATE (klocta3(NNN))
ALLOCATE (klocta4(NNN))
ALLOCATE (klocta1b(NNN))
ALLOCATE (klocta2b(NNN))
ALLOCATE (klocta3b(NNN))
ALLOCATE (klocta4b(NNN))
ALLOCATE (klcubex1(NNN))
ALLOCATE (klcubex2(NNN))
ALLOCATE (klcubey1(NNN))
ALLOCATE (klcubey2(NNN))
ALLOCATE (klcubez1(NNN))
ALLOCATE (klcubez2(NNN))
ALLOCATE (plane1(3,NNN))
ALLOCATE (plane2(3,NNN))
ALLOCATE (plane3(3,NNN))
ALLOCATE (plane4(3,NNN))
ALLOCATE (plane5(3,NNN))
ALLOCATE (plane6(3,NNN))
ALLOCATE (plane7(3,NNN))

end subroutine
