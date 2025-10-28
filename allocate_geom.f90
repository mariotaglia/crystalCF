subroutine allocate_geom

! subroutine allocateell
!
! Allocates variables depedent on NNN for spherical particles
! NNN: Number of particles in the system
!
!
use system
use ellipsoid
use geometries
use COrotmod

!general
ALLOCATE (rotmatrix(3,3, NNN))
ALLOCATE (sigma(NNN))
ALLOCATE (eeps(NNN))
ALLOCATE (longp(NNN))
ALLOCATE (Rellf(3,NNN))
ALLOCATE (Rell(3,NNN))

! ellipsoid
ALLOCATE (Aell(3, NNNell))
ALLOCATE (AellS(3, NNNell))
ALLOCATE (AellL(3, NNNell))
ALLOCATE (AellX(3, NNNell))
ALLOCATE (AAA(3,3, NNNell))
ALLOCATE (AAAS(3,3, NNNell))
ALLOCATE (AAAL(3,3, NNNell))
ALLOCATE (AAAX(3,3, NNNell))
ALLOCATE (orient(3, NNNell))

! cuboctahedron
ALLOCATE (rotmatCO(3,3, NNNco))
ALLOCATE (Loctall(NNNco))
ALLOCATE (Lcubell(NNNco))
ALLOCATE (klocta1(NNNco))
ALLOCATE (klocta2(NNNco))
ALLOCATE (klocta3(NNNco))
ALLOCATE (klocta4(NNNco))
ALLOCATE (klocta1b(NNNco))
ALLOCATE (klocta2b(NNNco))
ALLOCATE (klocta3b(NNNco))
ALLOCATE (klocta4b(NNNco))
ALLOCATE (klcubex1(NNNco))
ALLOCATE (klcubex2(NNNco))
ALLOCATE (klcubey1(NNNco))
ALLOCATE (klcubey2(NNNco))
ALLOCATE (klcubez1(NNNco))
ALLOCATE (klcubez2(NNNco))
ALLOCATE (plane1(3,NNNco))
ALLOCATE (plane2(3,NNNco))
ALLOCATE (plane3(3,NNNco))
ALLOCATE (plane4(3,NNNco))
ALLOCATE (plane5(3,NNNco))
ALLOCATE (plane6(3,NNNco))
ALLOCATE (plane7(3,NNNco))

end subroutine
