subroutine COrotation

use ellipsoid
use COrotMod

implicit none

integer i
real*8 nplane1(3), nplane2(3), nplane3(3), nplane4(3)
real*8 nplane5(3), nplane6(3), nplane7(3)
real*8 klocta1v(3), klocta2v(3)
real*8 klcubex1v(3), klcubex2v(3)
real*8 klcubey1v(3), klcubey2v(3)
real*8 klcubez1v(3), klcubez2v(3)

nplane1 = (/1.0,1.0,1.0/)
nplane2 = (/-1.0,1.0,1.0/)
nplane3 = (/1.0,-1.0,1.0/)
nplane4 = (/-1.0,-1.0,1.0/)
nplane5 = (/1.0,0.0,0.0/)
nplane6 = (/0.0,1.0,0.0/)
nplane7 = (/0.0,0.0,1.0/)

do i=1, NNN

! Plane x + y + z = k
plane1(:,i) = MATMUL(rotmatCO(:,:,i),nplane1)
! Plane - x + y + z = k
plane2(:,i) = MATMUL(rotmatCO(:,:,i),nplane2)
! Plane x - y + z = k
plane3(:,i) = MATMUL(rotmatCO(:,:,i),nplane3)
! Plane - x - y + z = k
plane4(:,i) = MATMUL(rotmatCO(:,:,i),nplane4)
! Plane x = k
plane5(:,i) = MATMUL(rotmatCO(:,:,i),nplane5)
! Plane y = k
plane6(:,i) = MATMUL(rotmatCO(:,:,i),nplane6)
! Plane z = k
plane7(:,i) = MATMUL(rotmatCO(:,:,i),nplane7)

! Plane constants

klocta1v(1) = 0.0
klocta1v(2) = 0.0
klocta1v(3) = Loctall(i)/2.0

klocta2v(1) = 0.0
klocta2v(2) = 0.0
klocta2v(3) = -Loctall(i)/2.0

klocta1(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klocta1v),plane1(:,i)) !VERIFICARRRRR
klocta1b(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klocta2v),plane1(:,i))
klocta2(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klocta1v),plane2(:,i))
klocta2b(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klocta2v),plane2(:,i))
klocta3(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klocta1v),plane3(:,i))
klocta3b(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klocta2v),plane3(:,i))
klocta4(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klocta1v),plane4(:,i))
klocta4b(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klocta2v),plane4(:,i))

klcubex1v(1) = Lcubell(i)/2.0
klcubex1v(2) = 0.0
klcubex1v(3) = 0.0
klcubex1(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klcubex1v),plane5(:,i))

klcubex2v(1) = -Lcubell(i)/2.0
klcubex2v(2) = 0.0
klcubex2v(3) = 0.0
klcubex2(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klcubex2v),plane5(:,i))

klcubey1v(1) = 0.0
klcubey1v(2) = Lcubell(i)/2.0
klcubey1v(3) = 0.0
klcubey1(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klcubey1v),plane6(:,i))

klcubey2v(1) = 0.0
klcubey2v(2) = -Lcubell(i)/2.0
klcubey2v(3) = 0.0
klcubey2(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klcubey2v),plane6(:,i))

klcubez1v(1) = 0.0
klcubez1v(2) = 0.0
klcubez1v(3) = Lcubell(i)/2.0
klcubez1(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klcubez1v),plane7(:,i))

klcubez2v(1) = 0.0
klcubez2v(2) = 0.0
klcubez2v(3) = -Lcubell(i)/2.0
klcubez2(i) = DOT_PRODUCT(MATMUL(rotmatCO(:,:,i),klcubez2v),plane7(:,i))

enddo

end subroutine COrotation
