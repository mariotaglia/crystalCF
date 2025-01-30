module solventchains
integer longsv
integer cuantassv
integer*2, allocatable :: pxsv(:,:)
integer*2, allocatable :: pysv(:,:)
integer*2, allocatable :: pzsv(:,:)
integer*2, allocatable :: ntranssv(:)
integer ingsv ! number of transs in current chain
real*8, ALLOCATABLE :: in1sv(:,:)  ! segment positions 
end module


module mparameters_monomer
integer N_poorsol ! number of different kais
integer N_monomer ! number of different monomer types
real*8, allocatable :: st_matrix(:,:) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
integer, allocatable :: hydroph(:) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent
real*8 interaction_00, interaction_11
endmodule mparameters_monomer

module branches
integer longb(3), longbb
integer branched
integer indexncha
endmodule

module system 
integer systemtype
integer transform_type
integer coordinate_system
integer vscan
real*8 delta
real*8 dx,dy,dz
real*8 cdiva
integer  dimx 
integer  dimy 
integer  dimz 
integer PBC(6)
integer vtkflag
integer flagmu
integer eqs ! number of set of equations 
endmodule

module ematrix
use system
real*8, allocatable :: volprot(:,:,:)
real*8, allocatable :: volprot1(:,:,:)
real*8, allocatable :: voleps(:,:,:)
real*8, allocatable :: voleps1(:,:,:)
integer, parameter :: maxvolx = 50000
real*8 volx(maxvolx)
real*8 com(maxvolx,3)
integer p0(maxvolx,3)
end module

module rotchain
use ematrix, only : maxvolx
real*8 rotangle(maxvolx)
endmodule

module channel
real*8 rchannel
real*8 originc(2)
real*8 sigmac, eepsc, sigmar
integer NBRUSH
integer RdimZ ! size of reservoirs in delta units
integer Nrings ! number of rings for systemtype = 42
real*8, allocatable :: ringpos(:) ! position along the pore
integer Npolx, Npoly
endmodule

module superellipse
real*8 sizeX, sizeY, pfactor
real*8 originc(2)
real*8 voriginc(3)
real*8 eepss, sigmas, sigmars
endmodule

module s2d
integer scx,scy,scz
endmodule

module mkinsol
double precision, allocatable :: pp(:)
endmodule

module montecarlo
real*8 free_energy
endmodule

module chainsdat
integer cuantas 
integer, allocatable :: newcuantas(:)
integer long
integer, allocatable :: segtype(:) ! sequence of the chain 
integer ncha 
real*8, ALLOCATABLE :: in1(:,:)  ! segment positions 
real*8, ALLOCATABLE :: posicion(:,:) ! posicion graft de la cadena ncha
real*8, ALLOCATABLE :: ngpol(:) ! posicion graft de la cadena ncha
integer, ALLOCATABLE :: cpp(:)
integer, ALLOCATABLE :: cppini(:)
integer maxcpp
real*8 lseg
real*8 lsegkai
integer readchains
integer ing ! number of transs in current chain
endmodule

module molecules
use system
real*8 vsol
real*8 vsol0
real*8 vsalt
real*8 benergy
endmodule

module cube
real*8 l_cube
real*8 c_cube(3)
integer l_pol
integer cubeR
endmodule

module mmask
real*8, allocatable :: mask(:,:,:)
endmodule

module kaist
integer hguess
real*8 hring
real*8 oval
integer nkp
real*8 kp
real*8 kps(1000)

integer nst
real*8 st
real*8 sts(100)
endmodule

module fields_fkfun
use system
use chainsdat
real*8, allocatable :: xtotal(:, :, :, :) ! xtotal para poor solvent
real*8 musolv ! solvent chem pot
real*8, allocatable :: q(:)
real*8, allocatable :: qsv(:,:,:)
real*8, allocatable :: sumprolnpro(:,:,:) ! for fe calc
real*8, allocatable :: sumprotrans(:,:,:) ! for fe calc
real*8, allocatable :: rhosv(:,:,:) 
real*8, allocatable :: sumtrans(:)
real*8 sumtranssv
real*8, allocatable :: pro(:,:)
real*8 prosv 
real*8, allocatable :: xh(:, :, :)
real*8 phisolv ! solvent density
real*8 shift
endmodule

module conformations
integer*2, allocatable :: px(:,:,:)
integer*2, allocatable :: py(:,:,:)
integer*2, allocatable :: pz(:,:,:)
integer*2, allocatable :: ntrans(:,:)
endmodule

module MPI
include 'mpif.h' ! librerias MPI
integer rank, size, ierr
integer flagsolver
endmodule

module kinsol
use system
integer iter
integer, parameter :: maxiters = 5000
integer *4 ier ! Kinsol error flag
integer *8 neq ! Kinsol number of equations
real*8 norma
real*8, ALLOCATABLE :: xflag(:) 
real*8, ALLOCATABLE :: xpar(:) 
endmodule

module const
real*8 pi 
real*8, parameter :: Na = 6.02d23 
integer seed
real*8 error  ! para comparar con la norma...
real*8 errel
integer itmax
integer infile
integer randominput
integer stdout
endmodule

module kai
integer Xulimit
real*8 cutoff
real*8, allocatable :: Xu(:,:,:)
real*8 sumXu
endmodule

module results
use system
real*8, allocatable :: avpol(:,:,:,:) ! avpol ix iy iz im
real*8, allocatable :: xtotalsum(:,:,:)
endmodule


module ellipsoid
integer NNN
real*8, allocatable :: rotmatrix(:,:,:)
real*8, allocatable :: Aell(:,:)
real*8, allocatable :: AellS(:,:)
real*8, allocatable :: AellL(:,:)
real*8, allocatable :: AellX(:,:)
real*8, allocatable :: AAA(:,:,:)
real*8, allocatable :: AAAS(:,:,:)
real*8, allocatable :: AAAL(:,:,:)
real*8, allocatable :: AAAX(:,:,:)
real*8, allocatable :: Rell(:,:)
real*8, allocatable :: Rellf(:,:)
real*8, allocatable :: orient(:,:)
real*8, allocatable :: sigma(:)
real*8, allocatable :: eeps(:)

! cubooctahedron only
real*8, allocatable :: Loctall(:)
real*8, allocatable :: Lcubell(:)
real*8, allocatable :: rotmatCO(:,:,:)

end module

module COrotMod
real*8, allocatable :: klocta1(:)
real*8, allocatable :: klocta2(:)
real*8, allocatable :: klocta3(:)
real*8, allocatable :: klocta4(:)
real*8, allocatable :: klocta1b(:)
real*8, allocatable :: klocta2b(:)
real*8, allocatable :: klocta3b(:)
real*8, allocatable :: klocta4b(:)
real*8, allocatable :: klcubex1(:)
real*8, allocatable :: klcubex2(:)
real*8, allocatable :: klcubey1(:)
real*8, allocatable :: klcubey2(:)
real*8, allocatable :: klcubez1(:)
real*8, allocatable :: klcubez2(:)
real*8, allocatable :: plane1(:,:)
real*8, allocatable :: plane2(:,:)
real*8, allocatable :: plane3(:,:)
real*8, allocatable :: plane4(:,:)
real*8, allocatable :: plane5(:,:)
real*8, allocatable :: plane6(:,:)
real*8, allocatable :: plane7(:,:)

end module

module transform
real*8 gama0
real*8 MAT(3,3)
real*8 TMAT(3,3)
real*8 IMAT(3,3)
endmodule

