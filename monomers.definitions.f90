subroutine monomer_definitions

use MPI
use mparameters_monomer

implicit none

N_poorsol = 1 ! number of different kais
N_monomer = 1

ALLOCATE (st_matrix(0:N_poorsol, 0:N_poorsol)) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
ALLOCATE (hydroph(0:N_monomer)) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent

st_matrix(0,0) = interaction_00
st_matrix(0,1) = sqrt(interaction_00*interaction_11)
st_matrix(1,0) = sqrt(interaction_00*interaction_11)
st_matrix(1,1) = interaction_11

hydroph(0) = 0
hydroph(1) = 1

end

