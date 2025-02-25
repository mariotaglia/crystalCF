# crystalCF

----  
  
# Keywords for init_params.txt  
  
----  
# General params

## !name *str*

Sets the binary name for the minimum aL estimation and data export names.

## !particules stoichiometric ratio
n1 *int*
n2 *int*

Sets the stoichiometric ratio between particles. n1 for part1 and n2 for part2.

## !radius part1
R1 *int*

Sets the larger particle radius in nm. 

## !list gamma *list*

Sets the approximate gamma in *list* to calculate the second particle's radius.

----  
# Binary params

## !list delta bin *list*

Sets the delta list for the binary system.

## !point sum dims bin.
gamma *real* *list[int]*

Given a gamma from gamma list, replace the sum in dim calculations from aL minimum. *int* can be any positive or negative integer. 
*int* =
0 dim estimation + 0 (dim estimated from aL estimation given delta)
1 dim estimation + 1
-1 dim estimation - 1

## !num cell bin
k_bin *int*

## !aL cell bin factor *real or expression*

----
# Particle params

## !cell part *str*

Generates a cell folder for each particle with each reference.

*str* =
fcc
bcc

## !list delta part *cell* *list*

## !num cell part
k_fcc *int*
k_bcc *int*

## !aL cell part factor
fcc *real or expression*
bcc *real or expression*

----
# Output flag

## !flag generate energy vs aL curves *True or False*
