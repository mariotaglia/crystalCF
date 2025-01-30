# crystalCF

MOLT-CF code for nanoparticle superlattices (originally vacuum_solv branch in crystal)



Keywords for DEFINITIONS.txt


## branched *int*

Defines if chains are linear or branched 

*int* =
0 : linear chains
1 : branched chains type 1: three ramifications following the initial segment of length long (see long below)

Expects to read a comment line and then 3 integers corresponding to the length of each of the three ramifications

2 : branched chains type 2: adds branches of one segment to part of the chain

Expects to read a comment line and then 2 integers corresponding to the first and last segment of the backbone (of length long) that have branches


## stdout *int*

Redirects output to fort.*int*
(use *int*=6 for standard output)


## sigmar *real*

Adds a random value of the surface coverage (can be used to generate noise for microphase separation)
*real* is the magnitude of the noise (in units of chains/nm^2)

### dimx *int*
### dimy *int*
### dimz *int*
#

Controls system size
*int* determines the size in units of delta


## delta *real*

Lattice discretization size
*real* is the lattice discretization size in nm

### dx *real*
### dy *real*
### dz *real*
#

Shift the lattice 

Use it to check that free energy differences are tranlationally invariant
*real* is the fraction of delta to shift the lattice


## gama *real*

Use a non-cubic cell with angle gama. 
*real* is the angle in degrees (e.g. use gama = 60.0 and cdiva = 1.633 for hexagonal)


## cdiva *real*

Non cubic cell.
Expands the c-axis of the cell by a factor *real*


## PBC *int1 int2 int3 int4 int5 int6*

Establishes the type of boundary conditions

*int1-int2* for x 
*int3-int4* for y 
*int5-int6* for z 

Possible values: 

2: Impenetrable wall
1: Periodic boundary
0: Bulk
3: Reflective boundary

## infile *int*

Controls initial guess

*int* =
0 : uses an homogeneous initial guess
1 : reads input from in.txt
-1 : reads input from in.txt, but does not solve (just does one iteration and writes to disk)
3 : same as 1 but mirrors the input in the x axis 


## nst *int*

Controls hydrophobicity

*int* is the number of hydrophobic cases to solve
expects a list with the hydrophobic strength following the "nst int" line


## benergy *real*

Energy of Gauche bonds
only works for “branched 0”

Adds an energy *real* for the gauche bonds (kBT per bond), real > 0 means that the gauche bonds are less stable than the trans ones.



## Xucutoff *real*

Cutoff of hydrophobic interactions
real is the cutoff radius of hydrophobic interactions in nm


## long *int*

Ligand chain length

*int* is the chain length pf the ligands


## lseg *real*

Segment length of the ligands
*real* is the segment length in nm
Used for chain generation


## cuantas  *int*

Number of conformations of the ligands
*int* is the number of conformations per grafting point


## vpol *real*

Segment volume
*real* is the volume of a segment in nm^3

## vtkflag *int*

Save vtk?

*int* =
1 : save vtk file
0 : do not save vtk file


### scx *int*
### scy *int*
### scz *int*
#

Supercell
Supercell for vtk file, *int* is the number of copies to expand the cell


## readchains *int*

Read conformations from file for faster initialization, only works for “branched 0”
*int* =
-1 : save chains to cadenas.dat and exit
1 : read chains from cadenas.dat
0 : generate chain conformations before running
		
## systemtype *int*

Determines the type of system

See below for the formated required for each input

*int* = 
1 : Nanoparticles with continuous brush
2 : Channel 3D with continuous brush
3 : Not in use
4 : Channel 3D with specific grafting points (includes reservoirs)
41: Same as 4, but with only one row of polymers at the middle of the channel
42: Same as 41, but with multiple rings
52: Same as 42, but with a cylinder instead of channel
6: Planar surface with polymers grafted in a rectangular array


## randominput *int*

Shifts the positions of grafting points to favor microphase separation for systemtype 
4 different shifts are implemented (see program)


## seed *int*

Random number generator for randominput only
*int* is the seed for the random number generator used for randominput = 1

## flagmu *int*

Decides a scan based on chemical potential or volume fraction of the solvent
*int* = 
0 scan volume fraction
1 scan chemical potential

see also nkp and vscan 


--------------------------------------------------------------------------------------------------------------------
vscan int	scan flag
	what to scan: 1, hamiltonian inception parameter kp, 2, attraction strength st. If 1 then the system will use the first st value.
hguess int	Hamiltonian inception for nanochannel	int is the rotational symmetry, use int = 0 for no inception
hring int	Hamiltonian inception for nanochannel	hring: int is the initial attractive zone’s distance to the wall (in xy plane)
oval real	Hamiltonian inception for nanochannel	oval: 3D shape of the growing zone 1 means spherical
nkp int	Strnght of Hamiltonian Inception loop. Followed by list.	int is the number of cases to stuy
--------------------------------------------------------------------------------------------------------------------


--------------------------------------------------------------------------------------------------------------------
# SYSTEMTYPE OPTIONS


systemtype = 1
Fixed format:
Coment line
Number of particles (N) 
Coment line
For each N, position x, position y, position z (N lines)
Coment line
For each N, radius x, radius y, radius z (N lines)
Coment line
For each N, 3x3 rotation matrix (Nx3 lines)
Coment line
For each N, surface coverage (N lines)
Coment line
For each N, charge (N lines)
Coment line
For each N, surface-polymer attraction (N lines)

systemtype 2:
Comment line
Radius of the channel (in nm)
Comment line
Surface coverage of the channel (in nm^-2)
Comment line
Surface charge of the channel inner wall (in charges/nm^-2_
Comment line
Surface-polymer attraction strength 


systemtype 41:
Comment line
Radius of the channel (in nm)
Comment line
Size of the reservoirs (in units of delta)
Comment line
Number of brushes in the tetha direction (program distributes in z direction to achieve same separation)
Comment line
Surface coverage of the channel (in nm^-2)
Comment line
Surface charge of the channel inner wall (in charges/nm^-2_
Comment line
Surface-polymer attraction strength 

systemtype 42:
Comment line
Radius of the channel (in nm)
Comment line
Size of the reservoirs (in units of delta)
Comment line
Number of brushes in the tetha direction (program distributes in z direction to achieve same separation)
Comment line
Number of rings in the z direcction
Comment line
Position of the rings in the z direction (in delta units)
Comment line
Surface charge of the channel inner wall (in charges/nm^-2_
Comment line
Surface-polymer attraction strength 



systemtype 52:
Comment line
Radius of the rod (in nm)
Comment line
Size of the reservoirs (in units of delta)
Comment line
Number of brushes in the tetha direction (program distributes in z direction to achieve same separation)
Comment line
Number of rings in the z direcction
Comment line
Position of the rings in the z direction (in delta units)
Comment line
Surface charge of the channel inner wall (in charges/nm^-2_
Comment line
Surface-polymer attraction strength 

systemtype 6:
Comment line
Number of polymers in x and y dimensions
Comment line
Surface-polymer attraction strength 



