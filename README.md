# crystalCF

MOLT-CF code for nanoparticle superlattices (originally vacuum_solv branch in crystal)



Keywords for DEFINITIONS.txt

----
# General options

## stdout *int*

Redirects output to fort.*int*
(use *int*=6 for standard output)

## seed *int*

Random number generator for randominput only
*int* is the seed for the random number generator used for randominput = 1

---
# Scan options

## infile *int*

Controls initial guess

*int* =
0 : uses an homogeneous initial guess
1 : reads input from in.txt
-1 : reads input from in.txt, but does not solve (just does one iteration and writes to disk)
3 : same as 1 but mirrors the input in the x axis 

## vscan *int*

Decides what to scan:

*int* =

1: Solvent volume fraction or chemical potential, see kp and flagmu
   It will use the first value of st
2: Attraction strength st, see nst
   It will use the first value of knp

## flagmu *int*

Decides a scan based on chemical potential or volume fraction of the solvent
*int* = 
0 scan volume fraction
1 scan chemical potential

see also nkp and vscan 

## nkp *real1* *real2* *real3*

Average solvent volume fraction in the system (if flagmu = 0) or chemical potential (if flagmu = 1)
Scans from *real1* to *real3* in steps of *real2*

## nst *int*

Controls hydrophobicity

*int* is the number of hydrophobic cases to solve
expects a list with the hydrophobic strength following the "nst int" line

---

# Output options

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

---
# Chain options

## benergy *real*

Energy of Gauche bonds
only works for “branched 0”

Adds an energy *real* for the gauche bonds (kBT per bond), real > 0 means that the gauche bonds are less stable than the trans ones.

## xucutoff *real*

Cutoff of hydrophobic interactions
real is the cutoff radius of hydrophobic interactions in nm

## long *int*

Ligand chain length

*int* is the chain length of the ligands

## cuantas  *int*

Number of conformations of the ligands
*int* is the number of conformations per grafting point

## lseg *real*

Segment length of the ligands
*real* is the segment length in nm
Used for chain generation

## lsegkai *real*

Segment length of the ligands
*real* is the segment length in nm
Used for LJ interactions

## longsv *int*

Solvent chain length

*int* is the chain length of the solvent

## cuantassv  *int*

Number of conformations of the solvent
*int* is the number of conformations per lattice site

## readchains *int*

Read conformations from file for faster initialization, only works for “branched 0”
*int* =
-1 : save chains to cadenas.dat and exit
1 : read chains from cadenas.dat
0 : generate chain conformations before running

## branched *int*

Defines if chains are linear or branched

*int* = 
0 : linear chains 
1 : branched chains type 1: three ramifications following the initial segment of length long (see long)
    Expects to read a comment line and then 3 integers corresponding to the length of each of the three ramifications
2 : branched chains type 2: adds branches of one segment to part of the chain
    Expects to read a comment line and then 2 integers corresponding to the first and last segment of the backbone (of length long) that have branches	

## randominput *int*

Shifts the positions of grafting points to favor microphase separation for systemtype 
4 different shifts are implemented (see program)

## vsol *real*

Segment volume
*real* is the volume of a segment in nm^3

## sigmar *real*

Adds a random value of the surface coverage (can be used to generate noise for microphase separation)
*real* is the magnitude of the noise (in units of chains/nm^2)

---
# Unit cell options

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

## transform_type *int*

Flag that defines the definition of the geometry of the unit cell:

*int* = 
1: Only for cubic, orthorhombic and monoclinic unit cells, expects to read:
Comment line
c/a *real*
Comment line
Angle(beta) *real*
2: For any cell, use a transformation matrix between tranformed cell to a cell of 90,90,90 angles, expects to read:
Comment line
3x3 Transformation matrix

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
	
## coordinate_system

Format of coordinates used to input NP positions:

1: real coordinates in nm 
2: fractional coordinates in the cell 
3: fractional x,y + real z

## systemtype *int*

Determines the type of system

See below for the fixed formated required for each input

*int* = 
1 : Nanoparticles with continuous brush  
2 : Channel 3D with continuous brush  
3 : Not in use  
4 : Channel 3D with specific grafting points (includes reservoirs)  
41: Same as 4, but with only one row of polymers at the middle of the channel  
42: Same as 41, but with multiple rings  
52: Same as 42, but with a cylinder instead of channel  
6: Planar surface with polymers grafted in a rectangular array  
7: Single cube  
80: Same as 2, but for a rod instead of a cylinder  
81: Superellipse (in x,y plane; infinite in z)  
9: Cubooctahedral NPs with continuous brush  

----
# Systemtype options:

systemtype 1: Nanoparticles with continuous brush
Coment line
Number of particles (N) 
Coment line
For each N, position x, position y, position z (N lines)
Coment line
For each N, radius x, radius y, radius z (N lines)
Coment line
For each N, 3x3 rotation matrix (3xN lines)
Coment line
For each N, surface coverage (N lines)
Coment line
For each N, surface-polymer attraction (N lines)

systemtype 2: Channel 3D with continuous brush
Comment line
Radius of the channel (in nm)
Comment line
Surface coverage of the channel (in nm^-2)
Comment line
Surface-polymer attraction strength 

systemtype 4: Channel 3D with specific grafting points (includes reservoirs)
Comment line
Radius of the channel (in nm)
Comment line
Size of the reservoirs (in units of delta)
Comment line
Number of brushes in the tetha direction (program distributes in z direction to achieve same separation)
Comment line
Surface coverage of the channel (in nm^-2)
Comment line
Surface-polymer attraction strength 

systemtype 41: Same as 4, but with only one row of polymers at the middle of the channel
Comment line
Radius of the channel (in nm)
Comment line
Size of the reservoirs (in units of delta)
Comment line
Number of brushes in the tetha direction (program distributes in z direction to achieve same separation)
Comment line
Surface coverage of the channel (in nm^-2)
Comment line
Surface-polymer attraction strength 

systemtype 42: Same as 41, but with multiple rings
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
Surface-polymer attraction strength 

systemtype 52: Same as 42, but with a cylinder instead of channel
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
Surface-polymer attraction strength 

systemtype 6: Planar surface with polymers grafted in a rectangular array
Comment line
Number of polymers in x and y dimensions
Comment line
Surface-polymer attraction strength 

systemtype 60: Channel (same as 42) + single NP
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
Surface-polymer attraction strength 
Coment line
Position x, position y, position z of the NP
Coment line
radius x, radius y, radius z of the NP
Coment line
3x3 rotation matrix (3 lines) of the NP
Coment line
Surface coverage of the NP
Coment line
Surface-polymer attraction of the NP

systemtype 7: Single cube
Comment line
Cube edge in nm
Comment line
x,y,z position of the center of cuve
Comment line
cubeR = use 0 for a whole cube, 1 for 1/3
Comment line
Surface-polymer attraction strength 

systemtype 80:
Comment line
Radius of the channel (in nm)
Comment line
Surface coverage of the channel (in nm^-2)
Comment line
Surface-polymer attraction strength 

systemtype 81: Superellipse
Comment line
Size in x and y directions
Comment line
pfactor (curvadure of superellipse corner)
Comment line
Surface-polymer attraction strength 
Comment line
Surface coverage (in nm^-2)
Comment line
Surface-polymer attraction strength 

systemtype 9: Cubooctahedral nanoparticles with continuous brush 
Coment line
Number of particles (N) 
Coment line
For each N, position x, position y, position z (N lines)
Coment line
For each N, octahedron size in nm (N lines)
Coment line
For each N, truncating cube size in nm (N lines)
Coment line
For each N, 3x3 rotation matrix (3xN lines)
Coment line
For each N, surface coverge (N lines)
Coment line
For each N, surface-polymer attraction (N lines)

