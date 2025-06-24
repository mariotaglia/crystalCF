subroutine readinput

        
! subroutine readinput
!
! Reads input from DEFINITIONS.txt
! Checks that all required inputs are entered
! In some cases, assigns defaults values if the input is absent
!

use molecules, only : benergy, vsol0
use const, only : infile, randominput, seed, seed_lig, stdout
use MPI
use ellipsoid
use chainsdat
use transform
use system
use kai
use kaist
use s2d
use channel
use superellipse
use branches
use cube
use solventchains
use mparameters_monomer
use clusters
implicit none

! Input related variables
character (len=100)  buffer,label
integer pos
integer, parameter :: fh = 15
integer ios
integer line
integer i, j
character(len=50) :: filename = 'DEFINITIONS.txt'
character basura
integer ndi
real*8 ndr
real*8 kpini, kpfin, kpstep, ikp
integer flag

stdout = 6

! not defined variables, change if any variable can take the value
ndi = -1e5
ndr = -1.0d10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check validity of input
!

seed = ndi
seed_lig = ndi
PBC = ndi
branched = ndi
sigmar = ndr
flagmu = ndi
vscan = ndi
scx = ndi
scy = ndi
scz = ndi
vtkflag = ndi
systemtype = ndi
dimx = ndi
dimy = ndi
dimz = ndi
long = ndi
longsv = ndi
cuantas = ndi
cuantassv = ndi
readchains = ndi
infile = ndi
interaction_00 = ndi
interaction_11 = ndi
randominput = ndi
cutoff = ndr
lseg = ndr
lsegkai = ndr
nst = ndi
delta = ndr
dx = ndr
dy = ndr
dz = ndr
cdiva = ndr
vsol0 = ndr
gama0 = ndr
benergy = ndr
coordinate_system = ndi
transform_type = ndi
dumpcluster = ndi
cluster_same = ndi
cutoffcluster = ndr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control file variables

line = 0
ios = 0

open(fh, file=filename)

if(rank.eq.0)write(stdout,*) 'parser:', 'Reading parameters from ', filename

! ios is negative  if an end of record condition is encountered or if
! an endfile condition was detected.  It is positive  if an error was
! detected.  ios is zero otherwise.

do while (ios == 0)

 read(fh, '(A)', iostat=ios) buffer
 if (ios == 0) then
 line = line + 1

! Find the first instance of whitespace.  Split label and data.

 pos = scan(buffer, ' ')

 label = buffer(1:pos)
 buffer = buffer(pos+1:)

 select case (label)

 case ('PBC') ! periodic boundary conditions
   read(buffer, *, iostat=ios) PBC(1),PBC(2),PBC(3),PBC(4),PBC(5),PBC(6)
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

   do j = 1,5,2
    if((PBC(j).eq.1).and.(PBC(j+1).ne.1)) then 
      write(stdout,*) 'parser:', 'Error in PBC'
      stop
    endif
    if((PBC(j+1).eq.1).and.(PBC(j).ne.1)) then
      write(stdout,*) 'parser:', 'Error in PBC'
      stop
    endif
   enddo

 case ('dumpcluster') ! dump cluster information: use 0 for no dumping, any other number N for clusters of size N
   read(buffer, *, iostat=ios) dumpcluster
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('cluster_same') ! 1: treat all particles as equivalent for clustering
   read(buffer, *, iostat=ios) cluster_same
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('cutoffcluster') ! dump cluster cutoff
   read(buffer, *, iostat=ios) cutoffcluster
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('seed') ! random seed 
   read(buffer, *, iostat=ios) seed
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('seed_lig') ! random seed_lig
   read(buffer, *, iostat=ios) seed_lig
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('stdout') ! output device, default 6 (stdout)
   read(buffer, *, iostat=ios) stdout
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('vtkflag') ! save vtk?
   read(buffer, *, iostat=ios) vtkflag
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('flagmu') ! flagmu = 0, scan solvent volume fraction; flagmu = 1, scan solvent chemical potential
   read(buffer, *, iostat=ios) flagmu
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('branched') ! used for branched chains 
   read(buffer, *, iostat=ios) branched
    
 if(branched.ne.0) then
   if(rank.eq.0)write(stdout,*) 'parser: Branched chains not longer supported because of unequal chain lengths. Stopping'
   call endall
 endif
        
  if(branched.eq.1) then
   read(fh, *) basura
   read(fh, *)longb(1), longb(2), longb(3)
 endif

 if(branched.eq.2) then
   read(fh, *) basura
   read(fh, *)longb(1), longb(2)
 endif

 case ('randominput') ! randomly shift explicit GP 
   read(buffer, *, iostat=ios) randominput
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('readchains') ! read chains from file?
   read(buffer, *, iostat=ios) readchains
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

   if(readchains.ne.0) then
     if(rank.eq.0)write(stdout,*) 'Readchains different to zero no longer supported because of unequal chain lengths. Stopping'
     call endall
   endif

 case ('dimx') ! cell dimensions
   read(buffer, *, iostat=ios) dimx
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
 case ('dimy')
   read(buffer, *, iostat=ios) dimy
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dimz')
   read(buffer, *, iostat=ios) dimz
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('scx') ! supercell for export
   read(buffer, *, iostat=ios) scx
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
 case ('scy') 
   read(buffer, *, iostat=ios) scy
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
 case ('scz')
   read(buffer, *, iostat=ios) scz
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('delta') ! discretization
   read(buffer, *, iostat=ios) delta
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dx') ! displacement of unit cell
   read(buffer, *, iostat=ios) dx
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
 case ('dy')
   read(buffer, *, iostat=ios) dy
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
 case ('dz')
   read(buffer, *, iostat=ios) dz
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('coordinate_system') ! 1: real coordinates 2: fractional coordinates in the cell 3: fractional x,y + real z
   read(buffer, *, iostat=ios) coordinate_system
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('transform_type')
   read(buffer, *, iostat=ios) transform_type
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

   select case (transform_type) ! Unit cell transformation

    case(1) ! modifies just 1 axis (c) and 1 angle (gama) 
     read(fh, *) basura
     read(fh, *) cdiva ! norm c' basis vector with respect bais vector a'
     read(fh, *) basura
     read(fh, *) gama0 ! angle between a' and b' vector
    case(2) ! Use the transformation matrix between tranformed cell to a cell of 90,90,90 angles
     read(fh, *) basura
     do j = 1,3 ! read transformation matrix
     read(fh, *) MAT(j,1), MAT(j,2), MAT(j,3)
     enddo

    endselect

case ('long') ! ligand chain length
   read(buffer, *, iostat=ios) long
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('longsv') ! solvent chain length
   read(buffer, *, iostat=ios) longsv
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('cuantas') ! ligand conformations
   read(buffer, *, iostat=ios) cuantas
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('cuantassv') ! solvent conformations
   read(buffer, *, iostat=ios) cuantassv
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('lseg') ! segment length for ligand and solvent chains
   read(buffer, *, iostat=ios) lseg
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('lsegkai') ! segment length for LJ attractive term
   read(buffer, *, iostat=ios) lsegkai
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('vsol') ! bead volume
   read(buffer, *, iostat=ios) vsol0
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('benergy') ! enegy of gauche bonds
   read(buffer, *, iostat=ios) benergy
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('vscan') ! vscan: decides what to scan 
   read(buffer, *, iostat=ios) vscan
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('sigmar') ! add random fluctuations to density
   read(buffer, *, iostat=ios) sigmar
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('infile') ! read initial guess from file
   read(buffer, *, iostat=ios) infile 
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
 
 case ('interaction_00') ! interaction between solvent beads to be use in st_matrix
   read(buffer, *, iostat=ios) interaction_00 
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('interaction_11') ! interaction between ligand beads to be used in st_matrix
   read(buffer, *, iostat=ios) interaction_11
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

case ('nkp') ! solvent volume fraction or chemical potential, depending on flagmu
  read(buffer, *, iostat=ios) kpini, kpstep, kpfin
  if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  ! Primero contar cuántos puntos hay
  if (kpini <= kpfin) then
     nkp = int((kpfin - kpini)/kpstep) + 1
  else
     nkp = int((kpini - kpfin)/kpstep) + 1
  endif

  if (allocated(kps)) deallocate(kps)
  allocate(kps(nkp))

  ikp = kpini
  if (kpini <= kpfin) then
    do i = 1, nkp
      kps(i) = ikp
      ikp = ikp + kpstep
    end do
  else
    do i = 1, nkp
      kps(i) = ikp
      ikp = ikp - kpstep
    end do
  endif

 case ('nst') ! number of st cases
   read(buffer, *, iostat=ios) nst
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   do i = 1, nst
   read(fh,*)sts(i)
   enddo 

 case ('Xucutoff') ! cut off for LJ attractions in nm
   read(buffer, *, iostat=ios) cutoff
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('systemtype')
   read(buffer, *, iostat=ios) systemtype
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

   select case (systemtype) ! TYPE OF SYSTEM
                            ! TYPE = 1 is nanoparticle crystal
                            ! TYPE = 2 is 3D channel
                            ! TYPE = 3 is a 3D channel with chains at specific conditions

    case(2, 80) ! channel, cylinder
     read(fh, *) basura
     read(fh, *) rchannel
     read(fh, *) basura
     read(fh, *) sigmac
     read(fh, *) basura
     read(fh, *) eepsc

   case(81) ! superellipse
     read(fh, *) basura
     read(fh, *) sizeX, sizeY
     read(fh, *) basura
     read(fh, *) pfactor
     read(fh, *) basura
     read(fh, *) sigmas
     read(fh, *) basura
     read(fh, *) eepss

    case(3, 4, 41)
     read(fh, *) basura
     read(fh, *) rchannel
     read(fh, *) basura
     read(fh, *) RdimZ
     read(fh, *) basura
     read(fh, *) NBRUSH ! number of brushes in the tetha direction
     read(fh, *) basura
     read(fh, *) eepsc

    case(6) ! planar surface
     read(fh, *) basura
     read(fh, *) Npolx, Npoly ! number of polymers in x and y
     read(fh, *) basura
     read(fh, *) eepsc
     NNN = 0 ! no particles

    case(7) ! cube
     read(fh, *) basura
     read(fh, *) l_cube ! lado del cubo
     read(fh, *) basura
     read(fh, *) c_cube(1), c_cube(2), c_cube(3) !posicion del centro del cubo
     read(fh, *) basura
     read(fh, *) l_pol ! numero de polimeros por lado (integer)
     read(fh, *) basura
     read(fh, *) cubeR! cubo entero = 0; tercio del cubo = 1
     read(fh, *) basura
     read(fh, *) eepsc
     NNN = 0

    case(9) !cuboctahedron
     
     read(fh, *) basura
     read(fh, *) NNN

     call allocateellCO

     read(fh, *) basura
     do j = 1, NNN
        read(fh, *) Rellf(1,j), Rellf(2,j), Rellf(3,j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'pos to',  Rellf(1,j), Rellf(2,j), Rellf(3,j)
     enddo

     read(fh, *) basura
     do j = 1, NNN
        read(fh, *) Loctall(j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'octahedron size to',  Loctall(j)
     enddo
     
     read(fh, *) basura
     do j = 1, NNN
        read(fh, *) Lcubell(j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'cube size to',  Lcubell(j)
     enddo
     
     read(fh, *) basura
     do j = 1, NNN
      read(fh, *) rotmatCO(1,1,j), rotmatCO(1,2,j), rotmatCO(1,3,j)
      read(fh, *) rotmatCO(2,1,j), rotmatCO(2,2,j), rotmatCO(2,3,j)
      read(fh, *) rotmatCO(3,1,j), rotmatCO(3,2,j), rotmatCO(3,3,j)
     enddo

     read(fh, *) basura
     do j = 1, NNN
        read(fh, *) sigma(j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'surface coverage to', sigma(j)
     enddo

     read(fh, *) basura
     do j = 1, NNN
        read(fh, *) eeps(j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'hydrophobicity to', eeps(j)
     enddo

     read(fh, '(A)') basura
     if(basura.eq."") then ! use default
       longp = long
       if(rank.eq.0)write(stdout,*) 'parser:','Using default chain lenght',long
     else ! read from list        
       do j = 1, NNN
         read(fh, *) longp(j)
       if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'chain lenght to', longp(j)
       enddo
     endif

    call COrotation

    case(42, 52) ! 42: channel, 52: rod
     read(fh, *) basura
     read(fh, *) rchannel
     read(fh, *) basura
     read(fh, *) RdimZ
     read(fh, *) basura
     read(fh, *) NBRUSH ! number of brushes in the tetha direction
     read(fh, *) basura
     read(fh, *) Nrings

     allocate (ringpos(Nrings))

     read(fh, *) basura
     do i = 1, Nrings
       read(fh, *) ringpos(i)
     enddo
      ringpos = ringpos - 0.5
    
     read(fh, *) basura
     read(fh, *) eepsc

    case(60) ! channel + particle
     read(fh, *) basura
     read(fh, *) rchannel
     read(fh, *) basura
     read(fh, *) RdimZ
     read(fh, *) basura
     read(fh, *) NBRUSH ! number of brushes in the tetha direction
     read(fh, *) basura
     read(fh, *) Nrings

     allocate (ringpos(Nrings))

     read(fh, *) basura
     do i = 1, Nrings
       read(fh, *) ringpos(i)
     enddo
      ringpos = ringpos - 0.5
    
     NNN = 1 ! only one particle

     call allocateell

     read(fh, *) basura  
     do j = 1, NNN
     read(fh, *) Rellf(1,j), Rellf(2,j), Rellf(3,j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'pos to',  Rellf(1,j), Rellf(2,j), Rellf(3,j)
     enddo

     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) Aell(1,j), Aell(2,j), Aell(3,j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'axis to',  Aell(1,j), Aell(2,j), Aell(3,j)
     enddo
     read(fh, *) basura

     do j = 1, NNN
     read(fh, *) rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
     read(fh, *) rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
     read(fh, *) rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     if(rank.eq.0) then
         write(stdout,*) 'parser:','Set particle',j,'rotation to:'
         write(stdout,*) 'parser:', rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
         write(stdout,*) 'parser:', rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
         write(stdout,*) 'parser:', rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     endif
     enddo

     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) eeps(j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'hydrophobicity to', eeps(j)
     enddo

     read(fh, '(A)') basura
     if(basura.eq."") then ! use default
       longp = long
       if(rank.eq.0)write(stdout,*) 'parser:','Using default chain lenght',long
     else ! read from list        
       do j = 1, NNN
         read(fh, *) longp(j)
       if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'chain lenght to', longp(j)
       enddo
     endif


    case(1) 
     read(fh, *) basura
     read(fh, *)NNN

     if(NNN.ne.0) then

     call allocateell
     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) Rellf(1,j), Rellf(2,j), Rellf(3,j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'pos to',  Rellf(1,j), Rellf(2,j), Rellf(3,j)
     enddo
     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) Aell(1,j), Aell(2,j), Aell(3,j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'axis to',  Aell(1,j), Aell(2,j), Aell(3,j)
     enddo
     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
     read(fh, *) rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
     read(fh, *) rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     if(rank.eq.0) then
         write(stdout,*) 'parser:','Set particle',j,'rotation to:'
         write(stdout,*) 'parser:', rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
         write(stdout,*) 'parser:', rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
         write(stdout,*) 'parser:', rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     endif
     enddo

     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) sigma(j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'surface coverage to', sigma(j)
     enddo

     read(fh, *) basura
     do j = 1, NNN
     read(fh, *) eeps(j)
     if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'hydrophobicity to', eeps(j)
     enddo

     read(fh, '(A)') basura
     if(basura.eq."") then ! use default
       longp = long
       if(rank.eq.0)write(stdout,*) 'parser:','Using default chain lenght',long
     else ! read from list        
       do j = 1, NNN
         read(fh, *) longp(j)
       if(rank.eq.0)write(stdout,*) 'parser:','Set particle',j,'chain lenght to', longp(j)
       enddo
     endif


     endif ! NNN
endselect
endselect

endif

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check validity of input
! 

if(systemtype.eq.2) then
 if((cdiva.ne.1.0).or.(gama0.ne.90.0)) then
  write(stdout,*) 'Channel works only for cdiva = 1 and gama0 = 90.0... ending'
  call MPI_FINALIZE(ierr) ! finaliza MPI
  stop
 endif
endif

if (branched.eq.1) then
 longbb = long
 long = longbb + longb(1) + longb(2) + longb(3)
endif 

if (branched.eq.2) then
 longbb = long
 long = long + longb(2) - longb(1)
endif 


if(vtkflag.eq.ndi)call stopundef('vtkflag')
if(dimx.eq.ndi)call stopundef('dimx')
if(scx.eq.ndi)call stopundef('scx')
if(scy.eq.ndi)call stopundef('scy')
if(scz.eq.ndi)call stopundef('scz')
if(dimy.eq.ndi)call stopundef('dimy')
if(dimz.eq.ndi)call stopundef('dimz')
if(ncha.eq.ndi)call stopundef('ncha')
if(long.eq.ndi)call stopundef('long')
if(longsv.eq.ndi)call stopundef('longsv')
if(cuantas.eq.ndi)call stopundef('cuantas')
if(cuantassv.eq.ndi)call stopundef('cuantassv')
if(infile.eq.ndi)call stopundef('infile')


if(interaction_00.eq.ndi) then
        interaction_00=1.0
        print*, 'interaction_00 undefined, use default value of', interaction_00
endif        
if(interaction_11.eq.ndi) then
        interaction_11=1.0
        print*, 'interaction_11 undefined, use default value of',interaction_11
endif        

if(cutoff.eq.ndr)call stopundef('Xucutoff')
if(readchains.eq.ndi)call stopundef('readchains')
if(systemtype.eq.ndi)call stopundef('systemtype')
if(nst.eq.ndi)call stopundef('nst')

if(delta.eq.ndr)call stopundef('delta')
if(dx.eq.ndr)call stopundef('dx')
if(dy.eq.ndr)call stopundef('dy')
if(dz.eq.ndr)call stopundef('dz')
if(lseg.eq.ndr)call stopundef('lseg')
if(lsegkai.eq.ndr)call stopundef('lsegkai')

if(vsol0.eq.ndr)call stopundef('vsol')
if(benergy.eq.ndr)call stopundef('benergy')
if(transform_type.eq.ndi)call stopundef('transform_type')
if(transform_type.eq.1)then
 if(gama0.eq.ndr)call stopundef('gama')
 if(cdiva.eq.ndr)call stopundef('cdiva')
endif
if(transform_type.eq.2)then
 do i=1,3
  do j=1,3
   if(MAT(i,j).eq.ndr)call stopundef('MAT')
  enddo
 enddo
endif

if(seed.eq.ndi) then
   seed = 938121
   print*, 'seed undefined, used default:', seed
endif

if(seed_lig.eq.ndi) then
   seed_lig = 14258825
   print*, 'seed_lig undefined, used default:', seed_lig
endif

if(PBC(1).eq.ndi)call stopundef('PBC')

if(flagmu.eq.ndi) then
   flagmu = 0
   print*, 'flagmu undefined, used default:', flagmu
endif

if(branched.eq.ndi) then
   branched = 0
   print*, 'branched undefined, used default:', branched
endif

if(sigmar.eq.ndr) then
   sigmar = 0.0
   print*, 'sigmar undefined, used default:', sigmar
endif

if(randominput.eq.ndi) then
   randominput = 0
   print*, 'randominput undefined, used default:', randominput
endif

if(dumpcluster.eq.ndi) then
   dumpcluster = 0
   print*, 'dumpcluster undefined, used default:', dumpcluster
endif

if(cluster_same.eq.ndi) then
   cluster_same = 0
   print*, 'cluster_same undefined, used default:', cluster_same
endif


if(cutoffcluster.eq.ndr) then
   cutoffcluster = 0.0
   print*, 'cutoffcluster undefined, used default:', cutoffcluster
endif



if(coordinate_system.eq.ndi)call stopundef('coordinate_system')


!!!!!!!! Find different lenghts for NPs !!!!!!!!!!!!!

nlongdif = 0
do j = 1, NNN
flag = 1
  do i = 1, nlongdif
     if(longp(j).eq.longdif(i)) then ! found chain lenght in previous particle 
         flag = 0 
     endif
  enddo
  if (flag.eq.1) then ! need to add a new chain lenght
     nlongdif = nlongdif + 1
     longdif(nlongdif) = longp(j)
     if(long.lt.longp(j))long=longp(j) ! long is updated to the longest chain, important for array allocation
  endif
enddo

end subroutine

subroutine stopundef(namevar)
use const
character(len=*) :: namevar
write(stdout,*) 'parser:', 'Variable ', namevar, ' is undefined '
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
end

