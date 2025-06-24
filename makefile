TARGET = crystalCF

#SRC  = modules.f90 SPmain.f90 parser.f90 init.f90 allocation.f90 allocateell.f90 3D.f90 cadenas.f90 cadenasMK.f90 fe.f90  fkfun.f90  kai.f90  kinsol.f90  pxs.f90  savetodisk.f90 rands.f90 ellipsoid.f90 dielectric.f90 monomers.definitions-onck.f90 chains.definitions.f90 sphere.f90 kapfromfile.f90

# ELECTRO
#SRC = modules.f90 SPmain.f90 channel.f90 PBC.f90 parser.f90 init.f90 allocation.f90 allocatencha.f90 allocateell.f90 allocateellCO.f90 3D.f90  allocatecpp.f90  cadenas.f90 cadenas_b.f90 cadenas_b2.f90  fe.f90  fkfun.f90  kai.f90  kinsol.f90  pxs.f90  savetodisk.f90 rands.f90 ellipsoid.f90 dielectric.f90 transform.f90 testsystem.f90 testsystemc.f90 testsystem_cube.f90 testsystemr.f90 monomers.definitions.f90 chains.definitions.f90 channel-part.f90 cube.f90 puntas.f90 cuboctahedron.f90 testsystem_cuboctahedron.f90 COrotation.f90 cylinder.f90 testsystem_cylinder.f90 superellipse.f90 testsystem_superellipse.f90

SRC = modules.f90 SPmain.f90 channel.f90 PBC.f90 parser.f90 init.f90 allocation.f90 allocatencha.f90 allocateell.f90 allocateellCO.f90 3D.f90  allocatecpp.f90  cadenas.f90 cadenas_b.f90 cadenas_b2.f90  dumpcluster.f90 fe.f90  fkfun.f90  kai.f90  kinsol.f90  pxs.f90  savetodisk.f90 rands.f90 ellipsoid.f90 transform.f90 testsystem.f90 testsystemc.f90 testsystem_cube.f90 testsystemr.f90 monomers.definitions.f90 chains.definitions.f90 channel-part.f90 cube.f90 puntas.f90 cuboctahedron.f90 testsystem_cuboctahedron.f90 COrotation.f90 cylinder.f90 testsystem_cylinder.f90 superellipse.f90 testsystem_superellipse.f90 pxssv.f90


HOST=$(shell hostname)
$(info HOST is ${HOST})

LFLAGS = -lm -L/projects/p31819/lib/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial      -Wl,-rpath,/projects/p31819/lib/lib


ifeq ($(HOST),leomisso)
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
endif

ifeq ($(HOST), santiago-desktop)
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so.1  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
endif


# some definitions
SHELL = /bin/bash
##FFLAGS= -Wunused -fbacktrace -fbounds-check # -O3
##FFLAGS= -Wunused -O3 #-fbacktrace -fbounds-check -O3
FFLAGS= -g -O3 -fallow-argument-mismatch

GIT_VERSION := $(shell git describe --abbrev=6 --dirty --always --tags)
GFLAGS=-cpp -D_VERSION=\"$(GIT_VERSION)\"

FF = mpif77 #${F90}
VER = ~/bin/crystalCF

all:	$(TARGET)

$(TARGET): $(SRC:.f90=.o)
	$(FF) -o $(TARGET) $(SRC:.f90=.o) $(LFLAGS) $(GFLAGS)
#	cp $(TARGET) $(VER)

$(SRC:.f90=.o): $(SRC)
	${FF} -c ${FFLAGS}  $(SRC) $(LFLAGS) $(GFLAGS)

install: all
	cp $(TARGET) $(VER)

clean:	
	@rm -f $(SRC:.f90=.o) $(SRC:.f90=.d) $(TARGET) *~

realclean: clean
	@rm -f .depend

depend dep:
	@$(FF)  $(CFLAGS) -MM $(SRC) > .depend 

ifeq (.depend, $(wildcard .depend))
include .depend
endif

















































