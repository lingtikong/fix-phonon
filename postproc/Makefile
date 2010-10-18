.SUFFIXES : .o .cpp
# compiler and flags
CC     = g++
LINK   = $(CC)
CFLAGS = -O $(DEBUG) $(UFLAG)
#
OFLAGS = -O $(DEBUG)
INC    = $(LPKINC) $(TCINC) $(SPGINC)
LIB    = $(LPKLIB) $(TCLIB) $(SPGLIB)
#
# cLapack library needed
LPKINC = -I/opt/clapack/3.2.1/include
LPKLIB = -L/opt/clapack/3.2.1/lib -lclapack -lblas -lf2c -lm
#
# Tricubic library needed
TCINC = -I/opt/tricubic/1.0/include
TCLIB = -L/opt/tricubic/1.0/lib -ltricubic
#
# spglib, used to get the irreducible q-points
UFLAG  = -DUseSPG
SPGINC = -I/opt/spglib/0.7.1/include
SPGLIB = -L/opt/spglib/0.7.1/lib -lsymspg

# Debug flags
#DEBUG = -g
#====================================================================
BASE   = phana
# executable name
EXE    = $(BASE)
#====================================================================
# source and rules
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

#====================================================================
all:  ${EXE}

${EXE}: $(OBJ)
	$(LINK) $(OFLAGS) $(OBJ) $(LIB) -o $@

clean: 
	rm -f *.o *~ *.mod ${EXE}

tar:
	rm -f ${BASE}.tar; tar -czvf ${BASE}.tar.gz *.cpp  *.h *.f *.f90 Makefile README

.f.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.f90.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<
.cpp.o:
	$(CC) $(CFLAGS) $(INC) -c $<
