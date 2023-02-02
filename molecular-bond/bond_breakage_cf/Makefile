FC = gfortran

.SUFFIXES : .f90 .F90

TARGET = a.out

SRC = prm_list.f90 \
	  read_prm.f90 \
	  xdr.F90 \
      coarse_grain.f90 \
	  read_dump.f90 \
	  calc_poscm.f90 \
	  main.f90

OBJ = prm_list.o \
	  read_prm.o \
	  xdr.o \
      coarse_grain.o \
	  read_dump.o \
	  calc_poscm.o \
	  main.o


$(TARGET): $(OBJ)
	$(FC) -fopenmp -o $@ $(OBJ) -L$(HOME)/local/xdrfile-1.1.4-gfortran -lxdrfile

.f90.o:
	$(FC) -c $<

.F90.o:
	$(FC) -c $<

all: clean $(TARGET)

clean:
	rm -f *.o *.mod

