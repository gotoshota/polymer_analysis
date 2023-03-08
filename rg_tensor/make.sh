#! /bin/bash

FC='gfortran'
#FC='ifort'
Target='main.exe'

${FC} -c prm_list.f90
${FC} -c read_prm.f90
${FC} -c xdr.F90
${FC} -c read_dump.f90
${FC} -c calc_poscm.f90
${FC} -c calc_inertia_tensor.f90
${FC} -c main.f90
${FC} -o ${Target} *.o -llapack -lblas -L${HOME}/local/xdrfile-1.1.4-gfortran -lxdrfile
#${FC} -o ${Target} *.o -mkl
#ifort -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback}  *.o 
#ifort -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback *.o
rm *.o *.mod
