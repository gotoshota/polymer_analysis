#! /bin/csh

set FC = 'ifort'
set FC = 'gfortran'
set Target = 'main.exe'

${FC} -c prm_list.f90
${FC} -c read_prm.f90
${FC} -c xdr.F90
${FC} -c read_dump.f90
${FC} -c calc_poscm.f90
${FC} -c main.f90
${FC} -o ${Target} *.o -lxdrfile
#ifort -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback}  *.o 

rm *.o *.mod
