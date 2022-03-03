#! /bin/csh

set FC = 'gfortran'
set FC = 'ifort'
set Target = 'main.exe'

${FC} -c prm_list.f90
${FC} -c read_prm.f90
${FC} -c read_dump.f90
${FC} -c calc_poscm.f90
${FC} -c calc_inertia_tensor.f90
${FC} -c main.f90
#${FC} -c rdf_cm.f90
#${FC} -c rdf_polymer.f90
${FC} -o ${Target} *.o  -mkl
#${FC} -o ${Target} *.o  -llapack -lblas
#ifort -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback}  *.o 

rm *.o *.mod
