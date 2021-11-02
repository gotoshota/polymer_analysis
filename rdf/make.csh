#! /bin/csh

set FC = 'gfortran'
set Target = 'main.exe'

${FC} -c prm_list.f90
${FC} -c read_prm.f90
${FC} -c read_dump.f90
#${FC} -c rdf.f90
${FC} -c rdf_cm.f90
#${FC} -c rdf_polymer.f90
${FC} -o ${Target} *.o 
#ifort -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback}  *.o 

rm *.o *.mod
