#! /bin/sh

FC='ifort'
#FC='gfortran'
Target='main.exe'

${FC} -c prm_list.f90
${FC} -c read_prm.f90
${FC} -c read_dump.f90
${FC} -c main.f90
${FC} -o ${Target} *.o 
#ifort -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback}  *.o 

rm *.o *.mod
