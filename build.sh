#!/bin/bash

module load git cmake openmpi

LAMMPS_LIB=/lammps/build/liblammps_linux.a
LAMMPS_SRC=/lammps/src
CXX=mpicxx
FCC=mpifort

# rm -r CMakeCache.txt CMakeFiles/ Makefile bin/ cmake_install.cmake modules/ src/

cmake -DLAMMPS_LIB=$LAMMPS_LIB \
      -DLAMMPS_INCLUDE=$LAMMPS_SRC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_Fortran_COMPILER=$FCC \
      -DCMAKE_BUILD_TYPE=DEBUG \
      ..

make -j 4