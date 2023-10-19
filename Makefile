#Change LMPLIB with directory with LAMMPS static library
LMPSRC = /home/achavira/software/lammps/src/
LMPLIB = /home/achavira/software/lammps/build/liblammps_linux.a
LMPCWRAP = ./LAMMPS-wrapper.o
LMPFLIB = ./LAMMPS.o

OBJS = LAMMPS.o LAMMPS-wrapper.o BINS.o DEFINE.o LMP.o INIT.o FUNCT.o DIMER.o MAIN.o AKMC.o
FLAGS = -O3 -g -fopenmp -std=legacy -ffree-form #-Wall -fcheck=all
LIBS = -lgomp -lstdc++ -lm -ldl -ljpeg -lpng

opld: $(OBJS)
	mpifort $(FLAGS) -o  $@  $^ $(LMPLIB) $(LIBS)
	@if [ ! -d Obj_files ]; then mkdir Obj_files; fi
	@mv *.o *.mod Obj_files
	
$(LMPCWRAP): LAMMPS-wrapper.cpp
	mpicxx -c -O3 -I $(LMPSRC) $< -o $@

$(LMPFLIB): LAMMPS.F90
	mpifort -c -O3 $< -o $@

%.o: %.f*
	mpifort $(FLAGS) -c $< -o $@ 

clean:
	rm -rf *.o *.mod Obj_files

format:
	fprettify --indent 2 --enable-decl --strict-indent --c-relations --whitespace-decl 1 *.f LMP.f90

