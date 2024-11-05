#Change LMPLIB with directory with LAMMPS static library
LMPSRC = /home/achavira/software/lammps/src/
LMPLIB = /home/achavira/software/lammps/build/liblammps_linux.a
LMPCWRAP = ./LAMMPS-wrapper.o
LMPFLIB = ./LAMMPS.o

OBJS = LAMMPS.o LAMMPS-wrapper.o bins.o define.o lmp.o queue.o helper.o ptm.o grain.o reflat.o feat.o init.o funct.o dimer.o main.o akmc.o
FLAGS = -O3 -g -fopenmp -std=legacy -ffree-form #-fcheck=all
LIBS = -lgomp -lstdc++ -lm -ldl -lmpi_cxx

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

