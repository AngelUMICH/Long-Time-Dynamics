

OBJS = DEFINE.o INIT.o FUNCT.o FORCE.o DIMER.o MAIN.o OPT.o AKMC.o FINAL.o

FLAGS = -O3 -fopenmp -std=legacy

%.o: %.f
	mpifort $(FLAGS) -o $@ -c $<

seakmc.exe: $(OBJS)
	mpifort $(FLAGS) -Wall -fcheck=all  -o  $@  $^
	@if [ ! -d Obj_files ]; then mkdir Obj_files; fi
	@mv *.o *.mod Obj_files

clean:
	rm -f *.o *.mod 

