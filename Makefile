#Change LMPLIB with directory with LAMMPS static library
LMPSRC = /home/achavira/software/lammps/src/
LMPLIB = /home/achavira/software/lammps/build/liblammps_linux.a
LMPCWRAP = LAMMPS-wrapper.o
LMPFLIB = LAMMPS.o

GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)

BUILD_DIR = Obj_files
OBJS_names = LAMMPS.o LAMMPS-wrapper.o version.o bins.o define.o lmp.o queue.o clusters.o init.o funct.o dimer.o main.o akmc.o
OBJS = $(OBJS_names:%=$(BUILD_DIR)/%)
FLAGS = -O3 -g -fopenmp -std=legacy -ffree-form -I$(BUILD_DIR) #-fcheck=all
LIBS = -lgomp -lstdc++ -lm -ldl -lmpi -lmpi_cxx

TARGET = opld

# Rule to build executable
$(TARGET): $(OBJS)
	mpifort $(FLAGS) -o  $@  $^ $(LMPLIB) $(LIBS)
	mv *.mod $(BUILD_DIR)
	
# Rule to build LAMMPS interfaces object file
$(BUILD_DIR)/$(LMPCWRAP): LAMMPS-wrapper.cpp | $(BUILD_DIR)
	mpicxx -c -O3 -I $(LMPSRC) $< -o $@

$(BUILD_DIR)/$(LMPFLIB): LAMMPS.F90 | $(BUILD_DIR)
	mpifort -c -O3 $< -o $@

# Rule to build version identifier object file
$(BUILD_DIR)/version.o: version.f90
	mpifort -c $< -o $@

version.f90:
	@echo "subroutine version" > version.f90
	@echo " print*,' git version   : $(GIT_VERSION)'" >> version.f90
	@echo "end subroutine" >> version.f90

# Rule to build rest of the object files
$(BUILD_DIR)/%.o: %.f90 | $(BUILD_DIR)
	mpifort $(FLAGS) -c $< -o $@ 

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Rule to delete everything and start from scratch
clean:
	rm -rf *.o *.mod version.f90 $(BUILD_DIR) $(TARGET)

