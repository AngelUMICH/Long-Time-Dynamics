# Defect Rate-Based Long-Time Dynamics (DRLD)

## Building Code

### Requirements

The LAMMPS **static** library must be built with the support of OpenMPI and must include any necessary packages needed to run MD such as the MANYBODY pakage (enables support for eam force styles). Note the paths to the library and source.

### Next Steps

After building the LAMMPS dependancy CMake will be used to compile the DRLD project. The project is only supported operate in parallel with the OpenMPI paradigm, so be sure to provide valid mpicxx and mpifort compilers. Additionally, the following CMake flags are required to build:

- LAMMPS_LIB <path to static library>
- LAMMPS_INCLUDE <path to LAMMPS src>

### How to Use OPLD

Essentially, there are at least 4 files required for the program to run correctly. These files describe parameters, approximate defect migration energies, a LAMMPS script which provides the atomic system, and a reference lattice system in xyz format. User might need to suppply additional files to describe atom potential or read atomic system from data file for example.

Use the MPI executable launcher with the name of the DRLD executable as the only argument.

### Examples

Explore the examples folder to get a sense of how to use the DRLD code.