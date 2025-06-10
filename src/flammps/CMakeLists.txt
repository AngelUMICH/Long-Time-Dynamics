# LAMMPS interface module CMakeLists.txt

# Create LAMMPS Fortran module library
add_library(lammps_mod OBJECT LAMMPS.F90)
set_source_files_properties(LAMMPS.F90 PROPERTIES COMPILE_FLAGS "-ffree-line-length-none")

# Make sure module files are generated in the central module directory
set_target_properties(lammps_mod PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Export the module directory to parent scope so other targets can find it
set(LAMMPS_MODULE_DIR ${CMAKE_Fortran_MODULE_DIRECTORY} PARENT_SCOPE)

# Create LAMMPS wrapper object library
add_library(lammpswrap_obj OBJECT LAMMPS-wrapper.cpp)
target_include_directories(lammpswrap_obj PRIVATE ${LAMMPS_INCLUDE})
