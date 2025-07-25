# Core module CMakeLists.txt

# Find all Fortran source files in this directory
file(GLOB CORE_SOURCES "*.f90" "*.F90")

# Add the sources to a library that can be linked
add_library(core_lib OBJECT ${CORE_SOURCES})

# Apply specific compiler flags to core sources
target_compile_options(core_lib PRIVATE -fopenmp -std=legacy -ffree-form -ffree-line-length-none)

# Apply special flags for specific files if needed
set_source_files_properties(version.F90 PROPERTIES COMPILE_FLAGS "-cpp")

# Handle file-specific dependencies
set_source_files_properties(lmp.f90 PROPERTIES OBJECT_DEPENDS 
    "${CMAKE_Fortran_MODULE_DIRECTORY}/lammps.mod"
)

# Add explicit dependencies
add_dependencies(core_lib lammps_mod basin_lib)

# Make sure core_lib can find all required module files
target_include_directories(core_lib PRIVATE 
    ${CMAKE_Fortran_MODULE_DIRECTORY}  # Use the central module directory
    ${CMAKE_CURRENT_SOURCE_DIR}/../basin  # Path to basin source directory
)

# Set module directory for this library's modules
set_target_properties(core_lib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
