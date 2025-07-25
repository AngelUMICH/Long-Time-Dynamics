cmake_minimum_required(VERSION 3.12)
project(LongTimeDynamics LANGUAGES Fortran CXX)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Debug, Release)" FORCE)
endif()

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Make sure the module directory exists
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Find MPI package
find_package(MPI REQUIRED Fortran)

# Require user to specify LAMMPS static library path
set(LAMMPS_LIB "" CACHE FILEPATH "Path to the LAMMPS static library (e.g., liblammps_linux.a)")
if(NOT LAMMPS_LIB)
  message(FATAL_ERROR "You must set LAMMPS_LIB to the path of your LAMMPS static library (e.g., -DLAMMPS_LIB=/path/to/liblammps_linux.a)")
endif()

# Require user to specify LAMMPS include directory
set(LAMMPS_INCLUDE "" CACHE PATH "Path to LAMMPS source include directory (e.g., /path/to/lammps/src)")
if(NOT LAMMPS_INCLUDE)
  message(FATAL_ERROR "You must set LAMMPS_INCLUDE to the path of your LAMMPS source directory (e.g., -DLAMMPS_INCLUDE=/path/to/lammps/src)")
endif()

# Define the Git version hash for compiled code
execute_process(
  COMMAND git describe --abbrev=5 --dirty --always --tags
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_VERSION
  RESULT_VARIABLE GIT_RESULT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)

if(NOT GIT_RESULT EQUAL 0)
  set(GIT_VERSION "unknown")
endif()

# Add subdirectories in order of dependencies
# Basin modules contain basic data structures
add_subdirectory(src/basin)
# LAMMPS module directory must be built before core
add_subdirectory(src/flammps)
# Core modules depend on both basin and flammps
add_subdirectory(src/core)

# Create the main executable
add_executable(opld 
  $<TARGET_OBJECTS:basin_lib>
  $<TARGET_OBJECTS:lammps_mod>
  $<TARGET_OBJECTS:core_lib> 
  $<TARGET_OBJECTS:lammpswrap_obj>
)

# Include directories for the main target
target_include_directories(opld PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core
    ${CMAKE_CURRENT_SOURCE_DIR}/src/basin
    ${CMAKE_CURRENT_SOURCE_DIR}/src/flammps
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    ${LAMMPS_INCLUDE}
)

# Make sure the executable depends on all component libraries
add_dependencies(opld basin_lib lammps_mod core_lib lammpswrap_obj)

# Set link libraries for the main executable
target_link_libraries(opld PRIVATE 
    MPI::MPI_Fortran 
    ${LAMMPS_LIB} 
    gomp 
    stdc++ 
    m 
    dl
)

# Compiler options for main executable
target_compile_options(opld PRIVATE 
    -fopenmp 
    -std=legacy 
    -ffree-form
    -ffree-line-length-none
)

target_compile_definitions(opld PRIVATE GIT_VERSION_F=\"${GIT_VERSION}\")

# Add configuration-specific compiler flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(opld PRIVATE -g -O0 -Wall -fbounds-check)
    target_compile_options(basin_lib PRIVATE -g -O0 -Wall -fbounds-check)
    target_compile_options(core_lib PRIVATE -g -O0 -Wall -fbounds-check)
    target_compile_options(lammps_mod PRIVATE -g -O0 -Wall -fbounds-check)
else()
    target_compile_options(opld PRIVATE -O3)
    target_compile_options(core_lib PRIVATE -O3)
    target_compile_options(lammps_mod PRIVATE -O3)
endif()