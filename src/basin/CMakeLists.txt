# Basin module CMakeLists.txt

# Find all Fortran source files in this directory
file(GLOB BASIN_SOURCES "*.f90" "*.F90")

# Add the sources to a library that can be linked
add_library(basin_lib OBJECT ${BASIN_SOURCES})

# Apply specific compiler flags to basin sources
target_compile_options(basin_lib PRIVATE -fopenmp -std=legacy -ffree-form -ffree-line-length-none)

# Manually specify module dependencies to ensure proper build order
# This ensures modules are built in the correct order based on dependencies
set_source_files_properties(state.f90 PROPERTIES DEPENDS "kind.f90")
set_source_files_properties(event.f90 PROPERTIES DEPENDS "state.f90")
set_source_files_properties(tranevent.f90 PROPERTIES DEPENDS "event.f90")
set_source_files_properties(event_array.f90 PROPERTIES DEPENDS "event.f90 tranevent.f90")
set_source_files_properties(basin.f90 PROPERTIES DEPENDS "state.f90 event_array.f90")
set_source_files_properties(network.f90 PROPERTIES DEPENDS "basin.f90")

# Set module directory for this library's modules
set_target_properties(basin_lib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Make modules visible to other components that depend on them
target_include_directories(basin_lib PRIVATE 
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)
