cmake_minimum_required(VERSION 3.10)
project(Sparse_Interp)

################################################################################

### Configuration
set(SPARSE_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${THIRD_PARTY_DIR}/Catch2/contrib)
set(SPARSE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/sparse_interp)

include(Warnings)
include(UseColors)
include(${PROJECT_NAME}Dependencies)
include(PrependCurrentPath)
include(Sparse_InterpUtils)

igl_include(core)
igl_include(predicates)
igl_include(restricted triangle)

add_subdirectory(${SPARSE_EXTERNAL}/openmesh/OpenMesh-11.0.0)
add_subdirectory(${SPARSE_EXTERNAL}/TinyAD)

############################
# B spline library.
# igl -> splines
include(B_spline/SOURCE.cmake)
prepend_current_path(SPLINE_SOURCES)
SI_copy_headers(${SPLINE_SOURCES})
SI_set_source_group(${SPLINE_SOURCES})  
add_library(splines ${SPLINE_SOURCES})
target_include_directories(splines PUBLIC ${PROJECT_BINARY_DIR}/include)
target_link_libraries(splines PUBLIC
    igl::core
    igl_restricted::triangle 
    Eigen3::Eigen 
    OpenMeshCore 
    OpenMeshTools
)
target_compile_definitions(splines PUBLIC
    SI_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/app/meshes_new/")
###############################

# the comparison methods
include(comparisons/compare.cmake)
add_library(compare  ${COMPARE_SOURCES})
SI_copy_headers(${COMPARE_SOURCES})
target_include_directories(compare PUBLIC ${PROJECT_BINARY_DIR}/include)
target_link_libraries(compare PUBLIC splines)

################################

set(FUNCTIONS_SOURCES 
    app/functions.hpp)
add_library(functions ${FUNCTIONS_SOURCES})
SI_copy_headers(${FUNCTIONS_SOURCES})
target_include_directories(functions PUBLIC ${PROJECT_BINARY_DIR}/include)
target_link_libraries(functions PUBLIC splines)

################################


add_executable(surfGui_bin
    src/main.cpp
    src/methods.cpp
)

target_include_directories(surfGui_bin PUBLIC 
    ${PROJECT_BINARY_DIR}/include
    ${SPARSE_EXTERNAL}/TinyAD/include
    )
target_include_directories(splines PRIVATE
    ${CMAKE_BINARY_DIR}/_deps/predicates-src
)

add_library(polyscope::polyscope ALIAS polyscope)
add_library(TinyAD::TinyAD ALIAS TinyAD)

target_link_libraries(surfGui_bin 
    PUBLIC 
        splines
        compare
        igl::core
        igl::predicates
        polyscope::polyscope
        TinyAD::TinyAD
)
