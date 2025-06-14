cmake_minimum_required(VERSION 3.10)
project(Sparse_Interp)


################################################################################

# libigl
# option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
# option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"              OFF)
# option(LIBIGL_WITH_CGAL              "Use CGAL"                     OFF)
# option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
# option(LIBIGL_WITH_CORK              "Use Cork"                     OFF)
# option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
# option(LIBIGL_WITH_LIM               "Use LIM"                      OFF)
# option(LIBIGL_WITH_MATLAB            "Use Matlab"                   OFF)
# option(LIBIGL_WITH_MOSEK             "Use MOSEK"                    OFF)
# option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   OFF)
# option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     OFF)
# option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    OFF)
# option(LIBIGL_WITH_PNG               "Use PNG"                      OFF)
# option(LIBIGL_WITH_PYTHON            "Use Python"                   OFF)
# option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
# option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 ON)
# option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"            OFF)
# option(LIBIGL_WITH_XML               "Use XML"                      OFF)
# option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         ON)

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
#igl_include(opengl)
#igl_include(glfw)
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
#target_sources(splines PRIVATE ${SPLINE_SOURCES})
target_include_directories(splines PUBLIC ${PROJECT_BINARY_DIR}/include)
target_link_libraries(splines PUBLIC
    igl::core
    #igl::predicates
    igl_restricted::triangle 
    Eigen3::Eigen 
    OpenMeshCore 
    OpenMeshTools
)
target_compile_definitions(splines PUBLIC
    SI_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/app/meshes_new/")
#target_compile_options(splines PRIVATE -g)
###############################

# the comparison methods
include(comparisons/compare.cmake)
add_library(compare  ${COMPARE_SOURCES})
SI_copy_headers(${COMPARE_SOURCES})
target_include_directories(compare PUBLIC ${PROJECT_BINARY_DIR}/include)
target_link_libraries(compare PUBLIC splines)
#target_compile_options(compare PRIVATE -g)

################################

# add_executable(${PROJECT_NAME}_bin
#     app/main.cpp 
# 	app/test.cpp
#   )
# target_include_directories(${PROJECT_NAME}_bin PUBLIC ${PROJECT_BINARY_DIR}/include)
# target_link_libraries(${PROJECT_NAME}_bin 
#     PUBLIC
#         splines 
#         igl::core 
# )
# target_link_libraries(${PROJECT_NAME}_bin PUBLIC compare)
# target_compile_options(${PROJECT_NAME}_bin PRIVATE -g)

# add_executable(example_bin
#     app/example.cpp
#   )
# target_include_directories(example_bin PUBLIC ${PROJECT_BINARY_DIR}/include)
# target_link_libraries(example_bin 
#     PUBLIC 
#         splines 
#         igl::core
# )

add_executable(surfGui_bin
    src/main.cpp
    src/methods.cpp
)

target_include_directories(surfGui_bin PUBLIC 
    ${PROJECT_BINARY_DIR}/include
    ${PROJECT_SOURCE_DIR}/app
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
#target_compile_options(surfGui_bin PRIVATE -g)
