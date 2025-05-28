# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.


# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(${PROJECT_NAME}DownloadExternal)
#include(polyscope_patcher)
include(FetchContent)

################################################################################
# Required libraries
################################################################################

# Eigen
if(NOT TARGET Eigen3::Eigen)
  sparse_interp_download_eigen()
  add_library(${PROJECT_NAME}_eigen INTERFACE)
  target_include_directories(${PROJECT_NAME}_eigen SYSTEM INTERFACE
    $<BUILD_INTERFACE:${SPARSE_EXTERNAL}/eigen>
    $<INSTALL_INTERFACE:include>
  )
  set_property(TARGET ${PROJECT_NAME}_eigen PROPERTY EXPORT_NAME Eigen3::Eigen)
  add_library(Eigen3::Eigen ALIAS ${PROJECT_NAME}_eigen)
  # Set Eigen directory environment variable (needed for EVCTCD)
  set(ENV{EIGEN3_INCLUDE_DIR} "${SPARSE_EXTERNAL}/eigen/")
endif()


  # libigl for timing
# if(NOT TARGET igl::core)
#   # sparse_interp_download_libigl()

#   #   # Import libigl targets
#   #   list(APPEND CMAKE_MODULE_PATH "${SPARSE_EXTERNAL}/libigl/cmake")
#   #   include(libigl)
  
#   endif()
if(NOT TARGET igl::core)
    FetchContent_Declare(
        libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG v2.6.0
    )
endif()
FetchContent_MakeAvailable(libigl)

# Polyscope
if(NOT TARGET polyscope::polyscope)
    FetchContent_Declare(
        polyscope
        GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
        GIT_TAG v2.4.0
    )
endif()
FetchContent_MakeAvailable(polyscope)
# Patch target clash with libigl
#FetchContent_GetProperties(polyscope)
#if (NOT polyscope_POPULATED)
  #FetchContent_Populate(polyscope)
  #patch_poly(glad polyscope_glad "${polyscope_SOURCE_DIR}/deps/glad/src")
  #patch_poly(glfw polyscope_glfw "${polyscope_SOURCE_DIR}/deps/glfw/src")
#endif()
#add_subdirectory(${polyscope_SOURCE_DIR} ${polyscope_BINARY_DIR})


# install openmesh
set(OM_FILE "${CMAKE_CURRENT_SOURCE_DIR}/external/openmesh.zip")
set(OM_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/openmesh" )
if(EXISTS "${OM_FILE}" OR EXISTS "${OM_PATH}")
  message("LSC: OpenMesh source file exists")
else()
  if(NOT EXISTS "${OM_FILE}")
    message("LSC: downloading OpenMesh")
    file(DOWNLOAD https://www.graphics.rwth-aachen.de/media/openmesh_static/Releases/11.0/OpenMesh-11.0.0.zip ${OM_FILE})
  endif()
endif()
if(NOT EXISTS "${OM_PATH}")
    message("LSC: OpenMesh unzipping")
    file(ARCHIVE_EXTRACT INPUT ${OM_FILE} DESTINATION ${OM_PATH})
    message("LSC: OpenMesh is unzipped")
else()
    message("LSC: OpenMesh is already unzipped")
endif()
message("LSC: OpenMesh file: ${OM_PATH}/OpenMesh-11.0.0/")

