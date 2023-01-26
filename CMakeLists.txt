# Required CMake Version
###############
cmake_minimum_required(VERSION 2.8)
# Project setup
###############
PROJECT(Schroedinger)
set(binname viennaio)
set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)
set(VERSION_PATCH 0)
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})
#ADD_DEFINITIONS(-DVIENNAMC_VERSION="${VERSION}")
# Default build type and set options for ccmake
##############
SET(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
# Default compile options
#################
# SET(CMAKE_CXX_FLAGS "-Wall -pedantic")
SET(CMAKE_CXX_FLAGS "-Wall")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 ")
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-omit-frame-pointer -fsanitize=address")
SET(CMAKE_CXX_STANDARD 14)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fmax-errors=1")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ferror-limit=1")
endif()

# Specify general include directory
################
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/ )
# Specify that the output executable is placed in the root build directory
################
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
# Specify that the libraries are placed in the root build directory
################
set(LIBRARY_OUTPUT_PATH    ${CMAKE_SOURCE_DIR}/lib/Schroedinger)
LINK_DIRECTORIES(${LIBRARY_OUTPUT_PATH})



#SET(LIBRARIES ${PYTHON_LIBRARIES})
# Subdirectories
################
#add_subdirectory(src)
# Add Sources
################
#aux_source_directory("./src" SRC_BASE)
#aux_source_directory("./src/Test" SRC_BASE)
# add executable and link
#ADD_EXECUTABLE(${binname}-${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} ${SRC_BASE})
#IF(VIENNAMC_STATIC_BUILD)
#  TARGET_LINK_LIBRARIES(${binname}-${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} ${LIBRARIES}  "-static")
#ELSE()
#  TARGET_LINK_LIBRARIES(${binname}-${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} ${LIBRARIES} )
#  ENDIF(VIENNAMC_STATIC_BUILD)
# Find prerequisites
####################
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(pybind11/include/pybind11)
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/ )
add_subdirectory(pybind11)



find_package(OpenMP REQUIRED)
add_executable(parallelSchroedinger parallelSchroedinger.cpp)
target_link_libraries(parallelSchroedinger PRIVATE OpenMP::OpenMP_CXX)
target_include_directories(parallelSchroedinger PRIVATE  ${CMAKE_SOURCE_DIR}/src/)

add_executable(serialSchroedinger serialSchroedinger.cpp)
target_link_libraries(serialSchroedinger PRIVATE OpenMP::OpenMP_CXX)
target_include_directories(serialSchroedinger PRIVATE  ${CMAKE_SOURCE_DIR}/src/)











set(SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pybind11/pyIO.cpp
    )
pybind11_add_module(Schroedinger ${SRC})
target_link_libraries(Schroedinger PRIVATE OpenMP::OpenMP_CXX)
