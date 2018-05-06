# Top level CMakeLists.txt file for FIBERRULES

# Require CMake 2.6
cmake_minimum_required(VERSION 2.6)

#------------------------------------------------------------------------------
# Set project name and version number

project(FIBERRULES)

set(FIBERRULES_VERSION_MAJOR "0")
set(FIBERRULES_VERSION_MINOR "1")
set(FIBERRULES_VERSION_MICRO "0")
set(FIBERRULES_VERSION "${FIBERRULES_VERSION_MAJOR}.${FIBERRULES_VERSION_MINOR}.${FIBERRULES_VERSION_MICRO}")
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER) 

# Set special link option, see `cmake --help-policy CMP0003`
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif()

# Set verbose output while testing CMake
set(CMAKE_VERBOSE_MAKEFILE 1)

# Set location of our FindFoo.cmake modules
set(FIBERRULES_CMAKE_DIR "${FIBERRULES_SOURCE_DIR}/cmake" CACHE INTERNAL "")
set(CMAKE_MODULE_PATH "${FIBERRULES_CMAKE_DIR}/modules")

#------------------------------------------------------------------------------
# Options

option(FIBERRULES_ENABLE_PYTHON "Enable Python extensions." ON)
option(FIBERRULES_ENABLE_OPENMP "Enable OpenMP." ON)
option(FIBERRULES_ENABLE_ARMADILLO "Enable ARMADILLO." ON)
option(FIBERRULES_WITH_LIBRARY_VERSION "Build with library version information." ON)
option(FIBERRULES_ENABLE_UNIT_TESTS "Build unit tests." ON)
option(BUILD_SHARED_LIBS "Build FIBERRULES with shared libraries." ON)
option(CMAKE_SKIP_RPATH "Do not add runtime paths when using shared libraries." ON)
option(CMAKE_INSTALL_RPATH_USE_LINK_PATH "Add paths to linker search and installed rpath." OFF)
option(CMAKE_USE_RELATIVE_PATHS "Use relative paths in makefiles and projects." OFF)

#------------------------------------------------------------------------------
# Compiler flags

# Add ansi (c++98)
if (CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi")
endif()

# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug Developer MinSizeRel Release RelWithDebInfo." FORCE)
endif()

# Check for some compiler flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-pipe HAS_PIPE)
if (HAS_PIPE)
  set(FIBERRULES_CXX_DEVELOPER_FLAGS "-pipe ${FIBERRULES_CXX_DEVELOPER_FLAGS}")
endif()

# Add some strict compiler checks
CHECK_CXX_COMPILER_FLAG("-Wall -Werror -pedantic" HAS_PEDANTIC)
if (HAS_PEDANTIC)
  set(FIBERRULES_CXX_DEVELOPER_FLAGS "-Wall -Werror -pedantic ${FIBERRULES_CXX_DEVELOPER_FLAGS}")
endif()

CHECK_CXX_COMPILER_FLAG("-Wshadow -Wextra" HAS_EXTRA)
if (HAS_EXTRA)
  # FIXME: Remove shadow. Is not liked by ARMADILLO
  set(FIBERRULES_CXX_DEVELOPER_FLAGS "-Wextra ${FIBERRULES_CXX_DEVELOPER_FLAGS}")
endif()

# Debug flags
CHECK_CXX_COMPILER_FLAG(-g HAS_DEBUG)
if (HAS_DEBUG)
  set(FIBERRULES_CXX_DEVELOPER_FLAGS "-g ${FIBERRULES_CXX_DEVELOPER_FLAGS}")
endif()

CHECK_CXX_COMPILER_FLAG(-O2 HAS_O2_OPTIMISATION)
if (HAS_O2_OPTIMISATION)
  set(FIBERRULES_CXX_DEVELOPER_FLAGS "-O2 ${FIBERRULES_CXX_DEVELOPER_FLAGS}")
endif()

# Set 'Developer' build type flags
set(CMAKE_CXX_FLAGS_DEVELOPER "${FIBERRULES_CXX_DEVELOPER_FLAGS}" CACHE STRING
  "Flags used by the compiler during development." FORCE)

# FIXME: Do we want to add -DDEBUG to RelWithDebInfo?
# Add debug definitions
if (CMAKE_BUILD_TYPE STREQUAL "Developer" OR CMAKE_BUILD_TYPE STREQUAL "Debug")
  list(APPEND FIBERRULES_CXX_DEFINITIONS "-DDEBUG")
endif()

# Disable some warnings for Intel compilers
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd161,488,654,858,1125,1170")
endif()

# Set system include flags to get around CMake bug on OSX
# See http://public.kitware.com/Bug/print_bug_page.php?bug_id=10837
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
endif()

#------------------------------------------------------------------------------
# Installation of FIBERRULES library

# Set FIBERRULES install sub-directories
set(FIBERRULES_BIN_DIR "bin" CACHE PATH "Binary installation directory.")
set(FIBERRULES_LIB_DIR "lib" CACHE PATH "Library installation directory.")
set(FIBERRULES_INCLUDE_DIR "include" CACHE PATH "C/C++ header installation directory.")
set(FIBERRULES_PKGCONFIG_DIR "lib/pkgconfig" CACHE PATH "pkg-config file installation directory.")
set(FIBERRULES_SHARE_DIR "share/fiberrules" CACHE PATH "Shared data installation directory.")
#set(FIBERRULES_MAN_DIR "share/man" CACHE PATH "Manual page installation directory.")
#set(FIBERRULES_DOC_DIR "${FIBERRULES_SHARE_DIR}/doc" CACHE PATH "FIBERRULES Documentation directory.")

# Append the library version information to the library target properties
if (FIBERRULES_WITH_LIBRARY_VERSION)
  string(REPLACE "+" "" FIBERRULES_LIBRARY_VERSION ${FIBERRULES_VERSION})
  # This setting of SOVERSION assumes that any API change
  # will increment either the minor or major version number.
  set(FIBERRULES_LIBRARY_PROPERTIES ${FIBERRULES_LIBRARY_PROPERTIES}
    VERSION ${FIBERRULES_LIBRARY_VERSION}
    SOVERSION ${FIBERRULES_VERSION_MAJOR}.${FIBERRULES_VERSION_MINOR}
  )
endif()

# Check dependencies
include(cmake/dependencies.cmake)

# Add source directory
add_subdirectory(fiberrules)

if (ENABLE_PYTHON_EXTENSION_MODULE)
  install(DIRECTORY
    ${CMAKE_SOURCE_DIR}/site-packages/fiberrules
    DESTINATION ${FIBERRULES_INSTALL_PYTHON_MODULE_DIR}
    USE_SOURCE_PERMISSIONS
    COMPONENT RuntimeLibraries
    )
endif()

#------------------------------------------------------------------------------
# If enabling unit testing add test subdirectory

if (FIBERRULES_ENABLE_UNIT_TESTS)
  enable_testing()
#  add_subdirectory(test)
  
endif()

#------------------------------------------------------------------------------
# Add "make uninstall" target

configure_file(
  "${FIBERRULES_CMAKE_DIR}/templates/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

# Generate config files
#include(cmake/generate_config_files.cmake)

