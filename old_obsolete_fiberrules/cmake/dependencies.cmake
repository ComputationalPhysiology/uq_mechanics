# CMake help file to add dependencies of fiberrules

#------------------------------------------------------------------------------
# Run tests to find required packages

find_package(PythonInterp)
if (FIBERRULES_ENABLE_PYTHON)
  find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT)

  # If Python is found, check for NumPy and SWIG
  if (PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND)
    find_package(NumPy REQUIRED)

    find_package(SWIG REQUIRED)
    if (${SWIG_VERSION} LESS 2.0)
      message(FATAL_ERROR " FIBERRULES requires SWIG version 2.0 or greater. You have version ${SWIG_VERSION}. Set FIBERRULES_ENABLE_PYTHON to False or install correct SWIG version.")
    endif()
    include(UseSWIG)
    set(PYTHON_FOUND TRUE)

    # Set numpy version define
    set(FIBERRULES_PYTHON_DEFINITIONS -DNUMPY_VERSION_MAJOR=${NUMPY_VERSION_MAJOR} -DNUMPY_VERSION_MINOR=${NUMPY_VERSION_MINOR} -DNUMPY_VERSION_MICRO=${NUMPY_VERSION_MICRO})

    # Only set define for none depricated API for NUMPY version 1.7 and larger
    if(NUMPY_VERSION VERSION_GREATER 1.6.2)
      set(FIBERRULES_PYTHON_DEFINITIONS ${FIBERRULES_PYTHON_DEFINITIONS} -DNPY_NO_DEPRECATED_API=NPY_${NUMPY_VERSION_MAJOR}_${NUMPY_VERSION_MINOR}_API_VERSION )
    endif()

  endif()
endif()

# Look for Boost (scoped_array support)
#set(Boost_ADDITIONAL_VERSIONS 1.43 1.43.0 1.44 1.44.0 1.45 1.45.0 1.46 1.46.0 
#    1.46.1 1.47 1.47.0 1.48 1.48.0 1.49 1.49.0 1.50 1.50.0)
#set(BOOST_ROOT $ENV{BOOST_DIR})
#find_package(Boost 1.36 QUIET)

# Check for OpenMP
if (FIBERRULES_ENABLE_OPENMP)
  find_package(OpenMP)
  include(CheckOpenMP)
  check_openmp_unsigned_int_loop_control_variable(OPENMP_UINT_TEST_RUNS)
  if (NOT OPENMP_UINT_TEST_RUNS)
    set(OPENMP_FOUND FALSE)
  endif()
endif()

# Check for LAPACK
if (FIBERRULES_ENABLE_LAPACK)
  find_package(LAPACK)
  find_package(LAPACKHeader)
endif()

# Check for ARMADILLO
if (FIBERRULES_ENABLE_ARMADILLO)
  find_package(Armadillo)
endif()

#------------------------------------------------------------------------------
# Get installation path for Python modules

# Get Python module path from distutils
if (PYTHONINTERP_FOUND)

  if (NOT DEFINED FIBERRULES_INSTALL_PYTHON_EXT_DIR)
    # Get path for platform-dependent Python modules (since we install a binary libary)
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
      OUTPUT_VARIABLE FIBERRULES_INSTALL_PYTHON_EXT_DIR
      )
    # Strip off CMAKE_INSTALL_PREFIX (is added later by CMake)
    string(REGEX REPLACE "${CMAKE_INSTALL_PREFIX}(/|\\\\)([^ ]*)" "\\2"
      FIBERRULES_INSTALL_PYTHON_EXT_DIR "${FIBERRULES_INSTALL_PYTHON_EXT_DIR}")
    set(FIBERRULES_INSTALL_PYTHON_EXT_DIR ${FIBERRULES_INSTALL_PYTHON_EXT_DIR}
      CACHE PATH "Python extension module installation directory.")
  endif()

  if (NOT DEFINED FIBERRULES_INSTALL_PYTHON_MODULE_DIR)
    # Get path for pure Python modules
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=False, prefix='${CMAKE_INSTALL_PREFIX}'))"
      OUTPUT_VARIABLE FIBERRULES_INSTALL_PYTHON_MODULE_DIR
      )
    # Strip off CMAKE_INSTALL_PREFIX (is added later by CMake)
    string(REGEX REPLACE "${CMAKE_INSTALL_PREFIX}(/|\\\\)([^ ]*)" "\\2"
      FIBERRULES_INSTALL_PYTHON_MODULE_DIR "${FIBERRULES_INSTALL_PYTHON_MODULE_DIR}")
    set(FIBERRULES_INSTALL_PYTHON_MODULE_DIR ${FIBERRULES_INSTALL_PYTHON_MODULE_DIR}
      CACHE PATH "Python module installation directory.")
  endif()
endif (PYTHONINTERP_FOUND)

#------------------------------------------------------------------------------
# Add include directories and libs of required packages
# Boost
#list(APPEND FIBERRULES_DEP_SYSTEM_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
#
# OpenMP
if (FIBERRULES_ENABLE_OPENMP AND OPENMP_FOUND)
  list(APPEND FIBERRULES_CXX_DEFINITIONS "-DHAS_OPENMP")
  set(FIBERRULES_CXX_FLAGS "${FIBERRULES_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

  if (MINGW)
    list(APPEND FIBERRULES_LINK_FLAGS "${OpenMP_CXX_FLAGS}")
  endif()

endif()

# LAPACK
if (FIBERRULES_ENABLE_LAPACK AND LAPACKHEADER_FOUND AND LAPACK_FOUND)

  list(APPEND FIBERRULES_DEP_INCLUDE_DIRECTORIES "${LAPACK_INCLUDE_DIRS}")
  list(APPEND FIBERRULES_CXX_DEFINITIONS "-DHAS_LAPACK")

  list(APPEND FIBERRULES_TARGET_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
  list(APPEND FIBERRULES_LINK_FLAGS "${LAPACK_LINKER_FLAGS}")
endif()

# Armadillo
if (FIBERRULES_ENABLE_ARMADILLO AND ARMADILLO_FOUND)

  list(APPEND FIBERRULES_DEP_INCLUDE_DIRECTORIES "${ARMADILLO_INCLUDE_DIRS}")
  list(APPEND FIBERRULES_CXX_DEFINITIONS "-DHAS_ARMADILLO")

  list(APPEND FIBERRULES_TARGET_LINK_LIBRARIES "${ARMADILLO_LIBRARIES}")
  list(APPEND FIBERRULES_LINK_FLAGS "${ARMADILLO_LINK_FLAGS}")
endif()

#------------------------------------------------------------------------------
# Set compiler flags, include directories and library dependencies

# Add compiler include directories
include_directories(${FIBERRULES_SOURCE_DIR} ${FIBERRULES_DEP_INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${FIBERRULES_DEP_SYSTEM_INCLUDE_DIRECTORIES})

# Add CXX defintions
add_definitions(${FIBERRULES_CXX_DEFINITIONS})
add_definitions(-DFIBERRULES_VERSION="${FIBERRULES_VERSION}")

# Add flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FIBERRULES_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${FIBERRULES_LINK_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${FIBERRULES_LINK_FLAGS}")

#------------------------------------------------------------------------------
# If everything needed for the Python extension module is found

if (PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND AND SWIG_FOUND) # AND Boost_FOUND
    set(ENABLE_PYTHON_EXTENSION_MODULE TRUE)
endif()
