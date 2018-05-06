project(swig)

set(SWIG_MODULE_NAME cpp)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Remove '-Werror', '-Wall' and 'pedantic' flags (if present) when compiling SWIG-generated files
string(REGEX REPLACE "-Wall" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REGEX REPLACE "-Wall" " " CMAKE_CXX_FLAGS_DEVELOPER "${CMAKE_CXX_FLAGS_DEVELOPER}")
string(REGEX REPLACE "-Wshadow" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REGEX REPLACE "-Wshadow" " " CMAKE_CXX_FLAGS_DEVELOPER "${CMAKE_CXX_FLAGS_DEVELOPER}")
string(REGEX REPLACE "-Werror=format-security" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REGEX REPLACE "-Werror=format-security" " " CMAKE_CXX_FLAGS_DEVELOPER "${CMAKE_CXX_FLAGS_DEVELOPER}")
string(REGEX REPLACE "-Werror" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REGEX REPLACE "-Werror" " " CMAKE_CXX_FLAGS_DEVELOPER "${CMAKE_CXX_FLAGS_DEVELOPER}")
string(REGEX REPLACE "-pedantic" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REGEX REPLACE "-pedantic" " " CMAKE_CXX_FLAGS_DEVELOPER "${CMAKE_CXX_FLAGS_DEVELOPER}")

# Set common include directoris
include_directories(${FIBERRULES_SOURCE_DIR} ${FIBERRULES_SOURCE_DIR}/fiberrules/swig ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIR})

# Set SWIG flags 
set(CMAKE_SWIG_FLAGS
  -module ${SWIG_MODULE_NAME}
  -O
#  -shadow
#  -modern
#  -modernargs
#  -fastdispatch
#  -fvirtual
#  -nosafecstrings
#  -noproxydel
#  -fastproxy
#  -fastinit
#  -fastunpack
#  -fastquery
#  -nobuildnone
  -Iinclude/swig
  ${FIBERRULES_CXX_DEFINITIONS}
  ${FIBERRULES_PYTHON_DEFINITIONS}
  )

add_definitions(${FIBERRULES_PYTHON_DEFINITIONS})

# Get all SWIG interface files

# All header files are dependencies
file(GLOB FIBERRULES_HEADERS ${CMAKE_CURRENT_SOURCE_DIR} ../*.h)
file(GLOB GLOBAL_INTERFACE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.i)

list(APPEND FIBERRULES_SWIG_DEPENDENCIES ${FIBERRULES_HEADERS})
list(APPEND FIBERRULES_SWIG_DEPENDENCIES ${GLOBAL_INTERFACE_FILES})

# This prevents swig being run unnecessarily
set_source_files_properties(fiberrules.i PROPERTIES SWIG_MODULE_NAME ${SWIG_MODULE_NAME})

# Tell CMake SWIG has generated a C++ file
set_source_files_properties(fiberrules.i PROPERTIES CPLUSPLUS ON)

# Generate SWIG files in 
set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR})

# Tell CMake which SWIG interface files should be checked for changes when recompile
set(SWIG_MODULE_${SWIG_MODULE_NAME}_EXTRA_DEPS ${FIBERRULES_SWIG_DEPENDENCIES})

# Tell CMake to run SWIG on fiberrules.i and to link against libfiberrules
swig_add_module(${SWIG_MODULE_NAME} python fiberrules.i)
swig_link_libraries(${SWIG_MODULE_NAME} fiberrules ${PYTHON_LIBRARIES})

# Install Python .py files
get_target_property(SWIG_MODULE_LOCATION ${SWIG_MODULE_${SWIG_MODULE_NAME}_REAL_NAME} LOCATION)

install(FILES
  ${SWIG_MODULE_LOCATION} ${CMAKE_CURRENT_BINARY_DIR}/${SWIG_MODULE_NAME}.py
  DESTINATION ${FIBERRULES_INSTALL_PYTHON_MODULE_DIR}/fiberrules
  COMPONENT RuntimeLibraries
  )

# Install all SWIG interface files
install(FILES ${GLOBAL_INTERFACE_FILES}
  DESTINATION 
  ${FIBERRULES_INCLUDE_DIR}/fiberrules/swig
  COMPONENT Development)

