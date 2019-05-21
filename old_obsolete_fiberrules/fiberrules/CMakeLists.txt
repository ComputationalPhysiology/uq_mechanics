set(FIBERRULES_H fiberrules.h)

#------------------------------------------------------------------------------
# FIBERRULES source directories

# Headers in this directory
file(GLOB HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)

# Sources in this directory
file(GLOB SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

#------------------------------------------------------------------------------
# Install header files

install(FILES ${FIBERRULES_H} DESTINATION ${FIBERRULES_INCLUDE_DIR} COMPONENT Development)

foreach(_HEADER ${HEADERS})
  list(FIND ${FIBERRULES_H} ${_HEADER} INDEX_EXCLUDES)
  if (${INDEX_EXCLUDES} LESS 0)
    install(FILES ${_HEADER} DESTINATION ${FIBERRULES_INCLUDE_DIR}/fiberrules COMPONENT Development)
  endif()
endforeach()

#------------------------------------------------------------------------------
# Define libraries
add_library(fiberrules ${FIBERRULES_H} ${HEADERS} ${SOURCES})
set_target_properties(fiberrules PROPERTIES ${FIBERRULES_LIBRARY_PROPERTIES})

# Strip leading and trailing whitespaces
string(STRIP FIBERRULES_TARGET_LINK_LIBRARIES "${FIBERRULES_TARGET_LINK_LIBRARIES}")

# Add FIBERRULES target libraries
target_link_libraries(fiberrules ${FIBERRULES_TARGET_LINK_LIBRARIES})

#------------------------------------------------------------------------------
# SWIG

if (ENABLE_PYTHON_EXTENSION_MODULE)

  # Add SWIG subdirectory
  add_subdirectory(swig)

endif()

#------------------------------------------------------------------------------
# Install

install(TARGETS fiberrules
  RUNTIME DESTINATION ${FIBERRULES_LIB_DIR} COMPONENT RuntimeExecutables
  LIBRARY DESTINATION ${FIBERRULES_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${FIBERRULES_LIB_DIR} COMPONENT Development
  )
