# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/rocio/Desktop/fiberrules

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/rocio/Desktop/fiberrules/build

# Include any dependencies generated for this target.
include fiberrules/CMakeFiles/fiberrules.dir/depend.make

# Include the progress variables for this target.
include fiberrules/CMakeFiles/fiberrules.dir/progress.make

# Include the compile flags for this target's objects.
include fiberrules/CMakeFiles/fiberrules.dir/flags.make

fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o: fiberrules/CMakeFiles/fiberrules.dir/flags.make
fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o: ../fiberrules/carptools.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/rocio/Desktop/fiberrules/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o"
	cd /home/rocio/Desktop/fiberrules/build/fiberrules && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/fiberrules.dir/carptools.cpp.o -c /home/rocio/Desktop/fiberrules/fiberrules/carptools.cpp

fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fiberrules.dir/carptools.cpp.i"
	cd /home/rocio/Desktop/fiberrules/build/fiberrules && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/rocio/Desktop/fiberrules/fiberrules/carptools.cpp > CMakeFiles/fiberrules.dir/carptools.cpp.i

fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fiberrules.dir/carptools.cpp.s"
	cd /home/rocio/Desktop/fiberrules/build/fiberrules && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/rocio/Desktop/fiberrules/fiberrules/carptools.cpp -o CMakeFiles/fiberrules.dir/carptools.cpp.s

fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o.requires:
.PHONY : fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o.requires

fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o.provides: fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o.requires
	$(MAKE) -f fiberrules/CMakeFiles/fiberrules.dir/build.make fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o.provides.build
.PHONY : fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o.provides

fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o.provides.build: fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o

fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o: fiberrules/CMakeFiles/fiberrules.dir/flags.make
fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o: ../fiberrules/fiberrulestools.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/rocio/Desktop/fiberrules/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o"
	cd /home/rocio/Desktop/fiberrules/build/fiberrules && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o -c /home/rocio/Desktop/fiberrules/fiberrules/fiberrulestools.cpp

fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fiberrules.dir/fiberrulestools.cpp.i"
	cd /home/rocio/Desktop/fiberrules/build/fiberrules && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/rocio/Desktop/fiberrules/fiberrules/fiberrulestools.cpp > CMakeFiles/fiberrules.dir/fiberrulestools.cpp.i

fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fiberrules.dir/fiberrulestools.cpp.s"
	cd /home/rocio/Desktop/fiberrules/build/fiberrules && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/rocio/Desktop/fiberrules/fiberrules/fiberrulestools.cpp -o CMakeFiles/fiberrules.dir/fiberrulestools.cpp.s

fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o.requires:
.PHONY : fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o.requires

fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o.provides: fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o.requires
	$(MAKE) -f fiberrules/CMakeFiles/fiberrules.dir/build.make fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o.provides.build
.PHONY : fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o.provides

fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o.provides.build: fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o

# Object files for target fiberrules
fiberrules_OBJECTS = \
"CMakeFiles/fiberrules.dir/carptools.cpp.o" \
"CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o"

# External object files for target fiberrules
fiberrules_EXTERNAL_OBJECTS =

fiberrules/libfiberrules.so.0.1.0: fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o
fiberrules/libfiberrules.so.0.1.0: fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o
fiberrules/libfiberrules.so.0.1.0: fiberrules/CMakeFiles/fiberrules.dir/build.make
fiberrules/libfiberrules.so.0.1.0: fiberrules/CMakeFiles/fiberrules.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library libfiberrules.so"
	cd /home/rocio/Desktop/fiberrules/build/fiberrules && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fiberrules.dir/link.txt --verbose=$(VERBOSE)
	cd /home/rocio/Desktop/fiberrules/build/fiberrules && $(CMAKE_COMMAND) -E cmake_symlink_library libfiberrules.so.0.1.0 libfiberrules.so.0.1 libfiberrules.so

fiberrules/libfiberrules.so.0.1: fiberrules/libfiberrules.so.0.1.0

fiberrules/libfiberrules.so: fiberrules/libfiberrules.so.0.1.0

# Rule to build all files generated by this target.
fiberrules/CMakeFiles/fiberrules.dir/build: fiberrules/libfiberrules.so
.PHONY : fiberrules/CMakeFiles/fiberrules.dir/build

fiberrules/CMakeFiles/fiberrules.dir/requires: fiberrules/CMakeFiles/fiberrules.dir/carptools.cpp.o.requires
fiberrules/CMakeFiles/fiberrules.dir/requires: fiberrules/CMakeFiles/fiberrules.dir/fiberrulestools.cpp.o.requires
.PHONY : fiberrules/CMakeFiles/fiberrules.dir/requires

fiberrules/CMakeFiles/fiberrules.dir/clean:
	cd /home/rocio/Desktop/fiberrules/build/fiberrules && $(CMAKE_COMMAND) -P CMakeFiles/fiberrules.dir/cmake_clean.cmake
.PHONY : fiberrules/CMakeFiles/fiberrules.dir/clean

fiberrules/CMakeFiles/fiberrules.dir/depend:
	cd /home/rocio/Desktop/fiberrules/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rocio/Desktop/fiberrules /home/rocio/Desktop/fiberrules/fiberrules /home/rocio/Desktop/fiberrules/build /home/rocio/Desktop/fiberrules/build/fiberrules /home/rocio/Desktop/fiberrules/build/fiberrules/CMakeFiles/fiberrules.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : fiberrules/CMakeFiles/fiberrules.dir/depend
