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
include fiberrules/swig/CMakeFiles/_cpp.dir/depend.make

# Include the progress variables for this target.
include fiberrules/swig/CMakeFiles/_cpp.dir/progress.make

# Include the compile flags for this target's objects.
include fiberrules/swig/CMakeFiles/_cpp.dir/flags.make

fiberrules/swig/fiberrulesPYTHON_wrap.cxx: ../fiberrules/swig
fiberrules/swig/fiberrulesPYTHON_wrap.cxx: ../fiberrules/swig/../fiberrules.h
fiberrules/swig/fiberrulesPYTHON_wrap.cxx: ../fiberrules/swig/../types.h
fiberrules/swig/fiberrulesPYTHON_wrap.cxx: ../fiberrules/swig/../fiberrulestools.h
fiberrules/swig/fiberrulesPYTHON_wrap.cxx: ../fiberrules/swig/../carptools.h
fiberrules/swig/fiberrulesPYTHON_wrap.cxx: ../fiberrules/swig/pre.i
fiberrules/swig/fiberrulesPYTHON_wrap.cxx: ../fiberrules/swig/version.i
fiberrules/swig/fiberrulesPYTHON_wrap.cxx: ../fiberrules/swig/exceptions.i
fiberrules/swig/fiberrulesPYTHON_wrap.cxx: ../fiberrules/swig/typemaps.i
fiberrules/swig/fiberrulesPYTHON_wrap.cxx: ../fiberrules/swig/post.i
fiberrules/swig/fiberrulesPYTHON_wrap.cxx: ../fiberrules/swig/fiberrules.i
fiberrules/swig/fiberrulesPYTHON_wrap.cxx: ../fiberrules/swig/fiberrules.i
	$(CMAKE_COMMAND) -E cmake_progress_report /home/rocio/Desktop/fiberrules/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Swig source"
	cd /home/rocio/Desktop/fiberrules/build/fiberrules/swig && /usr/bin/cmake -E make_directory /home/rocio/Desktop/fiberrules/build/fiberrules/swig
	cd /home/rocio/Desktop/fiberrules/build/fiberrules/swig && /usr/bin/swig2.0 -python -module cpp -O -Iinclude/swig -DHAS_OPENMP -DNUMPY_VERSION_MAJOR=1 -DNUMPY_VERSION_MINOR=13 -DNUMPY_VERSION_MICRO=3 -DNPY_NO_DEPRECATED_API=NPY_1_13_API_VERSION -outdir /home/rocio/Desktop/fiberrules/build/fiberrules/swig -c++ -I/home/rocio/Desktop/fiberrules -I/home/rocio/Desktop/fiberrules -I/home/rocio/Desktop/fiberrules/fiberrules/swig -I/usr/include/python2.7 -I/usr/include/i386-linux-gnu/python2.7 -I/usr/local/lib/python2.7/dist-packages/numpy/core/include -o /home/rocio/Desktop/fiberrules/build/fiberrules/swig/fiberrulesPYTHON_wrap.cxx /home/rocio/Desktop/fiberrules/fiberrules/swig/fiberrules.i

fiberrules/swig/cpp.py: fiberrules/swig/fiberrulesPYTHON_wrap.cxx

fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o: fiberrules/swig/CMakeFiles/_cpp.dir/flags.make
fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o: fiberrules/swig/fiberrulesPYTHON_wrap.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /home/rocio/Desktop/fiberrules/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o"
	cd /home/rocio/Desktop/fiberrules/build/fiberrules/swig && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o -c /home/rocio/Desktop/fiberrules/build/fiberrules/swig/fiberrulesPYTHON_wrap.cxx

fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.i"
	cd /home/rocio/Desktop/fiberrules/build/fiberrules/swig && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/rocio/Desktop/fiberrules/build/fiberrules/swig/fiberrulesPYTHON_wrap.cxx > CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.i

fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.s"
	cd /home/rocio/Desktop/fiberrules/build/fiberrules/swig && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/rocio/Desktop/fiberrules/build/fiberrules/swig/fiberrulesPYTHON_wrap.cxx -o CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.s

fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o.requires:
.PHONY : fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o.requires

fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o.provides: fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o.requires
	$(MAKE) -f fiberrules/swig/CMakeFiles/_cpp.dir/build.make fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o.provides.build
.PHONY : fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o.provides

fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o.provides.build: fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o

# Object files for target _cpp
_cpp_OBJECTS = \
"CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o"

# External object files for target _cpp
_cpp_EXTERNAL_OBJECTS =

fiberrules/swig/_cpp.so: fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o
fiberrules/swig/_cpp.so: fiberrules/swig/CMakeFiles/_cpp.dir/build.make
fiberrules/swig/_cpp.so: fiberrules/libfiberrules.so.0.1.0
fiberrules/swig/_cpp.so: /usr/lib/i386-linux-gnu/libpython2.7.so
fiberrules/swig/_cpp.so: fiberrules/swig/CMakeFiles/_cpp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared module _cpp.so"
	cd /home/rocio/Desktop/fiberrules/build/fiberrules/swig && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/_cpp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
fiberrules/swig/CMakeFiles/_cpp.dir/build: fiberrules/swig/_cpp.so
.PHONY : fiberrules/swig/CMakeFiles/_cpp.dir/build

fiberrules/swig/CMakeFiles/_cpp.dir/requires: fiberrules/swig/CMakeFiles/_cpp.dir/fiberrulesPYTHON_wrap.cxx.o.requires
.PHONY : fiberrules/swig/CMakeFiles/_cpp.dir/requires

fiberrules/swig/CMakeFiles/_cpp.dir/clean:
	cd /home/rocio/Desktop/fiberrules/build/fiberrules/swig && $(CMAKE_COMMAND) -P CMakeFiles/_cpp.dir/cmake_clean.cmake
.PHONY : fiberrules/swig/CMakeFiles/_cpp.dir/clean

fiberrules/swig/CMakeFiles/_cpp.dir/depend: fiberrules/swig/fiberrulesPYTHON_wrap.cxx
fiberrules/swig/CMakeFiles/_cpp.dir/depend: fiberrules/swig/cpp.py
	cd /home/rocio/Desktop/fiberrules/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rocio/Desktop/fiberrules /home/rocio/Desktop/fiberrules/fiberrules/swig /home/rocio/Desktop/fiberrules/build /home/rocio/Desktop/fiberrules/build/fiberrules/swig /home/rocio/Desktop/fiberrules/build/fiberrules/swig/CMakeFiles/_cpp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : fiberrules/swig/CMakeFiles/_cpp.dir/depend

