# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/egor/mipt-sle/tridiagonal_matrix

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/egor/mipt-sle/tridiagonal_matrix/build

# Include any dependencies generated for this target.
include CMakeFiles/tester.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/tester.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/tester.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tester.dir/flags.make

CMakeFiles/tester.dir/tester.cpp.o: CMakeFiles/tester.dir/flags.make
CMakeFiles/tester.dir/tester.cpp.o: ../tester.cpp
CMakeFiles/tester.dir/tester.cpp.o: CMakeFiles/tester.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/egor/mipt-sle/tridiagonal_matrix/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/tester.dir/tester.cpp.o"
	/usr/local/gcc-12.2.0/bin/g++-12.2 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tester.dir/tester.cpp.o -MF CMakeFiles/tester.dir/tester.cpp.o.d -o CMakeFiles/tester.dir/tester.cpp.o -c /home/egor/mipt-sle/tridiagonal_matrix/tester.cpp

CMakeFiles/tester.dir/tester.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tester.dir/tester.cpp.i"
	/usr/local/gcc-12.2.0/bin/g++-12.2 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/egor/mipt-sle/tridiagonal_matrix/tester.cpp > CMakeFiles/tester.dir/tester.cpp.i

CMakeFiles/tester.dir/tester.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tester.dir/tester.cpp.s"
	/usr/local/gcc-12.2.0/bin/g++-12.2 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/egor/mipt-sle/tridiagonal_matrix/tester.cpp -o CMakeFiles/tester.dir/tester.cpp.s

# Object files for target tester
tester_OBJECTS = \
"CMakeFiles/tester.dir/tester.cpp.o"

# External object files for target tester
tester_EXTERNAL_OBJECTS =

tester: CMakeFiles/tester.dir/tester.cpp.o
tester: CMakeFiles/tester.dir/build.make
tester: lib/libgtest_main.a
tester: ../tridiag/build/libtridiag.a
tester: lib/libgtest.a
tester: CMakeFiles/tester.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/egor/mipt-sle/tridiagonal_matrix/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable tester"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tester.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/tester.dir/build: tester
.PHONY : CMakeFiles/tester.dir/build

CMakeFiles/tester.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tester.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tester.dir/clean

CMakeFiles/tester.dir/depend:
	cd /home/egor/mipt-sle/tridiagonal_matrix/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/egor/mipt-sle/tridiagonal_matrix /home/egor/mipt-sle/tridiagonal_matrix /home/egor/mipt-sle/tridiagonal_matrix/build /home/egor/mipt-sle/tridiagonal_matrix/build /home/egor/mipt-sle/tridiagonal_matrix/build/CMakeFiles/tester.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/tester.dir/depend

