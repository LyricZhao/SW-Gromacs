# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build

# Utility rule file for cppcheck.

# Include the progress variables for this target.
include tests/CMakeFiles/cppcheck.dir/progress.make

tests/CMakeFiles/cppcheck:
	cd /home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/tests && /usr/bin/cmake -E echo "cppcheck was not found by CMake. Rerun CMake specifying CPPCHECK_EXECUTABLE."

cppcheck: tests/CMakeFiles/cppcheck
cppcheck: tests/CMakeFiles/cppcheck.dir/build.make
.PHONY : cppcheck

# Rule to build all files generated by this target.
tests/CMakeFiles/cppcheck.dir/build: cppcheck
.PHONY : tests/CMakeFiles/cppcheck.dir/build

tests/CMakeFiles/cppcheck.dir/clean:
	cd /home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/cppcheck.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/cppcheck.dir/clean

tests/CMakeFiles/cppcheck.dir/depend:
	cd /home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5 /home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/tests /home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build /home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/tests /home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/tests/CMakeFiles/cppcheck.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/cppcheck.dir/depend
