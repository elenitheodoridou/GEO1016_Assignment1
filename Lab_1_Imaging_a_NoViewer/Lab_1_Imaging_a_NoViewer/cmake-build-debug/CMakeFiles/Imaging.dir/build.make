# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Imaging.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/Imaging.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Imaging.dir/flags.make

CMakeFiles/Imaging.dir/main.cpp.o: CMakeFiles/Imaging.dir/flags.make
CMakeFiles/Imaging.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Imaging.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Imaging.dir/main.cpp.o -c /mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/main.cpp

CMakeFiles/Imaging.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Imaging.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/main.cpp > CMakeFiles/Imaging.dir/main.cpp.i

CMakeFiles/Imaging.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Imaging.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/main.cpp -o CMakeFiles/Imaging.dir/main.cpp.s

# Object files for target Imaging
Imaging_OBJECTS = \
"CMakeFiles/Imaging.dir/main.cpp.o"

# External object files for target Imaging
Imaging_EXTERNAL_OBJECTS =

bin/Imaging: CMakeFiles/Imaging.dir/main.cpp.o
bin/Imaging: CMakeFiles/Imaging.dir/build.make
bin/Imaging: lib/libeasy3d_core.a
bin/Imaging: lib/libeasy3d_util.a
bin/Imaging: lib/lib3rd_glog.a
bin/Imaging: CMakeFiles/Imaging.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bin/Imaging"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Imaging.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Imaging.dir/build: bin/Imaging
.PHONY : CMakeFiles/Imaging.dir/build

CMakeFiles/Imaging.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Imaging.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Imaging.dir/clean

CMakeFiles/Imaging.dir/depend:
	cd /mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer /mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer /mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/cmake-build-debug /mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/cmake-build-debug /mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/cmake-build-debug/CMakeFiles/Imaging.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Imaging.dir/depend

