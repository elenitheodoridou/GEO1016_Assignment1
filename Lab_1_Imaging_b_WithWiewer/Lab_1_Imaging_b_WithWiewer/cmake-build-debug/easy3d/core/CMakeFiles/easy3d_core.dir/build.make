# CMAKE generated file: DO NOT EDIT!
# Generated by "NMake Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

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

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE
NULL=nul
!ENDIF
SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2021.2.3\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2021.2.3\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

# Include any dependencies generated for this target.
include easy3d\core\CMakeFiles\easy3d_core.dir\depend.make
# Include the progress variables for this target.
include easy3d\core\CMakeFiles\easy3d_core.dir\progress.make

# Include the compile flags for this target's objects.
include easy3d\core\CMakeFiles\easy3d_core.dir\flags.make

easy3d\core\CMakeFiles\easy3d_core.dir\graph.cpp.obj: easy3d\core\CMakeFiles\easy3d_core.dir\flags.make
easy3d\core\CMakeFiles\easy3d_core.dir\graph.cpp.obj: ..\easy3d\core\graph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/graph.cpp.obj"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\easy3d_core.dir\graph.cpp.obj /FdCMakeFiles\easy3d_core.dir\easy3d_core.pdb /FS -c C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\graph.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\graph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/graph.cpp.i"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe > CMakeFiles\easy3d_core.dir\graph.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\graph.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\graph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/graph.cpp.s"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_core.dir\graph.cpp.s /c C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\graph.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\kdtree.cpp.obj: easy3d\core\CMakeFiles\easy3d_core.dir\flags.make
easy3d\core\CMakeFiles\easy3d_core.dir\kdtree.cpp.obj: ..\easy3d\core\kdtree.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/kdtree.cpp.obj"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\easy3d_core.dir\kdtree.cpp.obj /FdCMakeFiles\easy3d_core.dir\easy3d_core.pdb /FS -c C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\kdtree.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\kdtree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/kdtree.cpp.i"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe > CMakeFiles\easy3d_core.dir\kdtree.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\kdtree.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\kdtree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/kdtree.cpp.s"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_core.dir\kdtree.cpp.s /c C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\kdtree.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\point_cloud.cpp.obj: easy3d\core\CMakeFiles\easy3d_core.dir\flags.make
easy3d\core\CMakeFiles\easy3d_core.dir\point_cloud.cpp.obj: ..\easy3d\core\point_cloud.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/point_cloud.cpp.obj"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\easy3d_core.dir\point_cloud.cpp.obj /FdCMakeFiles\easy3d_core.dir\easy3d_core.pdb /FS -c C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\point_cloud.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\point_cloud.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/point_cloud.cpp.i"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe > CMakeFiles\easy3d_core.dir\point_cloud.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\point_cloud.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\point_cloud.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/point_cloud.cpp.s"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_core.dir\point_cloud.cpp.s /c C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\point_cloud.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\surface_mesh.cpp.obj: easy3d\core\CMakeFiles\easy3d_core.dir\flags.make
easy3d\core\CMakeFiles\easy3d_core.dir\surface_mesh.cpp.obj: ..\easy3d\core\surface_mesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/surface_mesh.cpp.obj"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\easy3d_core.dir\surface_mesh.cpp.obj /FdCMakeFiles\easy3d_core.dir\easy3d_core.pdb /FS -c C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\surface_mesh.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\surface_mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/surface_mesh.cpp.i"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe > CMakeFiles\easy3d_core.dir\surface_mesh.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\surface_mesh.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\surface_mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/surface_mesh.cpp.s"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_core.dir\surface_mesh.cpp.s /c C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\surface_mesh.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\manifold_builder.cpp.obj: easy3d\core\CMakeFiles\easy3d_core.dir\flags.make
easy3d\core\CMakeFiles\easy3d_core.dir\manifold_builder.cpp.obj: ..\easy3d\core\manifold_builder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/manifold_builder.cpp.obj"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\easy3d_core.dir\manifold_builder.cpp.obj /FdCMakeFiles\easy3d_core.dir\easy3d_core.pdb /FS -c C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\manifold_builder.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\manifold_builder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/manifold_builder.cpp.i"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe > CMakeFiles\easy3d_core.dir\manifold_builder.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\manifold_builder.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\manifold_builder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/manifold_builder.cpp.s"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_core.dir\manifold_builder.cpp.s /c C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core\manifold_builder.cpp
<<
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

# Object files for target easy3d_core
easy3d_core_OBJECTS = \
"CMakeFiles\easy3d_core.dir\graph.cpp.obj" \
"CMakeFiles\easy3d_core.dir\kdtree.cpp.obj" \
"CMakeFiles\easy3d_core.dir\point_cloud.cpp.obj" \
"CMakeFiles\easy3d_core.dir\surface_mesh.cpp.obj" \
"CMakeFiles\easy3d_core.dir\manifold_builder.cpp.obj"

# External object files for target easy3d_core
easy3d_core_EXTERNAL_OBJECTS =

lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\graph.cpp.obj
lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\kdtree.cpp.obj
lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\point_cloud.cpp.obj
lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\surface_mesh.cpp.obj
lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\manifold_builder.cpp.obj
lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\build.make
lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX static library ..\..\lib\easy3d_core.lib"
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	$(CMAKE_COMMAND) -P CMakeFiles\easy3d_core.dir\cmake_clean_target.cmake
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\lib.exe /nologo /machine:X86 /out:..\..\lib\easy3d_core.lib @CMakeFiles\easy3d_core.dir\objects1.rsp 
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug

# Rule to build all files generated by this target.
easy3d\core\CMakeFiles\easy3d_core.dir\build: lib\easy3d_core.lib
.PHONY : easy3d\core\CMakeFiles\easy3d_core.dir\build

easy3d\core\CMakeFiles\easy3d_core.dir\clean:
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core
	$(CMAKE_COMMAND) -P CMakeFiles\easy3d_core.dir\cmake_clean.cmake
	cd C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug
.PHONY : easy3d\core\CMakeFiles\easy3d_core.dir\clean

easy3d\core\CMakeFiles\easy3d_core.dir\depend:
	$(CMAKE_COMMAND) -E cmake_depends "NMake Makefiles" C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\easy3d\core C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core C:\Users\etheo\Documents\Geomatics\Q4\GEO1016\Assignment1\Lab_1_Imaging_b_WithWiewer\Lab_1_Imaging_b_WithWiewer\cmake-build-debug\easy3d\core\CMakeFiles\easy3d_core.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : easy3d\core\CMakeFiles\easy3d_core.dir\depend

