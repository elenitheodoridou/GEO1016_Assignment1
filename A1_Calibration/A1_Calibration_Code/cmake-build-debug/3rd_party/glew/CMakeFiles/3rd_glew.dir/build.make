﻿# CMAKE generated file: DO NOT EDIT!
# Generated by "NMake Makefiles" Generator, CMake Version 3.21

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
CMAKE_SOURCE_DIR = "C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code\cmake-build-debug"

# Include any dependencies generated for this target.
include 3rd_party\glew\CMakeFiles\3rd_glew.dir\depend.make
# Include any dependencies generated by the compiler for this target.
include 3rd_party\glew\CMakeFiles\3rd_glew.dir\compiler_depend.make

# Include the progress variables for this target.
include 3rd_party\glew\CMakeFiles\3rd_glew.dir\progress.make

# Include the compile flags for this target's objects.
include 3rd_party\glew\CMakeFiles\3rd_glew.dir\flags.make

3rd_party\glew\CMakeFiles\3rd_glew.dir\src\glew.c.obj: 3rd_party\glew\CMakeFiles\3rd_glew.dir\flags.make
3rd_party\glew\CMakeFiles\3rd_glew.dir\src\glew.c.obj: ..\3rd_party\glew\src\glew.c
3rd_party\glew\CMakeFiles\3rd_glew.dir\src\glew.c.obj: 3rd_party\glew\CMakeFiles\3rd_glew.dir\compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object 3rd_party/glew/CMakeFiles/3rd_glew.dir/src/glew.c.obj"
	cd C:\Users\T420\Desktop\2021-2~1\GEOMAT~1\PERIOD~4\GEO101~1\A1_CAL~1\A1_CAL~1\CMAKE-~1\3RD_PA~1\glew
	$(CMAKE_COMMAND) -E cmake_cl_compile_depends --dep-file=CMakeFiles\3rd_glew.dir\src\glew.c.obj.d --working-dir="C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code\cmake-build-debug\3rd_party\glew" --filter-prefix="Note: including file: " -- C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe @<<
 /nologo $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) /showIncludes /FoCMakeFiles\3rd_glew.dir\src\glew.c.obj /FdCMakeFiles\3rd_glew.dir\3rd_glew.pdb /FS -c "C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code\3rd_party\glew\src\glew.c"
<<
	cd C:\Users\T420\Desktop\2021-2~1\GEOMAT~1\PERIOD~4\GEO101~1\A1_CAL~1\A1_CAL~1\CMAKE-~1

3rd_party\glew\CMakeFiles\3rd_glew.dir\src\glew.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/3rd_glew.dir/src/glew.c.i"
	cd C:\Users\T420\Desktop\2021-2~1\GEOMAT~1\PERIOD~4\GEO101~1\A1_CAL~1\A1_CAL~1\CMAKE-~1\3RD_PA~1\glew
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe > CMakeFiles\3rd_glew.dir\src\glew.c.i @<<
 /nologo $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code\3rd_party\glew\src\glew.c"
<<
	cd C:\Users\T420\Desktop\2021-2~1\GEOMAT~1\PERIOD~4\GEO101~1\A1_CAL~1\A1_CAL~1\CMAKE-~1

3rd_party\glew\CMakeFiles\3rd_glew.dir\src\glew.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/3rd_glew.dir/src/glew.c.s"
	cd C:\Users\T420\Desktop\2021-2~1\GEOMAT~1\PERIOD~4\GEO101~1\A1_CAL~1\A1_CAL~1\CMAKE-~1\3RD_PA~1\glew
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\cl.exe @<<
 /nologo $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) /FoNUL /FAs /FaCMakeFiles\3rd_glew.dir\src\glew.c.s /c "C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code\3rd_party\glew\src\glew.c"
<<
	cd C:\Users\T420\Desktop\2021-2~1\GEOMAT~1\PERIOD~4\GEO101~1\A1_CAL~1\A1_CAL~1\CMAKE-~1

# Object files for target 3rd_glew
3rd_glew_OBJECTS = \
"CMakeFiles\3rd_glew.dir\src\glew.c.obj"

# External object files for target 3rd_glew
3rd_glew_EXTERNAL_OBJECTS =

lib\3rd_glew.lib: 3rd_party\glew\CMakeFiles\3rd_glew.dir\src\glew.c.obj
lib\3rd_glew.lib: 3rd_party\glew\CMakeFiles\3rd_glew.dir\build.make
lib\3rd_glew.lib: 3rd_party\glew\CMakeFiles\3rd_glew.dir\objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library ..\..\lib\3rd_glew.lib"
	cd C:\Users\T420\Desktop\2021-2~1\GEOMAT~1\PERIOD~4\GEO101~1\A1_CAL~1\A1_CAL~1\CMAKE-~1\3RD_PA~1\glew
	$(CMAKE_COMMAND) -P CMakeFiles\3rd_glew.dir\cmake_clean_target.cmake
	cd C:\Users\T420\Desktop\2021-2~1\GEOMAT~1\PERIOD~4\GEO101~1\A1_CAL~1\A1_CAL~1\CMAKE-~1
	cd C:\Users\T420\Desktop\2021-2~1\GEOMAT~1\PERIOD~4\GEO101~1\A1_CAL~1\A1_CAL~1\CMAKE-~1\3RD_PA~1\glew
	C:\PROGRA~2\MICROS~4\2019\COMMUN~1\VC\Tools\MSVC\1429~1.301\bin\Hostx86\x86\lib.exe /nologo /machine:X86 /out:..\..\lib\3rd_glew.lib @CMakeFiles\3rd_glew.dir\objects1.rsp 
	cd C:\Users\T420\Desktop\2021-2~1\GEOMAT~1\PERIOD~4\GEO101~1\A1_CAL~1\A1_CAL~1\CMAKE-~1

# Rule to build all files generated by this target.
3rd_party\glew\CMakeFiles\3rd_glew.dir\build: lib\3rd_glew.lib
.PHONY : 3rd_party\glew\CMakeFiles\3rd_glew.dir\build

3rd_party\glew\CMakeFiles\3rd_glew.dir\clean:
	cd C:\Users\T420\Desktop\2021-2~1\GEOMAT~1\PERIOD~4\GEO101~1\A1_CAL~1\A1_CAL~1\CMAKE-~1\3RD_PA~1\glew
	$(CMAKE_COMMAND) -P CMakeFiles\3rd_glew.dir\cmake_clean.cmake
	cd C:\Users\T420\Desktop\2021-2~1\GEOMAT~1\PERIOD~4\GEO101~1\A1_CAL~1\A1_CAL~1\CMAKE-~1
.PHONY : 3rd_party\glew\CMakeFiles\3rd_glew.dir\clean

3rd_party\glew\CMakeFiles\3rd_glew.dir\depend:
	$(CMAKE_COMMAND) -E cmake_depends "NMake Makefiles" "C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code" "C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code\3rd_party\glew" "C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code\cmake-build-debug" "C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code\cmake-build-debug\3rd_party\glew" "C:\Users\T420\Desktop\2021-2022\Geomatics\Periode 4\GEO1016 Photogrammetry and 3D Computer Vision\A1_Calibration\A1_Calibration_Code\cmake-build-debug\3rd_party\glew\CMakeFiles\3rd_glew.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : 3rd_party\glew\CMakeFiles\3rd_glew.dir\depend

