# Install script for directory: /mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/3rd_party/cminpack

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibraryx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/cmake-build-release/lib/lib3rd_cminpack.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcminpack_hdrsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/cminpack-" TYPE FILE FILES
    "/mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/3rd_party/cminpack/cminpack.h"
    "/mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/3rd_party/cminpack/minpack.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/mnt/c/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_a_NoViewer/Lab_1_Imaging_a_NoViewer/cmake-build-release/3rd_party/cminpack/cmake/cmake_install.cmake")

endif()

