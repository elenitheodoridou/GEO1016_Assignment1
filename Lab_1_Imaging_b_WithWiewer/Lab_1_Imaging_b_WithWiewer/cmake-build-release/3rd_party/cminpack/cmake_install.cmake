# Install script for directory: C:/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_b_WithWiewer/Lab_1_Imaging_b_WithWiewer/3rd_party/cminpack

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/Lab_1_Imaging_b_WithWiewer")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibraryx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_b_WithWiewer/Lab_1_Imaging_b_WithWiewer/cmake-build-release/lib/3rd_cminpack.lib")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcminpack_hdrsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/cminpack-" TYPE FILE FILES
    "C:/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_b_WithWiewer/Lab_1_Imaging_b_WithWiewer/3rd_party/cminpack/cminpack.h"
    "C:/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_b_WithWiewer/Lab_1_Imaging_b_WithWiewer/3rd_party/cminpack/minpack.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("C:/Users/etheo/Documents/Geomatics/Q4/GEO1016/Assignment1/Lab_1_Imaging_b_WithWiewer/Lab_1_Imaging_b_WithWiewer/cmake-build-release/3rd_party/cminpack/cmake/cmake_install.cmake")

endif()

