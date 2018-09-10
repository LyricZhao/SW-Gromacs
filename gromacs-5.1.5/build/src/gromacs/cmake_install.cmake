# Install script for directory: /home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/src/gromacs

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local/gromacs")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "0")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/gmxlib/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/mdlib/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/listed-forces/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/commandline/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/domdec/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/ewald/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/fft/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/linearalgebra/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/math/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/mdrunutility/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/random/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/onlinehelp/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/options/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/pbcutil/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/timing/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/topology/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/utility/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/fileio/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/swap/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/essentialdynamics/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/pulling/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/simd/cmake_install.cmake")
  INCLUDE("/home/export/base/CPC/cpc051/zhaocg/gromacs-5.1.5/build/src/gromacs/imd/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

