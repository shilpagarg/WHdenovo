# Install script for directory: /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithAsserts")
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
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../spades.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../metaspades.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../plasmidspades.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../rnaspades.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../truspades.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE FILE FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../spades_init.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades" TYPE FILE FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../VERSION")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades" TYPE DIRECTORY FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../test_dataset")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades" TYPE DIRECTORY FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../test_dataset_truspades")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades" TYPE DIRECTORY FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../test_dataset_plasmid")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades" TYPE FILE FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../manual.html")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades" TYPE FILE FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../README.md")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades" TYPE FILE FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../truspades_manual.html")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades" TYPE FILE FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../rnaspades_manual.html")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades" TYPE FILE FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../LICENSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades" TYPE FILE FILES "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/../GPLv2.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/ext/cmake_install.cmake")
  include("/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/ext_tools/cmake_install.cmake")
  include("/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/common/cmake_install.cmake")
  include("/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/projects/cmake_install.cmake")
  include("/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/spades_pipeline/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
