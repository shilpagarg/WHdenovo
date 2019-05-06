# Install script for directory: /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades/pyyaml3" TYPE FILE FILES
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/__init__.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/composer.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/constructor.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/cyaml.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/dumper.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/emitter.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/error.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/events.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/loader.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/nodes.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/parser.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/reader.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/representer.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/resolver.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/scanner.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/serializer.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/python_libs/pyyaml3/tokens.py"
    )
endif()

