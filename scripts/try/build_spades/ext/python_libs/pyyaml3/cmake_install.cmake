# Install script for directory: /n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/build_spades")
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
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/__init__.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/composer.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/constructor.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/cyaml.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/dumper.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/emitter.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/error.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/events.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/loader.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/nodes.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/parser.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/reader.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/representer.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/resolver.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/scanner.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/serializer.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/ext/src/python_libs/pyyaml3/tokens.py"
    )
endif()

