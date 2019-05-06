# Install script for directory: /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades/spades_pipeline" TYPE FILE FILES
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/hammer_logic.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/process_cfg.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/spades_logic.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/corrector_logic.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/support.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/options_storage.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/lucigen_nxmate.py"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades/spades_pipeline/truspades" TYPE FILE FILES
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/truspades/reference_construction.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/truspades/moleculo_filter_contigs.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/truspades/break_by_coverage.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/truspades/moleculo_postprocessing.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/truspades/barcode_extraction.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/truspades/generate_quality.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/truspades/id_generation.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/truspades/launch_options.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/truspades/string_dist_utils.py"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades/spades_pipeline/common" TYPE FILE FILES
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/common/alignment.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/common/parallel_launcher.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/common/sam_parser.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/common/SeqIO.py"
    "/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/spades_pipeline/common/__init__.py"
    )
endif()

