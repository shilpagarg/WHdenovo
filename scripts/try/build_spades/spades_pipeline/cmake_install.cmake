# Install script for directory: /n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades/spades_pipeline" TYPE FILE FILES
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/hammer_logic.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/process_cfg.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/spades_logic.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/dipspades_logic.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/corrector_logic.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/support.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/options_storage.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/lucigen_nxmate.py"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades/spades_pipeline/truspades" TYPE FILE FILES
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/truspades/reference_construction.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/truspades/moleculo_filter_contigs.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/truspades/break_by_coverage.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/truspades/moleculo_postprocessing.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/truspades/barcode_extraction.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/truspades/generate_quality.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/truspades/id_generation.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/truspades/launch_options.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/truspades/string_dist_utils.py"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades/spades_pipeline/common" TYPE FILE FILES
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/common/alignment.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/common/parallel_launcher.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/common/sam_parser.py"
    "/n/groups/church/shilpa/new/butter_ex/butter_ex/scripts/try/src/spades_pipeline/common/SeqIO.py"
    )
endif()

