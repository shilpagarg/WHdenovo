# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades

# Include any dependencies generated for this target.
include common/modules/coverage_model/CMakeFiles/coverage_model.dir/depend.make

# Include the progress variables for this target.
include common/modules/coverage_model/CMakeFiles/coverage_model.dir/progress.make

# Include the compile flags for this target's objects.
include common/modules/coverage_model/CMakeFiles/coverage_model.dir/flags.make

common/modules/coverage_model/CMakeFiles/coverage_model.dir/kmer_coverage_model.cpp.o: common/modules/coverage_model/CMakeFiles/coverage_model.dir/flags.make
common/modules/coverage_model/CMakeFiles/coverage_model.dir/kmer_coverage_model.cpp.o: /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/common/modules/coverage_model/kmer_coverage_model.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object common/modules/coverage_model/CMakeFiles/coverage_model.dir/kmer_coverage_model.cpp.o"
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/common/modules/coverage_model && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/coverage_model.dir/kmer_coverage_model.cpp.o -c /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/common/modules/coverage_model/kmer_coverage_model.cpp

common/modules/coverage_model/CMakeFiles/coverage_model.dir/kmer_coverage_model.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/coverage_model.dir/kmer_coverage_model.cpp.i"
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/common/modules/coverage_model && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/common/modules/coverage_model/kmer_coverage_model.cpp > CMakeFiles/coverage_model.dir/kmer_coverage_model.cpp.i

common/modules/coverage_model/CMakeFiles/coverage_model.dir/kmer_coverage_model.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/coverage_model.dir/kmer_coverage_model.cpp.s"
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/common/modules/coverage_model && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/common/modules/coverage_model/kmer_coverage_model.cpp -o CMakeFiles/coverage_model.dir/kmer_coverage_model.cpp.s

# Object files for target coverage_model
coverage_model_OBJECTS = \
"CMakeFiles/coverage_model.dir/kmer_coverage_model.cpp.o"

# External object files for target coverage_model
coverage_model_EXTERNAL_OBJECTS =

common/modules/coverage_model/libcoverage_model.a: common/modules/coverage_model/CMakeFiles/coverage_model.dir/kmer_coverage_model.cpp.o
common/modules/coverage_model/libcoverage_model.a: common/modules/coverage_model/CMakeFiles/coverage_model.dir/build.make
common/modules/coverage_model/libcoverage_model.a: common/modules/coverage_model/CMakeFiles/coverage_model.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libcoverage_model.a"
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/common/modules/coverage_model && $(CMAKE_COMMAND) -P CMakeFiles/coverage_model.dir/cmake_clean_target.cmake
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/common/modules/coverage_model && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/coverage_model.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
common/modules/coverage_model/CMakeFiles/coverage_model.dir/build: common/modules/coverage_model/libcoverage_model.a

.PHONY : common/modules/coverage_model/CMakeFiles/coverage_model.dir/build

common/modules/coverage_model/CMakeFiles/coverage_model.dir/clean:
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/common/modules/coverage_model && $(CMAKE_COMMAND) -P CMakeFiles/coverage_model.dir/cmake_clean.cmake
.PHONY : common/modules/coverage_model/CMakeFiles/coverage_model.dir/clean

common/modules/coverage_model/CMakeFiles/coverage_model.dir/depend:
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src/common/modules/coverage_model /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/common/modules/coverage_model /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/common/modules/coverage_model/CMakeFiles/coverage_model.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : common/modules/coverage_model/CMakeFiles/coverage_model.dir/depend
