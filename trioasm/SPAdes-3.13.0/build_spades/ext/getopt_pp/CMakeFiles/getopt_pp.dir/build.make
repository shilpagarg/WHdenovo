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
include ext/getopt_pp/CMakeFiles/getopt_pp.dir/depend.make

# Include the progress variables for this target.
include ext/getopt_pp/CMakeFiles/getopt_pp.dir/progress.make

# Include the compile flags for this target's objects.
include ext/getopt_pp/CMakeFiles/getopt_pp.dir/flags.make

ext/getopt_pp/CMakeFiles/getopt_pp.dir/getopt_pp.cpp.o: ext/getopt_pp/CMakeFiles/getopt_pp.dir/flags.make
ext/getopt_pp/CMakeFiles/getopt_pp.dir/getopt_pp.cpp.o: /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/getopt_pp/getopt_pp.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object ext/getopt_pp/CMakeFiles/getopt_pp.dir/getopt_pp.cpp.o"
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/ext/getopt_pp && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/getopt_pp.dir/getopt_pp.cpp.o -c /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/getopt_pp/getopt_pp.cpp

ext/getopt_pp/CMakeFiles/getopt_pp.dir/getopt_pp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/getopt_pp.dir/getopt_pp.cpp.i"
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/ext/getopt_pp && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/getopt_pp/getopt_pp.cpp > CMakeFiles/getopt_pp.dir/getopt_pp.cpp.i

ext/getopt_pp/CMakeFiles/getopt_pp.dir/getopt_pp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/getopt_pp.dir/getopt_pp.cpp.s"
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/ext/getopt_pp && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/getopt_pp/getopt_pp.cpp -o CMakeFiles/getopt_pp.dir/getopt_pp.cpp.s

# Object files for target getopt_pp
getopt_pp_OBJECTS = \
"CMakeFiles/getopt_pp.dir/getopt_pp.cpp.o"

# External object files for target getopt_pp
getopt_pp_EXTERNAL_OBJECTS =

ext/getopt_pp/libgetopt_pp.a: ext/getopt_pp/CMakeFiles/getopt_pp.dir/getopt_pp.cpp.o
ext/getopt_pp/libgetopt_pp.a: ext/getopt_pp/CMakeFiles/getopt_pp.dir/build.make
ext/getopt_pp/libgetopt_pp.a: ext/getopt_pp/CMakeFiles/getopt_pp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libgetopt_pp.a"
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/ext/getopt_pp && $(CMAKE_COMMAND) -P CMakeFiles/getopt_pp.dir/cmake_clean_target.cmake
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/ext/getopt_pp && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/getopt_pp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ext/getopt_pp/CMakeFiles/getopt_pp.dir/build: ext/getopt_pp/libgetopt_pp.a

.PHONY : ext/getopt_pp/CMakeFiles/getopt_pp.dir/build

ext/getopt_pp/CMakeFiles/getopt_pp.dir/clean:
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/ext/getopt_pp && $(CMAKE_COMMAND) -P CMakeFiles/getopt_pp.dir/cmake_clean.cmake
.PHONY : ext/getopt_pp/CMakeFiles/getopt_pp.dir/clean

ext/getopt_pp/CMakeFiles/getopt_pp.dir/depend:
	cd /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/src /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/ext/src/getopt_pp /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/ext/getopt_pp /home/ec2-user/whdenovo-v0424/trioasm/SPAdes-3.13.0/build_spades/ext/getopt_pp/CMakeFiles/getopt_pp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ext/getopt_pp/CMakeFiles/getopt_pp.dir/depend
