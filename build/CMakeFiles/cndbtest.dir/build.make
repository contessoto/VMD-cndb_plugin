# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/d/Dev/VMD-cndb_plugin

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/d/Dev/VMD-cndb_plugin/build

# Include any dependencies generated for this target.
include CMakeFiles/cndbtest.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/cndbtest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cndbtest.dir/flags.make

CMakeFiles/cndbtest.dir/cndbtest.c.o: CMakeFiles/cndbtest.dir/flags.make
CMakeFiles/cndbtest.dir/cndbtest.c.o: ../cndbtest.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/Dev/VMD-cndb_plugin/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/cndbtest.dir/cndbtest.c.o"
	/usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/cndbtest.dir/cndbtest.c.o   -c /mnt/d/Dev/VMD-cndb_plugin/cndbtest.c

CMakeFiles/cndbtest.dir/cndbtest.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/cndbtest.dir/cndbtest.c.i"
	/usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /mnt/d/Dev/VMD-cndb_plugin/cndbtest.c > CMakeFiles/cndbtest.dir/cndbtest.c.i

CMakeFiles/cndbtest.dir/cndbtest.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/cndbtest.dir/cndbtest.c.s"
	/usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /mnt/d/Dev/VMD-cndb_plugin/cndbtest.c -o CMakeFiles/cndbtest.dir/cndbtest.c.s

CMakeFiles/cndbtest.dir/cndbplugin.c.o: CMakeFiles/cndbtest.dir/flags.make
CMakeFiles/cndbtest.dir/cndbplugin.c.o: ../cndbplugin.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/Dev/VMD-cndb_plugin/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/cndbtest.dir/cndbplugin.c.o"
	/usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/cndbtest.dir/cndbplugin.c.o   -c /mnt/d/Dev/VMD-cndb_plugin/cndbplugin.c

CMakeFiles/cndbtest.dir/cndbplugin.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/cndbtest.dir/cndbplugin.c.i"
	/usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /mnt/d/Dev/VMD-cndb_plugin/cndbplugin.c > CMakeFiles/cndbtest.dir/cndbplugin.c.i

CMakeFiles/cndbtest.dir/cndbplugin.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/cndbtest.dir/cndbplugin.c.s"
	/usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /mnt/d/Dev/VMD-cndb_plugin/cndbplugin.c -o CMakeFiles/cndbtest.dir/cndbplugin.c.s

# Object files for target cndbtest
cndbtest_OBJECTS = \
"CMakeFiles/cndbtest.dir/cndbtest.c.o" \
"CMakeFiles/cndbtest.dir/cndbplugin.c.o"

# External object files for target cndbtest
cndbtest_EXTERNAL_OBJECTS =

cndbtest: CMakeFiles/cndbtest.dir/cndbtest.c.o
cndbtest: CMakeFiles/cndbtest.dir/cndbplugin.c.o
cndbtest: CMakeFiles/cndbtest.dir/build.make
cndbtest: libcndb.a
cndbtest: /home/vinicius/anaconda3/lib/libhdf5.so
cndbtest: /home/vinicius/anaconda3/x86_64-conda_cos6-linux-gnu/sysroot/usr/lib/librt.so
cndbtest: /home/vinicius/anaconda3/x86_64-conda_cos6-linux-gnu/sysroot/usr/lib/libpthread.so
cndbtest: /home/vinicius/anaconda3/lib/libz.so
cndbtest: /home/vinicius/anaconda3/x86_64-conda_cos6-linux-gnu/sysroot/usr/lib/libdl.so
cndbtest: /home/vinicius/anaconda3/x86_64-conda_cos6-linux-gnu/sysroot/usr/lib/libm.so
cndbtest: /home/vinicius/anaconda3/lib/libhdf5_hl.so
cndbtest: CMakeFiles/cndbtest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/d/Dev/VMD-cndb_plugin/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable cndbtest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cndbtest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cndbtest.dir/build: cndbtest

.PHONY : CMakeFiles/cndbtest.dir/build

CMakeFiles/cndbtest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cndbtest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cndbtest.dir/clean

CMakeFiles/cndbtest.dir/depend:
	cd /mnt/d/Dev/VMD-cndb_plugin/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/d/Dev/VMD-cndb_plugin /mnt/d/Dev/VMD-cndb_plugin /mnt/d/Dev/VMD-cndb_plugin/build /mnt/d/Dev/VMD-cndb_plugin/build /mnt/d/Dev/VMD-cndb_plugin/build/CMakeFiles/cndbtest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cndbtest.dir/depend
