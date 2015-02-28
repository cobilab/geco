# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

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
CMAKE_SOURCE_DIR = /home/x/git/geco

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/x/git/geco

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running interactive CMake command-line interface..."
	/usr/bin/cmake -i .
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/x/git/geco/CMakeFiles /home/x/git/geco/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/x/git/geco/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named GeCo

# Build rule for target.
GeCo: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GeCo
.PHONY : GeCo

# fast build rule for target.
GeCo/fast:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/build
.PHONY : GeCo/fast

#=============================================================================
# Target rules for targets named GeDe

# Build rule for target.
GeDe: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GeDe
.PHONY : GeDe

# fast build rule for target.
GeDe/fast:
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/build
.PHONY : GeDe/fast

arith.o: arith.c.o
.PHONY : arith.o

# target to build an object file
arith.c.o:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/arith.c.o
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/arith.c.o
.PHONY : arith.c.o

arith.i: arith.c.i
.PHONY : arith.i

# target to preprocess a source file
arith.c.i:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/arith.c.i
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/arith.c.i
.PHONY : arith.c.i

arith.s: arith.c.s
.PHONY : arith.s

# target to generate assembly for a file
arith.c.s:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/arith.c.s
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/arith.c.s
.PHONY : arith.c.s

arith_aux.o: arith_aux.c.o
.PHONY : arith_aux.o

# target to build an object file
arith_aux.c.o:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/arith_aux.c.o
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/arith_aux.c.o
.PHONY : arith_aux.c.o

arith_aux.i: arith_aux.c.i
.PHONY : arith_aux.i

# target to preprocess a source file
arith_aux.c.i:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/arith_aux.c.i
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/arith_aux.c.i
.PHONY : arith_aux.c.i

arith_aux.s: arith_aux.c.s
.PHONY : arith_aux.s

# target to generate assembly for a file
arith_aux.c.s:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/arith_aux.c.s
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/arith_aux.c.s
.PHONY : arith_aux.c.s

bitio.o: bitio.c.o
.PHONY : bitio.o

# target to build an object file
bitio.c.o:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/bitio.c.o
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/bitio.c.o
.PHONY : bitio.c.o

bitio.i: bitio.c.i
.PHONY : bitio.i

# target to preprocess a source file
bitio.c.i:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/bitio.c.i
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/bitio.c.i
.PHONY : bitio.c.i

bitio.s: bitio.c.s
.PHONY : bitio.s

# target to generate assembly for a file
bitio.c.s:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/bitio.c.s
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/bitio.c.s
.PHONY : bitio.c.s

buffer.o: buffer.c.o
.PHONY : buffer.o

# target to build an object file
buffer.c.o:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/buffer.c.o
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/buffer.c.o
.PHONY : buffer.c.o

buffer.i: buffer.c.i
.PHONY : buffer.i

# target to preprocess a source file
buffer.c.i:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/buffer.c.i
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/buffer.c.i
.PHONY : buffer.c.i

buffer.s: buffer.c.s
.PHONY : buffer.s

# target to generate assembly for a file
buffer.c.s:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/buffer.c.s
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/buffer.c.s
.PHONY : buffer.c.s

common.o: common.c.o
.PHONY : common.o

# target to build an object file
common.c.o:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/common.c.o
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/common.c.o
.PHONY : common.c.o

common.i: common.c.i
.PHONY : common.i

# target to preprocess a source file
common.c.i:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/common.c.i
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/common.c.i
.PHONY : common.c.i

common.s: common.c.s
.PHONY : common.s

# target to generate assembly for a file
common.c.s:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/common.c.s
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/common.c.s
.PHONY : common.c.s

context.o: context.c.o
.PHONY : context.o

# target to build an object file
context.c.o:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/context.c.o
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/context.c.o
.PHONY : context.c.o

context.i: context.c.i
.PHONY : context.i

# target to preprocess a source file
context.c.i:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/context.c.i
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/context.c.i
.PHONY : context.c.i

context.s: context.c.s
.PHONY : context.s

# target to generate assembly for a file
context.c.s:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/context.c.s
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/context.c.s
.PHONY : context.c.s

geco.o: geco.c.o
.PHONY : geco.o

# target to build an object file
geco.c.o:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/geco.c.o
.PHONY : geco.c.o

geco.i: geco.c.i
.PHONY : geco.i

# target to preprocess a source file
geco.c.i:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/geco.c.i
.PHONY : geco.c.i

geco.s: geco.c.s
.PHONY : geco.s

# target to generate assembly for a file
geco.c.s:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/geco.c.s
.PHONY : geco.c.s

gede.o: gede.c.o
.PHONY : gede.o

# target to build an object file
gede.c.o:
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/gede.c.o
.PHONY : gede.c.o

gede.i: gede.c.i
.PHONY : gede.i

# target to preprocess a source file
gede.c.i:
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/gede.c.i
.PHONY : gede.c.i

gede.s: gede.c.s
.PHONY : gede.s

# target to generate assembly for a file
gede.c.s:
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/gede.c.s
.PHONY : gede.c.s

mem.o: mem.c.o
.PHONY : mem.o

# target to build an object file
mem.c.o:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/mem.c.o
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/mem.c.o
.PHONY : mem.c.o

mem.i: mem.c.i
.PHONY : mem.i

# target to preprocess a source file
mem.c.i:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/mem.c.i
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/mem.c.i
.PHONY : mem.c.i

mem.s: mem.c.s
.PHONY : mem.s

# target to generate assembly for a file
mem.c.s:
	$(MAKE) -f CMakeFiles/GeCo.dir/build.make CMakeFiles/GeCo.dir/mem.c.s
	$(MAKE) -f CMakeFiles/GeDe.dir/build.make CMakeFiles/GeDe.dir/mem.c.s
.PHONY : mem.c.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... GeCo"
	@echo "... GeDe"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... arith.o"
	@echo "... arith.i"
	@echo "... arith.s"
	@echo "... arith_aux.o"
	@echo "... arith_aux.i"
	@echo "... arith_aux.s"
	@echo "... bitio.o"
	@echo "... bitio.i"
	@echo "... bitio.s"
	@echo "... buffer.o"
	@echo "... buffer.i"
	@echo "... buffer.s"
	@echo "... common.o"
	@echo "... common.i"
	@echo "... common.s"
	@echo "... context.o"
	@echo "... context.i"
	@echo "... context.s"
	@echo "... geco.o"
	@echo "... geco.i"
	@echo "... geco.s"
	@echo "... gede.o"
	@echo "... gede.i"
	@echo "... gede.s"
	@echo "... mem.o"
	@echo "... mem.i"
	@echo "... mem.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

