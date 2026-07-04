# LBM Makefile
#
# Select the compiler toolchain with:   make TOOLCHAIN=<name>
#   gnu    (default)  GNU g++            - solid everywhere, on both AMD and Intel
#   clang             LLVM clang++       - comparable to GCC; AMD's AOCC is based on it
#   aocc              AMD AOCC clang++   - AMD's tuned LLVM (source AOCC env first)
#   intel             Intel oneAPI icpx  - sometimes fastest on Intel CPUs (source setvars.sh first)
# The compiler binary can be overridden explicitly, e.g.:  make TOOLCHAIN=gnu GXX=g++-12
#
# NOTE: what matters most on both vendors is building on the target machine so
# -march=native / -xHost picks up the right SIMD ISA (AVX2/AVX-512), and using
# the same toolchain for Cantera if possible - measure before assuming a vendor
# compiler wins.

TARGET= LBM
TOOLCHAIN ?= gnu

ifeq ($(TOOLCHAIN),gnu)
  GXX ?= g++-11
  CXX = $(GXX)
  ARCHFLAG = -march=native
  LTOFLAG  = -flto=auto
else ifeq ($(TOOLCHAIN),clang)
  GXX ?= clang++
  CXX = $(GXX)
  ARCHFLAG = -march=native
  LTOFLAG  = -flto=thin
else ifeq ($(TOOLCHAIN),aocc)
  GXX ?= clang++          # AOCC installs its own clang++; put it first in PATH
  CXX = $(GXX)
  ARCHFLAG = -march=native
  LTOFLAG  = -flto=thin
else ifeq ($(TOOLCHAIN),intel)
  GXX ?= icpx
  CXX = $(GXX)
  ARCHFLAG = -xHost
  LTOFLAG  = -flto
else
  $(error Unknown TOOLCHAIN '$(TOOLCHAIN)': use gnu, clang, aocc or intel)
endif

# -MMD -MP: header dependency tracking, so a plain 'make' after editing
#           defines.hpp rebuilds exactly the affected files (no 'make clean' needed)
CXXFLAGS=-Wall -O3 -ffast-math -fopenmp $(ARCHFLAG) $(LTOFLAG) -MMD -MP
SRCDIR= ./src
OBJDIR= ./obj
RESTART= ./restart
SOURCE= $(wildcard $(SRCDIR)/*.cpp)
OBJECT= $(addprefix $(OBJDIR)/, $(notdir $(SOURCE:.cpp=.o)))
DEPEND= $(OBJECT:.o=.d)
OUTPUT= out.$(TARGET).txt field*.vtr output0* rst stt output0*.vtm geometryflag.vtm geometryflag stop *.txt



$(TARGET): $(OBJECT)
	$(CXX) $(CXXFLAGS) $(shell pkg-config --cflags cantera) -o $(TARGET) $(OBJECT) $(shell pkg-config --libs cantera)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $< -I./src/headers

-include $(DEPEND)

.PHONY: clean
clean:
	rm -rf $(TARGET) $(OBJECT) $(DEPEND) $(OBJDIR) $(RESTART)
output_clean:
	rm -rf $(OUTPUT)
restart_clean:
	rm -rf *.dat
