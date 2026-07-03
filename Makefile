# LBM Makefile

TARGET= LBM
CXX= g++-11
# -flto=auto : link-time optimization, lets the compiler inline the equilibrium
#              functions across source files (measurable speedup of the hot loops)
# -MMD -MP   : generate header dependency files so incremental builds recompile
#              what is affected by a header change ("make clean" no longer required
#              after editing defines.hpp; a plain "make" rebuilds the right files)
CXXFLAGS=-Wall -O3 -ffast-math -fopenmp -march=native -flto=auto -MMD -MP
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
