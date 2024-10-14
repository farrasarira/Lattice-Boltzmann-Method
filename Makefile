# LBM Makefile

TARGET= LBM
CXX= g++-11
CXXFLAGS=-Wall -O3 -ffast-math -fopenmp -march=native
SRCDIR= ./src
OBJDIR= ./obj
RESTART= ./restart
SOURCE= $(wildcard $(SRCDIR)/*.cpp)
OBJECT= $(addprefix $(OBJDIR)/, $(notdir $(SOURCE:.cpp=.o)))
OUTPUT= out.$(TARGET).txt field*.vtr output0* rst stt output0*.vtm geometryflag.vtm geometryflag stop *.txt



$(TARGET): $(OBJECT)
	$(CXX) $(CXXFLAGS) $(shell pkg-config --cflags cantera) -o $(TARGET) $(OBJECT) $(shell pkg-config --libs cantera)
	export OMP_NUM_THREADS=1

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp 
	-mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $< -I./src/headers

.PHONY: clean
clean: 
	rm -rf $(TARGET) $(OBJECT) $(OUTPUT) $(OBJDIR) $(RESTART)
output_clean:
	rm -rf $(OUTPUT)
restart_clean:
	rm -rf *.dat