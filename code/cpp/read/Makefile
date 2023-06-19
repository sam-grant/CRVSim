CXX = g++
CXXFLAGS = -g -Wall -Wextra -std=c++11
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

SRC = Plotter.C
TARGET = Plotter.out
LINK = LinkDef.h

all: $(TARGET)

$(TARGET): $(SRC) $(LINK)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC) $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm -f $(TARGET)