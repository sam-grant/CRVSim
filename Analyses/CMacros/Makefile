C++ = g++

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

Target = CompCRV.out

all:$(Target)

CompCRV.out: CompCRV.C
	${C++} -o $@ CompCRV.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm *.out