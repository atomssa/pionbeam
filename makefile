ROOT_INC        = $(shell root-config --incdir)
ROOT_LIBS       = $(shell root-config --libs)
CXXFLAGS        = -Wall -std=c++0x
CXXFLAGS_GPROF  = -Wall -pg

gen: HBeam.h pion_generator.cpp
	g++ pion_generator.cpp -o gen $(CXXFLAGS) -I$(ROOT_INC) $(ROOT_LIBS) -lMinuit

gen_gprof: HBeam.h pion_generator.cpp
	g++ pion_generator.cpp -o gen_gprof $(CXXFLAGS_GPROF) -I$(ROOT_INC) $(ROOT_LIBS) -lMinuit

clean:
	rm -vf gen gen_prof
