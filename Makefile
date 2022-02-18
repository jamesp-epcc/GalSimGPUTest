#	g++ -o test1 test1.cpp -std=c++11 -I/home/james/lsst/GalSim/include -L/home/james/lsst/GalSim/build/shared_clib -lgalsim

test1: test1.cpp
	CC -o test1 test1.cpp -std=c++11 -mp=gpu -target-accel=nvidia80 -I/global/homes/j/jamesp/gpu/GalSim/include -L/global/homes/j/jamesp/gpu/GalSim/src -lgalsim
