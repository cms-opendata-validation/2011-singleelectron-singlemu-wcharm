#!/bin/bash

# compile code (produces two executables)
<<<<<<< HEAD
g++ -O2 -o wcharmMakeHist `root-config --cflags --libs` -lMathMore -std=c++11 wcharmMakeHist.cxx
g++ -O2 -o wcharmMakePlots `root-config --cflags --libs` -std=c++11 wcharmMakePlots.cxx
=======
g++ wcharmMakeHist.cxx -o wcharmMakeHist `root-config --cflags --libs` -lMathMore -std=c++11
g++ wcharmMakePlots.cxx -o wcharmMakePlots `root-config --cflags --libs` -std=c++11
>>>>>>> f8954ab947d2b688c29761e22d65e8e9227856cc

# create needed directories if do not exist yet
mkdir -p data mc hist plots
