#!/bin/bash

# compile code (produces two executables)
g++ wcharmMakeHist.cxx -o wcharmMakeHist `root-config --cflags --libs` -lMathMore -std=c++11
g++ wcharmMakePlots.cxx -o wcharmMakePlots `root-config --cflags --libs` -std=c++11

# create needed directories if do not exist yet
mkdir -p data mc hist plots
