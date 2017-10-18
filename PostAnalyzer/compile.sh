#!/bin/bash

# compile code (produces two executables)
g++ -o wcharmMakeHist `root-config --cflags --libs` -lMathMore -std=c++11 wcharmMakeHist.cxx
g++ -g -o wcharmMakePlots `root-config --cflags --libs` -std=c++11 wcharmMakePlots.cxx

# create needed directories if do not exist yet
mkdir -p data mc hist plots
