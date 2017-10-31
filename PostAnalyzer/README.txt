This part processes ROOT ntuples (output of Analyzer, see 
Analyzer/README.txt) to produce ROOT histograms (wcharmMakeHist.cxx) 
and final plots and numbers from them (wcharmMakePlots.cxx).
This is pure C++ and ROOT code: does not require CMSSW, works outside VM 
(although it can work on VM also, of course).

General description of contents (find further description inside the files):
   wcharmMakeHist.cxx: master file to produce histograms
   wcharm_eventReco.h: ttbar event reconstruction
   wcharm_selection.h: ttbar event selection
   wcharm_kinReco.h: kinematic reconstruction
   wcharm_tree.h: tree structure of input ROOT ntuples
   wcharm_settings.h: global settings (directory names)
   wcharmMakePlots.cxx: master file to produce final plots and numbers
   wcharm_plots.h: helper file for plotting

To run the analysis, make sure input ntuples are in place, for default 
directory structure you need to run from the root analysis directory:
mv Analyzer/ntuples-data PostAnalyzer/ntuples-data
mv Analyzer/ntuples-mc PostAnalyzer/ntuples-mc
then compile the code:
./compile.sh
and run two commands:
./wcharmMakeHist
./wcharmMakePlots

Also you could do only the last step (plotting) by using "reference" 
histograms produced with the full samples and available with the code 
(PostAnalyzer/hist-REF directory), for this modify wcharm_settings.h. 
Another application of the "reference" histograms could be for 
validation (produce new histograms and compare to the reference ones).
