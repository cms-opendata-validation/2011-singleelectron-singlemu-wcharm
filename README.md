# CMS measurements of associated W + c production at 7 TeV

Relevant CMS publications:
 * Measurement of associated W + c production in pp collisions at sqrt(s) = 7 TeV: JHEP 02 (2014) 013 [arXiv:1310.1138, SMP-12-002]
 * total cross section: JHEP 1608 (2016) 029 [arXiv:1603.02303,TOP-13-004]

For the general description of the analysis see also attached description-wcharm.pdf

There are two parts in this analysis:
 * Analyzer: ntuple production, requires CMSSW (the instructions assume that you will work on a VM properly contextualized for CMS, available from http://opendata.cern.ch/VM/CMS) and network connection; takes ~ X weeks to process the full data + MC samples and ~ XGB free space for the produced ntuples
 * PostAnalyzer: ntuple processing, produces final numbers and plots, standalone code (requires only gcc and ROOT); takes about X minutes

## Instructions how to run the analysis

First you need to create the working area (this step is only needed the first time you setup this program). You can create the working area for this analysis on the VM which has other instances of CMSSW, just keep them in different directories.
```
mkdir WorkingArea
cd ./WorkingArea
cmsrel CMSSW_5_3_32
cd ./CMSSW_5_3_32/src
cmsenv
git clone https://github.com/zenaiev/2011-singleelectron-singlemu-wcharm.git
scram b
cd 2011-singleelectron-singlemu-wcharm/Analyzer
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA FT_53_LV5_AN1
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1 START53_LV6A1
```
(no need to download data/MC input file lists and JSON: provided with the code)

Also the code in PostAnalyzer should be compiled:
```
cd PostAnalyzer
./compile.sh
```

## Running the analysis
Generally, the analysis steps are:
 * run Analyzer/run.sh (look inside first), this processes AOD files (CMS data stored at CERN server, several TB) and produces plain ROOT ntuple files (~XGB), takes ~ X weeks, extensive network access
 * move produced ntuples to PostAnalyzer directories (this step is manual on purpose, in order not to overwrite accidentally ntuples produced taking long time etc.)
 * run PostAnalyzer/ttbarMakeHist to process ROOT ntuples to create histograms (~X mins)
 * run PostAnalyzer/ttbarMakePlots to produce final plots from created histograms (X seconds)

Further description of these steps you can find Analyzer/README.txt and Postanalyzer/README.txt

### Authors: ###### Oleksandr Kot, Oleksandr Zenaiev
