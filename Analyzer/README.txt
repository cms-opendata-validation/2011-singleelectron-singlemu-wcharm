This part processes CMS input data in AOD (Analysis Object Data) format 
to produce plain ROOT ntuples with needed event contents, which are 
further processed by PostAnalyzer.

General description of contents (find further description inside the files):
   run.sh: this is the main script which you will run to process each 
           data and MC sample
   analyzer_cfg.py: standard CMSSW configuration file for cmsRun 
           (you can run the command 'cmsRun analyzer_cfg.py' to process 
           one input data file)
   src/Analyzer.cc: C++ analysis code (does basic event selection, 
           fill output ROOT ntuples)
   data/ mc/: directories with input file lists
   BuildFile.xml: standard CMSSW file for code compialtion 
           (if you have problems with linking, most likely you need 
           to change something in this file)
   python/analyzer_cfi.py: standard (automatically created) CMSSW file

To run the analysis, look at run.sh
