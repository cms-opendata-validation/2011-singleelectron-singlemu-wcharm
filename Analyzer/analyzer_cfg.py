import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.Utilities.FileUtils as FileUtils
import sys

# general stuff
process = cms.Process("Demo")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.suppressWarning = cms.untracked.vstring('DemoAnalyzer')
process.MessageLogger.categories.append('demo')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# set max events number
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000))

# ******* is it MC? *********
gen = 0
mc  = 1
if gen == 1 and mc == 0: 
  sys.exit("Error: gen = 1 requires mc = 1")
# ************************************************

# input and output
#inFileTest = 'file:root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/W1Jet_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/30000/0472896B-D713-E411-8C82-0025905A60D2.root'
inFileTest = 'file:/home/cms-opendata/cmsOpenDataFiles/00FB9565-01B4-E311-9233-00248C0BE014.root'
outFileTest = 'wcharmSelTmp.root'
#cacheSize = cms.untracked.uint32(100*1024*1024),
if len(sys.argv) < 4: print("Usage: cmsRun analyzer_cfg.py <input list> <output file>");inputList = inFileTest;outFile = outFileTest#; sys.exit("Wrong usage!")
else:                 inputList = FileUtils.loadListFromFile(sys.argv[2]);outFile = sys.argv[3]
if len(sys.argv) < 4: process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inputList))
else: process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(*inputList),
)
# ************************************************

# global tag from http://opendata.cern.ch/getting-started/CMS?ln=en
# Before should be done:
#    ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA FT_53_LV5_AN1
#    ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db FT_53_LV5_AN1_RUNA.db
if mc == 0:
  # global tag data
  process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')
  process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
else:
  # global tag MC
  process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db')
  process.GlobalTag.globaltag = 'START53_LV6A1::All'
# ************************************************

# apply JSON
if mc == 0:
  goodJSON = 'data/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt'
  myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')
  process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
  process.source.lumisToProcess.extend(myLumis)
# ************************************************

# Load jet correction services for all jet algoritms
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")

# run analyzer
process.demo = cms.EDAnalyzer('DemoAnalyzer',outFile = cms.string(outFile), gen = cms.int32(gen))
process.p = cms.Path(process.demo)
# ************************************************

# ************************************************
# ************************************************
# ************************************************
# print MC generated particles (only hard scattering)
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
#    src = cms.InputTag("genParticles"),
#    printP4 = cms.untracked.bool(True),
#    printPtEtaPhi = cms.untracked.bool(False),
#    printVertex = cms.untracked.bool(False)
#  )
#process.p = cms.Path( process.printDecay)
#  
# print MC generated particles
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
#                                   src = cms.InputTag("genParticles"),                                                                 
#                                   printP4 = cms.untracked.bool(True),
#                                   printPtEtaPhi = cms.untracked.bool(False),
#                                   printVertex = cms.untracked.bool(False),
#                                   printStatus = cms.untracked.bool(True),
#                                   printIndex = cms.untracked.bool(False),
#                                   status = cms.untracked.vint32( 3 )
#                                   )
#process.p = cms.Path(process.printDecay * process.printTree)
# ************************************************
