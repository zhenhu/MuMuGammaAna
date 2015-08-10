import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v0', '') 
#process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v1'  # or some other global tag depending on your CMSSW release and sample. 

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v6/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_2_1_Ajl.root')
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Rootuple.root'),
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load('MuMuGammaAna.Onia.MuMuGammaRootupler_cfi')
process.p = cms.Path(process.rootuple)

process.rootuple.upsilon_mass = cms.double(9.4603)
process.rootuple.triggerCuts = cms.uint32(36)      # 36 for upsilon, 73 for Jpsi
process.rootuple.isMC = cms.bool(False)                 # is mc?
process.rootuple.onia_mass_cuts = cms.vdouble(0.1,500)    # you may need to adjust this
