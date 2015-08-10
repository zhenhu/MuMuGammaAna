import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('MuMuGammaRootupler',
                          dimuons = cms.InputTag("onia2MuMuPAT"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          conversions = cms.InputTag("PhotonCandidates","conversions"), # "conversions"
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          onia_pdgid = cms.uint32(443),
                          onia_mass_cuts = cms.vdouble(0.1,500.0),
                          isMC = cms.bool(False),
                          OnlyBest = cms.bool(False),
                          OnlyGen = cms.bool(False),
								  upsilon_mass = cms.double(9.4603), #9.4603. 3.0969
								  triggerCuts = cms.uint32(36)  # 36 for Upsilon, 73 for Jpsi
                          )
