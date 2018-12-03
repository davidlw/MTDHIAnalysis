import FWCore.ParameterSet.Config as cms

mtdhiana = cms.EDAnalyzer('GenPIDMatchingJets',
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  GenJetCollection = cms.untracked.InputTag("ak4GenJets"),

  isETL = cms.untracked.bool(True),

  sigmaT = cms.untracked.double(0.03),
  nSigmaT = cms.untracked.double(3.0),

  minEtaJet = cms.untracked.double(0),
  maxEtaJet = cms.untracked.double(4.0),
  minPtJet = cms.untracked.double(0),
  maxPtJet = cms.untracked.double(10000.0),
)
