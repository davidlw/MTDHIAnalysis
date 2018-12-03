import FWCore.ParameterSet.Config as cms

mtdhiana = cms.EDAnalyzer('GenPIDMatching',
  GenParticleCollection = cms.untracked.InputTag("genParticles"),

  isETL = cms.untracked.bool(True),
  isLambdaC = cms.untracked.bool(False),

  sigmaT = cms.untracked.double(0.03),
  nSigmaT = cms.untracked.double(3.0)
)
