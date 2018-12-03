import FWCore.ParameterSet.Config as cms

process = cms.Process("mtdhi")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'file:/eos/cms/store/group/phys_heavyions/flowcorr/step1.root',

#'root://xrootd.cmsaf.mit.edu//store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_1.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_2.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_3.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_4.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_5.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_6.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_7.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_8.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_9.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_10.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_11.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_12.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_13.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_14.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_15.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_16.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_17.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_18.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_19.root',
#'/store/user/davidlw/Hydjet_Quenched_MinBias_5020GeV_PhaseII/GEN_SIM_102X_test_v2/180924_215458/0000/step1_20.root',
                ),
#secondaryFileNames = cms.untracked.vstring(
#)
                            )
process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('mtdhi_woETL.root')
                                   )

process.load("MTDHIAnalysis.GenPIDMatching.genpidmatching_cfi")

process.mtdhiana.isETL = cms.untracked.bool(False)

process.mtdhiana_20ps_3sigma = process.mtdhiana.clone()
process.mtdhiana_30ps_3sigma = process.mtdhiana.clone()
process.mtdhiana_50ps_3sigma = process.mtdhiana.clone()
process.mtdhiana_70ps_3sigma = process.mtdhiana.clone()

process.mtdhiana_20ps_2sigma = process.mtdhiana.clone()
process.mtdhiana_20ps_2sigma.nSigmaT = cms.untracked.double(2.0)
process.mtdhiana_30ps_2sigma = process.mtdhiana_20ps_2sigma.clone()
process.mtdhiana_50ps_2sigma = process.mtdhiana_20ps_2sigma.clone()
process.mtdhiana_70ps_2sigma = process.mtdhiana_20ps_2sigma.clone()

process.mtdhiana_20ps_1sigma = process.mtdhiana.clone()
process.mtdhiana_20ps_1sigma.nSigmaT = cms.untracked.double(1.0)
process.mtdhiana_30ps_1sigma = process.mtdhiana_20ps_1sigma.clone()
process.mtdhiana_50ps_1sigma = process.mtdhiana_20ps_1sigma.clone()
process.mtdhiana_70ps_1sigma = process.mtdhiana_20ps_1sigma.clone()

process.mtdhiana_20ps_3sigma.sigmaT = cms.untracked.double(0.02)
process.mtdhiana_30ps_3sigma.sigmaT = cms.untracked.double(0.03)
process.mtdhiana_50ps_3sigma.sigmaT = cms.untracked.double(0.05)
process.mtdhiana_70ps_3sigma.sigmaT = cms.untracked.double(0.07)

process.mtdhiana_20ps_2sigma.sigmaT = cms.untracked.double(0.02)
process.mtdhiana_30ps_2sigma.sigmaT = cms.untracked.double(0.03)
process.mtdhiana_50ps_2sigma.sigmaT = cms.untracked.double(0.05)
process.mtdhiana_70ps_2sigma.sigmaT = cms.untracked.double(0.07)

process.mtdhiana_20ps_1sigma.sigmaT = cms.untracked.double(0.02)
process.mtdhiana_30ps_1sigma.sigmaT = cms.untracked.double(0.03)
process.mtdhiana_50ps_1sigma.sigmaT = cms.untracked.double(0.05)
process.mtdhiana_70ps_1sigma.sigmaT = cms.untracked.double(0.07)

#process.p = cms.Path(process.mtdhiana_20ps_1sigma)
#process.p1 = cms.Path(process.mtdhiana_30ps_1sigma)
#process.p2 = cms.Path(process.mtdhiana_50ps_1sigma)
#process.p3 = cms.Path(process.mtdhiana_70ps_1sigma)

#process.p4 = cms.Path(process.mtdhiana_20ps_2sigma)
#process.p5 = cms.Path(process.mtdhiana_30ps_2sigma)
#process.p6 = cms.Path(process.mtdhiana_50ps_2sigma)
#process.p7 = cms.Path(process.mtdhiana_70ps_2sigma)

process.p8 = cms.Path(process.mtdhiana_20ps_3sigma)
process.p9 = cms.Path(process.mtdhiana_30ps_3sigma)
process.p10 = cms.Path(process.mtdhiana_50ps_3sigma)
process.p11 = cms.Path(process.mtdhiana_70ps_3sigma)
