# PYTHON configuration file for class: OffsetTreeMaker
# Author: H. Bandyopadhyay
# Date:  12 - June - 2022

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = readFiles
                  , inputCommands=cms.untracked.vstring('keep *',
                  'drop floatBXVector_gtStage2Digis_CICADAScore_RECO')
                 )

readFiles.extend( [
  '/store/data/Run2024J/PPRefZeroBiasPlusForward0/MINIAOD/PromptReco-v1/000/387/716/00000/00aec353-392a-4dc2-ab44-f5f2e6b2e4d7.root'

] );

isMC = cms.bool(False)

if isMC:
  OutputName = "MC_test"
#  eraName = "Summer20UL18_V2_MC"
  jetType_name = "AK4PFchs" # or "AK4PF"

  process.load( "Configuration.Geometry.GeometryIdeal_cff" )
  process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
  process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )
  from Configuration.AlCa.GlobalTag import GlobalTag
  process.GlobalTag = GlobalTag( process.GlobalTag, '140X_mcRun3_2024_realistic_v26' )

else:
  run = "2024J_test"
  OutputName = "_Run3_Data_"+run+"_v1"

#  eraName = "Winter22Run3"+"_RunD"+"_V2_DATA"
  jetType_name = "AK4PFchs" # or "AK4PF"

  process.load( "Configuration.Geometry.GeometryIdeal_cff" )
  process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
  process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )

  from Configuration.AlCa.GlobalTag import GlobalTag
  process.GlobalTag = GlobalTag( process.GlobalTag, '140X_dataRun3_Prompt_v4' )

  # ZeroBias Trigger
  process.HLTZeroBias =cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    #HLTPaths = cms.vstring('HLT_ZeroBias_part*','HLT_ZeroBias_v*'),
#    HLTPaths = cms.vstring('HLT_ZeroBias_v*'),
    HLTPaths = cms.vstring('HLT_PPRefZeroBias_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw = cms.bool(False)
  )

  #Beam Halo
  process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')

  #HCAL HBHE
  process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
  process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
  process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResultRun2Tight'),
    reverseDecision = cms.bool(False)
  )

process.pf = cms.EDAnalyzer("OffsetTreeMaker_mAOD",
    numSkip = cms.int32(1),                                        # change num, pufile and lumibxfile (parsepileup also), vetomap
    RootFileName = cms.string("Offset" + OutputName + ".root"),
    puFileName = cms.string("pileup_24.txt"),
    jetVetoMapFileName = cms.string("Winter24Prompt24_2024BCDEFGHI.root"),       
#    pubxFileName = cms.string("lumi-per-bx_Run2024I-part1.root"),
    mapName2 = cms.string("jetvetomap_all"),   
    isMC = isMC,
    writeCands = cms.bool(False),
   # trackTag = cms.InputTag("generalTracks"),
   # pfTag = cms.InputTag("particleFlow"),
   # pvTag = cms.InputTag("offlinePrimaryVertices"),
   # muTag = cms.InputTag("addPileupInfo"),
    rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoCTag = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoC0Tag = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
    rhoCCTag = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
   # pfJetTag = cms.InputTag("ak4PFJetsCHS"),
    #era = cms.string(eraName),
#    jet_type = cms.string(jetType_name),
    doL1L2L3Res = cms.bool(False),
    dojetVetoMap = cms.bool(True),
    writeEnergyDeposition = cms.bool(True),
    #miniAOD
    trackTag = cms.InputTag("lostTracks"),
    pfTag = cms.InputTag("packedPFCandidates"),
    pvTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muTag = cms.InputTag("slimmedAddPileupInfo"),
    pfJetTag = cms.InputTag("slimmedJets")
)

process.myseq = cms.Sequence( process.pf )

if isMC :
  process.p = cms.Path( process.myseq )
else:
  process.p = cms.Path( process.HLTZeroBias * 
#                        process.CSCTightHaloFilter *
#                        process.HBHENoiseFilterResultProducer *
#                        process.ApplyBaselineHBHENoiseFilter *
                        process.myseq )
