import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

process = cms.Process("tnp")

###################################################################
options = dict()
varOptions = VarParsing('analysis')
varOptions.register(
    "isMC",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Compute MC efficiencies"
    )

varOptions.parseArguments()

isolationDef = '(pfIsolationR04().sumChargedHadronPt + max(0., pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5*pfIsolationR04().sumPUPt))/pt()'
trackIsoDef = 'trackIso()/pt()'
options['HLTProcessName']          = "HLT"
options['MUON_COLL']               = "slimmedMuons"
options['MUON_CUTS']               = "((isTrackerMuon || isGlobalMuon) && abs(eta)<2.4 && pt>5)"
options['MUON_TAG_CUTS']           = "(userInt('isTightMuon')==1 && pt > 25 && abs(eta) < 2.1 && "+isolationDef+" < 0.15)"
options['MAXEVENTS']               = cms.untracked.int32(-1) 
options['OUTPUTEDMFILENAME']       = 'edmFile.root'
options['DEBUG']                   = cms.bool(False)
options['TAU_COLL']                = "slimmedTaus"
options['TAU_CUTS']                = "(tauID(\"decayModeFinding\") && abs(eta)<2.3 && pt>17)"

from PhysicsTools.TagAndProbe.treeMakerOptions_cfi import *

if (varOptions.isMC):
    #options['INPUT_FILE_NAME']     = '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/00000/00F0B3DC-211B-E611-A6A0-001E67248A39.root'
    options['INPUT_FILE_NAME']     = '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/000FF6AC-9F2A-E611-A063-0CC47A4C8EB0.root'
    options['OUTPUT_FILE_NAME']    = "TnPTree_mc_muon.root"
    options['TnPPATHS']            = cms.vstring()#"HLT_IsoTkMu20_v*")
    options['TnPHLTTagFilters']    = cms.vstring()#"hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring()
    options['GLOBALTAG']           = 'auto:run2_mc'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()
else:
    options['INPUT_FILE_NAME']     = "/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/02D9C19F-571A-E611-AD8E-02163E013732.root"
    options['OUTPUT_FILE_NAME']    = "TnPTree_data_muon.root"
    options['TnPPATHS']            = ["HLT_IsoTkMu22_v*",]
    options['TnPHLTTagFilters']    = ["hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09"]
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("")
    options['GLOBALTAG']           = 'auto:run2_data'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()

###################################################################

#setModules(process, options)
from PhysicsTools.TagAndProbe.treeContent_cfi import *

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options['GLOBALTAG'], '')


process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options['INPUT_FILE_NAME']),
                            eventsToProcess = options['EVENTSToPROCESS']
                            )

process.maxEvents = cms.untracked.PSet( input = options['MAXEVENTS'])

process.sampleInfo = cms.EDProducer("tnp::SampleInfoTree",
                                        #vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                        genInfo = cms.InputTag("generator")
                                        )


##########
### ID ###
##########
process.mID = cms.EDProducer(
    "MuonIdEmbedder",
    src = cms.InputTag(options['MUON_COLL']),
    vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
)
process.mPV = cms.EDProducer(
    "MuonIpEmbedder",
    src = cms.InputTag('mID'),
    vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamspotSrc = cms.InputTag("offlineBeamSpot"),
)

process.idEmbedSequence = cms.Sequence(process.mID*process.mPV)
muonSource = 'mPV'


############
### Tags ###
############
process.tagMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag(muonSource),
    cut = cms.string(options['MUON_TAG_CUTS']),
    filter = cms.bool(True)
)

process.tagMuonsTriggerMatched = cms.EDProducer("PatMuonTriggerCandProducer",
    filterNames = cms.vstring(options['TnPHLTTagFilters']),
    inputs      = cms.InputTag("tagMuons"),
    bits        = cms.InputTag('TriggerResults::HLT'),
    objects     = cms.InputTag('selectedPatTrigger'),
    dR          = cms.double(0.1),
    isAND       = cms.bool(True)
    )

process.probeMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag(muonSource),
    cut = cms.string(options['MUON_CUTS']), 
)

process.probeTaus = cms.EDFilter("PATTauRefSelector",
    src = cms.InputTag(options['TAU_COLL']),
    cut = cms.string(options['TAU_CUTS']), 
)

######################
### Trigger Probes ###
######################

# triggers
#IsoMu19_eta2p1_LooseIsoPFTau20
#IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1
#IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg
#IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1
#IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg
#IsoMu22
#IsoMu24
#IsoMu27
#IsoTkMu22
#IsoTkMu24
#IsoTkMu27
#Mu17_TrkIsoVVL_Mu8_TrkIsoVVL
#Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
#Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL
#Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ
#Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
#Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL
#Mu300
#Mu350
#Mu45_eta2p1
#Mu50
#Mu55
#Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL

process.probeTriggerSeq = cms.Sequence()

process.probeTriggersMu17Leg = cms.EDProducer("PatMuonTriggerCandProducer",
    filterNames = cms.vstring("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4", "hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17"),
    inputs      = cms.InputTag("probeMuons"),
    bits        = cms.InputTag('TriggerResults::HLT'),
    objects     = cms.InputTag('selectedPatTrigger'),
    dR          = cms.double(0.1),
    isAND       = cms.bool(True)
    )
process.probeTriggerSeq += process.probeTriggersMu17Leg

process.tagTriggersMu17Leg = process.probeTriggersMu17Leg.clone()
process.tagTriggersMu17Leg.inputs = cms.InputTag("tagMuonsTriggerMatched")
process.probeTriggerSeq += process.tagTriggersMu17Leg

process.probeTriggersMu8Leg = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu8Leg.filterNames = cms.vstring("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4", "hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8")
process.probeTriggerSeq += process.probeTriggersMu8Leg

process.probeTriggersTkMu8Leg = process.probeTriggersMu17Leg.clone()
process.probeTriggersTkMu8Leg.filterNames = cms.vstring("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4", "hltDiMuonGlbFiltered17TrkFiltered8")
process.probeTriggerSeq += process.probeTriggersTkMu8Leg

process.probeTriggersMu8ORTkMu8LegPre = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu8ORTkMu8LegPre.filterNames = cms.vstring("hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8", "hltDiMuonGlbFiltered17TrkFiltered8")
process.probeTriggersMu8ORTkMu8LegPre.isAND = cms.bool(False)
process.probeTriggerSeq += process.probeTriggersMu8ORTkMu8LegPre

process.probeTriggersMu8ORTkMu8Leg = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu8ORTkMu8Leg.inputs = cms.InputTag('probeTriggersMu8ORTkMu8LegPre')
process.probeTriggersMu8ORTkMu8Leg.filterNames = cms.vstring("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4", "hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4")
process.probeTriggersMu8ORTkMu8Leg.isAND = cms.bool(False)
process.probeTriggerSeq += process.probeTriggersMu8ORTkMu8Leg

process.probeTriggersMu8LegDZ = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu8LegDZ.filterNames = cms.vstring("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2")
process.probeTriggerSeq += process.probeTriggersMu8LegDZ

process.probeTriggersTkMu8LegDZ = process.probeTriggersMu17Leg.clone()
process.probeTriggersTkMu8LegDZ.filterNames = cms.vstring("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2")
process.probeTriggerSeq += process.probeTriggersTkMu8LegDZ

process.probeTriggersMu8ORTkMu8LegDZ = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu8ORTkMu8LegDZ.filterNames = cms.vstring("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2", "hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2")
process.probeTriggersMu8ORTkMu8LegDZ.isAND = cms.bool(False)
process.probeTriggerSeq += process.probeTriggersMu8ORTkMu8LegDZ

# electron muon
process.probeTriggersMu17Ele12MLeg = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu17Ele12MLeg.filterNames = cms.vstring("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23")
process.probeTriggerSeq += process.probeTriggersMu17Ele12MLeg

process.probeTriggersMu8Ele17MLeg = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu8Ele17MLeg.filterNames = cms.vstring("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23")
process.probeTriggerSeq += process.probeTriggersMu8Ele17MLeg

process.probeTriggersMu8Ele23MLeg = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu8Ele23MLeg.filterNames = cms.vstring("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23")
process.probeTriggerSeq += process.probeTriggersMu8Ele23MLeg

process.probeTriggersMu23Ele8MLeg = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu23Ele8MLeg.filterNames = cms.vstring("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23")
process.probeTriggerSeq += process.probeTriggersMu23Ele8MLeg

process.probeTriggersMu23Ele12MLeg = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu23Ele12MLeg.filterNames = cms.vstring("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23")
process.probeTriggerSeq += process.probeTriggersMu23Ele12MLeg

# muon tau
process.probeTriggersMu17Tau20MLegSingleL1 = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu17Tau20MLegSingleL1.filterNames = cms.vstring("hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersMu17Tau20MLegSingleL1
process.tagTriggersMu17Tau20MLegSingleL1 = process.probeTriggersMu17Tau20MLegSingleL1.clone()
process.tagTriggersMu17Tau20MLegSingleL1.inputs = cms.InputTag("tagMuonsTriggerMatched")
process.probeTriggerSeq += process.tagTriggersMu17Tau20MLegSingleL1

process.probeTriggersMu19Tau20MLegSingleL1 = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu19Tau20MLegSingleL1.filterNames = cms.vstring("hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersMu19Tau20MLegSingleL1
process.tagTriggersMu19Tau20MLegSingleL1 = process.probeTriggersMu19Tau20MLegSingleL1.clone()
process.tagTriggersMu19Tau20MLegSingleL1.inputs = cms.InputTag("tagMuonsTriggerMatched")
process.probeTriggerSeq += process.tagTriggersMu19Tau20MLegSingleL1

process.probeTriggersMu21Tau20MLegSingleL1 = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu21Tau20MLegSingleL1.filterNames = cms.vstring("hltL3fL1sSingleMu20erIorSingleMu22erL1f0L2f10QL3Filtered21Q")
process.probeTriggerSeq += process.probeTriggersMu21Tau20MLegSingleL1
process.tagTriggersMu21Tau20MLegSingleL1 = process.probeTriggersMu21Tau20MLegSingleL1.clone()
process.tagTriggersMu21Tau20MLegSingleL1.inputs = cms.InputTag("tagMuonsTriggerMatched")
process.probeTriggerSeq += process.tagTriggersMu21Tau20MLegSingleL1

process.probeTriggersMu17Tau20MLeg = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu17Tau20MLeg.filterNames = cms.vstring("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersMu17Tau20MLeg
process.tagTriggersMu17Tau20MLeg = process.probeTriggersMu17Tau20MLeg.clone()
process.tagTriggersMu17Tau20MLeg.inputs = cms.InputTag("tagMuonsTriggerMatched")
process.probeTriggerSeq += process.tagTriggersMu17Tau20MLeg

process.probeTriggersMu19Tau20MLeg = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu19Tau20MLeg.filterNames = cms.vstring("hltL3crIsoL1sMu18erTauJet20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersMu19Tau20MLeg
process.tagTriggersMu19Tau20MLeg = process.probeTriggersMu19Tau20MLeg.clone()
process.tagTriggersMu19Tau20MLeg.inputs = cms.InputTag("tagMuonsTriggerMatched")
process.probeTriggerSeq += process.tagTriggersMu19Tau20MLeg


# IsoMu
process.probeTriggersIsoMu18 = process.probeTriggersMu17Leg.clone()
process.probeTriggersIsoMu18.filterNames = cms.vstring("hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoMu18

process.probeTriggersIsoMu20 = process.probeTriggersMu17Leg.clone()
process.probeTriggersIsoMu20.filterNames = cms.vstring("hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoMu20

process.probeTriggersIsoMu22 = process.probeTriggersMu17Leg.clone()
process.probeTriggersIsoMu22.filterNames = cms.vstring("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoMu22

process.probeTriggersIsoMu22Eta2p1 = process.probeTriggersMu17Leg.clone()
process.probeTriggersIsoMu22Eta2p1.filterNames = cms.vstring("hltL3crIsoL1sSingleMu20erL1f0L2f10QL3f22QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoMu22Eta2p1

process.probeTriggersIsoMu24 = process.probeTriggersMu17Leg.clone()
process.probeTriggersIsoMu24.filterNames = cms.vstring("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoMu24

process.probeTriggersIsoMu27 = process.probeTriggersMu17Leg.clone()
process.probeTriggersIsoMu27.filterNames = cms.vstring("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoMu27

# IsoTkMu
process.probeTriggersIsoTkMu18 = process.probeTriggersMu17Leg.clone()
process.probeTriggersIsoTkMu18.filterNames = cms.vstring("hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoTkMu18

process.probeTriggersIsoTkMu20 = process.probeTriggersMu17Leg.clone()
process.probeTriggersIsoTkMu20.filterNames = cms.vstring("hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoTkMu20

process.probeTriggersIsoTkMu22 = process.probeTriggersMu17Leg.clone()
process.probeTriggersIsoTkMu22.filterNames = cms.vstring("hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoTkMu22

process.probeTriggersIsoTkMu22Eta2p1 = process.probeTriggersMu17Leg.clone()
process.probeTriggersIsoTkMu22Eta2p1.filterNames = cms.vstring("hltL3fL1sMu20erL1f0Tkf22QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoTkMu22Eta2p1

process.probeTriggersIsoTkMu24 = process.probeTriggersMu17Leg.clone()
process.probeTriggersIsoTkMu24.filterNames = cms.vstring("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoTkMu24

process.probeTriggersIsoTkMu27 = process.probeTriggersMu17Leg.clone()
process.probeTriggersIsoTkMu27.filterNames = cms.vstring("hltL3fL1sMu22Or25L1f0Tkf27QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoTkMu27

# Mu
process.probeTriggersMu45Eta2p1 = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu45Eta2p1.filterNames = cms.vstring("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered45e2p1Q")
process.probeTriggerSeq += process.probeTriggersMu45Eta2p1

process.probeTriggersMu50 = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu50.filterNames = cms.vstring("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q")
process.probeTriggerSeq += process.probeTriggersMu50

process.probeTriggersMu55 = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu55.filterNames = cms.vstring("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered55Q")
process.probeTriggerSeq += process.probeTriggersMu55

# ORs
process.probeTriggersIsoMu18ORIsoTkMu18 = process.probeTriggersMu8ORTkMu8Leg.clone()
process.probeTriggersIsoMu18ORIsoTkMu18.filterNames = cms.vstring("hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09","hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoMu18ORIsoTkMu18

process.probeTriggersIsoMu20ORIsoTkMu20 = process.probeTriggersMu8ORTkMu8Leg.clone()
process.probeTriggersIsoMu20ORIsoTkMu20.filterNames = cms.vstring("hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09","hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoMu20ORIsoTkMu20

process.probeTriggersIsoMu22ORIsoTkMu22 = process.probeTriggersMu8ORTkMu8Leg.clone()
process.probeTriggersIsoMu22ORIsoTkMu22.filterNames = cms.vstring("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09","hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoMu22ORIsoTkMu22

process.probeTriggersIsoMu22Eta2p1ORIsoTkMu22Eta2p1 = process.probeTriggersMu8ORTkMu8Leg.clone()
process.probeTriggersIsoMu22Eta2p1ORIsoTkMu22Eta2p1.filterNames = cms.vstring("hltL3crIsoL1sSingleMu20erL1f0L2f10QL3f22QL3trkIsoFiltered0p09","hltL3fL1sMu20erL1f0Tkf22QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoMu22Eta2p1ORIsoTkMu22Eta2p1

process.probeTriggersIsoMu24ORIsoTkMu24 = process.probeTriggersMu8ORTkMu8Leg.clone()
process.probeTriggersIsoMu24ORIsoTkMu24.filterNames = cms.vstring("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09","hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoMu24ORIsoTkMu24

process.probeTriggersIsoMu27ORIsoTkMu27 = process.probeTriggersMu8ORTkMu8Leg.clone()
process.probeTriggersIsoMu27ORIsoTkMu27.filterNames = cms.vstring("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09","hltL3fL1sMu22Or25L1f0Tkf27QL3trkIsoFiltered0p09")
process.probeTriggerSeq += process.probeTriggersIsoMu27ORIsoTkMu27

# "soups"
# isomu22 isotkmu22 45eta2p1 50
process.probeTriggersSingleMuSoup = process.probeTriggersMu8ORTkMu8Leg.clone()
process.probeTriggersSingleMuSoup.filterNames = cms.vstring("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09","hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09","hltL3fL1sMu22Or25L1f0L2f10QL3Filtered45e2p1Q","hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q")
process.probeTriggerSeq += process.probeTriggersSingleMuSoup

############
### Taus ###
############
process.probeDummy = cms.EDProducer("PatTauTriggerCandProducer",
    filterNames = cms.vstring(),
    inputs      = cms.InputTag("probeTaus"),
    bits        = cms.InputTag('TriggerResults::HLT'),
    objects     = cms.InputTag('selectedPatTrigger'),
    dR          = cms.double(0.1),
    isAND       = cms.bool(True)
    )


process.probeTriggersMu17Tau20TLegSingleL1 = process.probeDummy.clone()
process.probeTriggersMu17Tau20TLegSingleL1.filterNames = cms.vstring("hltOverlapFilterSingleIsoMu17LooseIsoPFTau20")
process.probeTriggerSeq += process.probeTriggersMu17Tau20TLegSingleL1

process.probeTriggersMu19Tau20TLegSingleL1 = process.probeDummy.clone()
process.probeTriggersMu19Tau20TLegSingleL1.filterNames = cms.vstring("hltOverlapFilterSingleIsoMu19LooseIsoPFTau20")
process.probeTriggerSeq += process.probeTriggersMu19Tau20TLegSingleL1

process.probeTriggersMu21Tau20TLegSingleL1 = process.probeDummy.clone()
process.probeTriggersMu21Tau20TLegSingleL1.filterNames = cms.vstring("hltOverlapFilterSingleIsoMu21LooseIsoPFTau20")
process.probeTriggerSeq += process.probeTriggersMu21Tau20TLegSingleL1

process.probeTriggersMu17Tau20TLeg = process.probeDummy.clone()
process.probeTriggersMu17Tau20TLeg.filterNames = cms.vstring("hltOverlapFilterIsoMu17LooseIsoPFTau20")
process.probeTriggerSeq += process.probeTriggersMu17Tau20TLeg

process.probeTriggersMu19Tau20TLeg = process.probeDummy.clone()
process.probeTriggersMu19Tau20TLeg.filterNames = cms.vstring("hltOverlapFilterIsoMu19LooseIsoPFTau20")
process.probeTriggerSeq += process.probeTriggersMu19Tau20TLeg


###################################################################
## TnP PAIRS
###################################################################

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuonsTriggerMatched@+ probeMuons@-"), # charge coniugate states are implied
    cut   = cms.string("40 < mass < 200")
)

process.tpPairsTau = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuonsTriggerMatched@+ probeTaus@-"), # charge coniugate states are implied
    cut   = cms.string("20 < mass < 200")
)

#process.tpPairsMCEmbedded = cms.EDProducer("pairMCInfoEmbedder",
#    input = cms.InputTag("tpPairs"),
#    leg1Matches = cms.InputTag("muMcMatch"),
#    leg2Matches = cms.InputTag("muMcMatch")
#)

process.muMcMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag(muonSource),
    distMin = cms.double(0.3),
    matched = cms.InputTag("prunedGenParticles"),
    checkCharge = cms.bool(True)
)

##############
### Pileup ###
##############
from SimGeneral.MixingModule.mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU_cfi import mix
from DevTools.TagAndProbe.pileup_cfi import currentPileup
pu_distribs = { "mc" : mix.input.nbPileupEvents.probValue }
data_pu_distribs = { "golden" : currentPileup }

process.pileupReweightingProducer = cms.EDProducer("PileupWeightProducer",
                                                   #hardcodedWeights = cms.untracked.bool(True),
                                                   pileupInfoTag    = cms.InputTag("slimmedAddPileupInfo"),
                                                   PileupMC = cms.vdouble(pu_distribs["mc"]),
                                                   PileupData = cms.vdouble(data_pu_distribs["golden"]),
                                                   )


##########################################################################
## TREE MAKER OPTIONS
##########################################################################
ZVariablesToStore = cms.PSet(
    eta = cms.string("eta"),
    abseta = cms.string("abs(eta)"),
    pt  = cms.string("pt"),
    mass  = cms.string("mass"),
    )   

ProbeVariablesToStore = cms.PSet(
    probe_eta             = cms.string("eta"),
    probe_abseta          = cms.string("abs(eta)"),
    probe_pt              = cms.string("pt"),
    probe_et              = cms.string("et"),
    probe_e               = cms.string("energy"),
    probe_q               = cms.string("charge"),
    probe_isoR04          = cms.string(isolationDef),
    probe_trackIso        = cms.string(trackIsoDef),
    probe_dz              = cms.string('userFloat("dz")'),
    probe_dxy             = cms.string('userFloat("dB2D")'),
    probe_isGlobalMuon    = cms.string('isGlobalMuon'),
    probe_isTrackerMuon   = cms.string('isTrackerMuon'),
    probe_matchedStations = cms.string('numberOfMatchedStations'),
    probe_bestTrackType   = cms.string('muonBestTrackType'),
    )

ProbeVariablesToStoreTau = cms.PSet(
    probe_eta             = cms.string("eta"),
    probe_abseta          = cms.string("abs(eta)"),
    probe_pt              = cms.string("pt"),
    probe_et              = cms.string("et"),
    probe_e               = cms.string("energy"),
    probe_q               = cms.string("charge"),
    )

TagFlagsToStore = cms.PSet(
    passingMu17                 = cms.InputTag("tagTriggersMu17Leg"),
    pasingMu17Tau20MLegSingleL1 = cms.InputTag("tagTriggersMu17Tau20MLegSingleL1"),
    pasingMu19Tau20MLegSingleL1 = cms.InputTag("tagTriggersMu19Tau20MLegSingleL1"),
    pasingMu21Tau20MLegSingleL1 = cms.InputTag("tagTriggersMu21Tau20MLegSingleL1"),
    pasingMu17Tau20MLeg         = cms.InputTag("tagTriggersMu17Tau20MLeg"),
    pasingMu19Tau20MLeg         = cms.InputTag("tagTriggersMu19Tau20MLeg"),
)

TagVariablesToStore = cms.PSet(
    eta    = cms.string("eta"),
    abseta = cms.string("abs(eta)"),
    pt     = cms.string("pt"),
    et     = cms.string("et"),
    e      = cms.string("energy"),
    q      = cms.string("charge"),
    )

CommonStuffForMuonProbe = cms.PSet(
    variables = cms.PSet(ProbeVariablesToStore),
    ignoreExceptions =  cms.bool (True),
    addRunLumiInfo   =  cms.bool (True),
    pileupInfoTag = cms.InputTag("slimmedAddPileupInfo"),
    addEventVariablesInfo   =  cms.bool(True),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    #pfMet = cms.InputTag(""),
    pairVariables =  cms.PSet(ZVariablesToStore),
    pairFlags     =  cms.PSet(
        mass60to120 = cms.string("60 < mass < 120")
        ),
    tagVariables   =  cms.PSet(TagVariablesToStore),
    tagFlags       =  cms.PSet(TagFlagsToStore),    
    )

CommonStuffForTauProbe = cms.PSet(
    variables = cms.PSet(ProbeVariablesToStoreTau),
    ignoreExceptions =  cms.bool (True),
    addRunLumiInfo   =  cms.bool (True),
    pileupInfoTag = cms.InputTag("slimmedAddPileupInfo"),
    addEventVariablesInfo   =  cms.bool(True),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    #pfMet = cms.InputTag(""),
    pairVariables =  cms.PSet(ZVariablesToStore),
    pairFlags     =  cms.PSet(
        mass60to120 = cms.string("60 < mass < 120")
        ),
    tagVariables   =  cms.PSet(TagVariablesToStore),
    tagFlags       =  cms.PSet(TagFlagsToStore),    
    )

#mcTruthCommonStuff = cms.PSet(
#    isMC = cms.bool(False),
#    tagMatches = cms.InputTag("muMcMatch"),
#    probeMatches = cms.InputTag("muMcMatch"),
#    motherPdgId = cms.vint32(22,23),
#    #motherPdgId = cms.vint32(443), # JPsi
#    #motherPdgId = cms.vint32(553), # Yupsilon
#    makeMCUnbiasTree = cms.bool(False),
#    checkMotherInUnbiasEff = cms.bool(False),
#    mcVariables = cms.PSet(
#        probe_eta = cms.string("eta"),
#        probe_abseta = cms.string("abs(eta)"),
#        probe_et  = cms.string("et"),
#        probe_e  = cms.string("energy"),
#        ),
#    mcFlags     =  cms.PSet(
#        probe_isPromptFinalState = cms.string("isPromptFinalState")
#        ),      
#    )
mcTruthCommonStuff = cms.PSet(
    isMC = cms.bool(True),
    #tagMatches = cms.InputTag("McMatchTag"),
    #motherPdgId = cms.vint32(),
    motherPdgId = cms.vint32(22,23),
    #motherPdgId = cms.vint32(443), # JPsi
    #motherPdgId = cms.vint32(553), # Yupsilon
    #makeMCUnbiasTree = cms.bool(False),
    #checkMotherInUnbiasEff = cms.bool(False),
    genParticles = cms.InputTag("prunedGenParticles"),
    useTauDecays = cms.bool(False),
    checkCharge = cms.bool(False),
    pdgId = cms.int32(13),
    mcVariables = cms.PSet(
        probe_eta = cms.string("eta"),
        probe_abseta = cms.string("abs(eta)"),
        probe_et  = cms.string("et"),
        probe_e  = cms.string("energy"),
        ),
    mcFlags     =  cms.PSet(
        probe_flag = cms.string("pt>0")
        ),
    )

if (not varOptions.isMC):
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(False)
        )


process.muonEffs = cms.EDAnalyzer("TagProbeFitTreeProducer",
    CommonStuffForMuonProbe, mcTruthCommonStuff,
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("Random2"),
    flags         = cms.PSet(
        passingLoose                          = cms.string("isLooseMuon"),
        passingMedium                         = cms.string("isMediumMuon"),
        passingTight                          = cms.string("userInt('isTightMuon')==1"), 
        passingMediumICHEP                    = cms.string("userInt('isMediumMuonICHEP')==1"), 
        passingIsoVeryLoose                   = cms.string(isolationDef+" < 0.4"),
        passingIsoLoose                       = cms.string(isolationDef+" < 0.25"),
        passingIsoTight                       = cms.string(isolationDef+" < 0.15"),
        passingTrackIso                       = cms.string(trackIsoDef+" < 0.4"),
        passingMu17                           = cms.InputTag("probeTriggersMu17Leg"),
        passingMu8                            = cms.InputTag("probeTriggersMu8Leg"),
        passingTkMu8                          = cms.InputTag("probeTriggersTkMu8Leg"),
        passingMu8ORTkMu8                     = cms.InputTag("probeTriggersMu8ORTkMu8Leg"),
        passingMu8DZ                          = cms.InputTag("probeTriggersMu8LegDZ"),
        passingTkMu8DZ                        = cms.InputTag("probeTriggersTkMu8LegDZ"),
        passingMu8ORTkMu8DZ                   = cms.InputTag("probeTriggersMu8ORTkMu8LegDZ"),
        #passingIsoMu18                        = cms.InputTag("probeTriggersIsoMu18"),
        #passingIsoMu20                        = cms.InputTag("probeTriggersIsoMu20"),
        passingIsoMu22                        = cms.InputTag("probeTriggersIsoMu22"),
        passingIsoMu24                        = cms.InputTag("probeTriggersIsoMu24"),
        passingIsoMu27                        = cms.InputTag("probeTriggersIsoMu27"),
        #passingIsoTkMu18                      = cms.InputTag("probeTriggersIsoTkMu18"),
        #passingIsoTkMu20                      = cms.InputTag("probeTriggersIsoTkMu20"),
        passingIsoTkMu22                      = cms.InputTag("probeTriggersIsoTkMu22"),
        passingIsoTkMu24                      = cms.InputTag("probeTriggersIsoTkMu24"),
        passingIsoTkMu27                      = cms.InputTag("probeTriggersIsoTkMu27"),
        #passingIsoMu18ORIsoTkMu18             = cms.InputTag("probeTriggersIsoMu18ORIsoTkMu18"),
        #passingIsoMu20ORIsoTkMu20             = cms.InputTag("probeTriggersIsoMu20ORIsoTkMu20"),
        passingIsoMu22ORIsoTkMu22             = cms.InputTag("probeTriggersIsoMu22ORIsoTkMu22"),
        passingIsoMu24ORIsoTkMu24             = cms.InputTag("probeTriggersIsoMu24ORIsoTkMu24"),
        passingIsoMu27ORIsoTkMu27             = cms.InputTag("probeTriggersIsoMu27ORIsoTkMu27"),
        passingMu45Eta2p1                     = cms.InputTag("probeTriggersMu45Eta2p1"),
        passingMu50                           = cms.InputTag("probeTriggersMu50"),
        passingMu55                           = cms.InputTag("probeTriggersMu55"),
        #passingMu17Ele12MLeg                  = cms.InputTag("probeTriggersMu17Ele12MLeg"),
        #passingMu8Ele17MLeg                   = cms.InputTag("probeTriggersMu8Ele17MLeg"),
        #passingMu8Ele23MLeg                   = cms.InputTag("probeTriggersMu8Ele23MLeg"),
        #passingMu23Ele8MLeg                   = cms.InputTag("probeTriggersMu23Ele8MLeg"),
        #passingMu23Ele12MLeg                  = cms.InputTag("probeTriggersMu23Ele12MLeg"),
        #passingMu17Tau20MLegSingleL1          = cms.InputTag("probeTriggersMu17Tau20MLegSingleL1"),
        passingMu19Tau20MLegSingleL1          = cms.InputTag("probeTriggersMu19Tau20MLegSingleL1"),
        passingMu21Tau20MLegSingleL1          = cms.InputTag("probeTriggersMu21Tau20MLegSingleL1"),
        passingMu17Tau20MLeg                  = cms.InputTag("probeTriggersMu17Tau20MLeg"),
        passingMu19Tau20MLeg                  = cms.InputTag("probeTriggersMu19Tau20MLeg"),
        passingSingleMuSoup                   = cms.InputTag("probeTriggersSingleMuSoup"),
    ),
    allProbes     = cms.InputTag("probeMuons"),
    )

process.tauEffs = cms.EDAnalyzer("TagProbeFitTreeProducer",
    CommonStuffForTauProbe, mcTruthCommonStuff,
    tagProbePairs = cms.InputTag("tpPairsTau"),
    arbitration   = cms.string("Random2"),
    flags         = cms.PSet(
        passingMu17Tau20TLegSingleL1          = cms.InputTag("probeTriggersMu17Tau20TLegSingleL1"),
        passingMu19Tau20TLegSingleL1          = cms.InputTag("probeTriggersMu19Tau20TLegSingleL1"),
        passingMu21Tau20TLegSingleL1          = cms.InputTag("probeTriggersMu21Tau20TLegSingleL1"),
        passingMu17Tau20TLeg                  = cms.InputTag("probeTriggersMu17Tau20TLeg"),
        passingMu19Tau20TLeg                  = cms.InputTag("probeTriggersMu19Tau20TLeg"),
    ),
    allProbes     = cms.InputTag("probeTaus"),
    )

process.tpPairSeq = cms.Sequence(
    process.tpPairs*
    process.tpPairsTau
)

if varOptions.isMC :
    process.tpPairSeq += process.muMcMatch
    #process.tpPairSeq += process.tpPairsMCEmbedded
    process.tpPairSeq += process.pileupReweightingProducer
    process.muonEffs.isMC = cms.bool(True)
    process.muonEffs.eventWeight   = cms.InputTag("generator")
    process.muonEffs.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")
    #setattr(process.muonEffs.pairVariables, 'mc_mass', cms.string("userFloat('mc_mass')"))
    #process.muonEffs.tagProbePairs = cms.InputTag("tpPairs")
    process.tauEffs.isMC = cms.bool(True)
    process.tauEffs.eventWeight   = cms.InputTag("generator")
    process.tauEffs.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")

#if not options.isMC :
#    import FWCore.PythonUtilities.LumiList as LumiList
#    process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/'+options['json']).getVLuminosityBlockRange()

process.p = cms.Path(
    process.sampleInfo *
    process.idEmbedSequence *
    (process.tagMuons + process.probeMuons + process.probeTaus) *
    (process.tagMuonsTriggerMatched + process.probeTriggerSeq) *
    process.tpPairSeq *
    process.muonEffs *
    process.tauEffs
    )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(options['OUTPUTEDMFILENAME']),
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
                               )
process.outpath = cms.EndPath(process.out)
if (not options['DEBUG']):
    process.outpath.remove(process.out)

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string(options['OUTPUT_FILE_NAME']),
    closeFileFast = cms.untracked.bool(True)
    )

