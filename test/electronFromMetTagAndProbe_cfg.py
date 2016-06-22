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

options['isMC']                    = varOptions.isMC
options['HLTProcessName']          = "HLT"
options['ELECTRON_COLL']           = "slimmedElectrons"
options['ELECTRON_CUTS']           = "(abs(eta)<2.5)"
options['ELECTRON_TAG_CUTS']       = "(abs(eta)<=2.5) && !(1.4442<=abs(eta)<=1.566) && pt >= 25.0"
options['SUPERCLUSTER_COLL']       = "reducedEgamma:reducedSuperClusters"
options['SUPERCLUSTER_CUTS']       = "abs(eta)<2.5 && !(1.4442< abs(eta) <1.566) && et>10.0"
options['MAXEVENTS']               = cms.untracked.int32(-1) 
options['OUTPUTEDMFILENAME']       = 'edmFile.root'
options['DEBUG']                   = cms.bool(False)
options['useAOD']                  = cms.bool(False)

from PhysicsTools.TagAndProbe.treeMakerOptions_cfi import *

if (varOptions.isMC):
    options['INPUT_FILE_NAME']     = '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/00000/00F0B3DC-211B-E611-A6A0-001E67248A39.root'
    options['OUTPUT_FILE_NAME']    = "TnPTree_mc_singleElectron.root"
    options['TnPPATHS']            = cms.vstring()#"HLT_Ele23_WPLoose_Gsf_v*")
    options['TnPHLTTagFilters']    = cms.vstring()#"hltEle23WPLooseGsfTrackIsoFilter")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring()
    options['GLOBALTAG']           = 'auto:run2_mc'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()
else:
    options['INPUT_FILE_NAME']     = "/store/data/Run2016B/MET/MINIAOD/PromptReco-v2/000/273/158/00000/06A9DFDA-201A-E611-858F-02163E0136F7.root"
    options['OUTPUT_FILE_NAME']    = "TnPTree_data_singleElectron.root"
    options['TnPPATHS']            = ["HLT_PFMET170_HBHECleaned_v*","HLT_MET200_v*"]
    options['TnPHLTTagFilters']    = []
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("")
    options['GLOBALTAG']           = 'auto:run2_data'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()

###################################################################

#setModules(process, options)

#############
### Setup ###
#############
process.sampleInfo = cms.EDProducer("tnp::SampleInfoTree",
                                    #vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    genInfo = cms.InputTag("generator")
                                    )

process.eleVarHelper = cms.EDProducer("PatElectronVariableHelper",
                                      probes = cms.InputTag(options['ELECTRON_COLL']),
                                      vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
                                      )

##################
### Select HLT ###
##################
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone()
process.hltFilter.throw = cms.bool(True)
process.hltFilter.HLTPaths = options['TnPPATHS']

##############
### Pileup ###
##############
from SimGeneral.MixingModule.mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU_cfi import mix
from DevTools.TagAndProbe.pileup_cfi import currentPileup
pu_distribs = { "mc" : mix.input.nbPileupEvents.probValue }
data_pu_distribs = { "golden" : currentPileup }

from PhysicsTools.TagAndProbe.pileupConfiguration_cfi import pileupProducer
process.pileupReweightingProducer = pileupProducer.clone()

process.pileupReweightingProducer.pileupMC = cms.vdouble(pu_distribs['mc'])
process.pileupReweightingProducer.PileupData = cms.vdouble(data_pu_distribs["golden"])

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

##########
### ID ###
##########

from PhysicsTools.TagAndProbe.electronIDModules_cfi import *
setIDs(process, options)

############
### Tags ###
############
process.tagElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag(options['ELECTRON_COLL']),
    cut = cms.string(options['ELECTRON_TAG_CUTS']),
    filter = cms.bool(True)
)

process.tagElectronsTriggerMatched = cms.EDProducer("PatElectronTriggerCandProducer",
    filterNames = cms.vstring(options['TnPHLTTagFilters']),
    inputs      = cms.InputTag("tagElectrons"),
    bits        = cms.InputTag('TriggerResults::HLT'),
    objects     = cms.InputTag('selectedPatTrigger'),
    dR          = cms.double(0.4),
    isAND       = cms.bool(True)
    )

process.goodElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag(options['ELECTRON_COLL']),
    cut = cms.string(options['ELECTRON_CUTS']),
)



######################
### Trigger probes ###
######################
process.goodElectronsMeasureHLT = cms.Sequence()

# single electron
process.goodElectronsMeasureHLTEle23 = cms.EDProducer("PatElectronTriggerCandProducer",
                                                filterNames = cms.vstring("hltEle23WPLooseGsfTrackIsoFilter"),
                                                inputs      = cms.InputTag("goodElectrons"),
                                                bits        = cms.InputTag('TriggerResults::HLT'),
                                                objects     = cms.InputTag('selectedPatTrigger'),
                                                dR          = cms.double(0.3),
                                                isAND       = cms.bool(False)
                                                )
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle23

process.goodElectronsMeasureHLTEle22Eta2p1 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle22Eta2p1.filterNames = cms.vstring("hltSingleEle22WPLooseGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle22Eta2p1

process.goodElectronsMeasureHLTEle24Eta2p1 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle24Eta2p1.filterNames = cms.vstring("hltSingleEle24WPLooseGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle24Eta2p1

process.goodElectronsMeasureHLTEle25Eta2p1 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle25Eta2p1.filterNames = cms.vstring("hltEle25erWPLooseGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle25Eta2p1

process.goodElectronsMeasureHLTEle27 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle27.filterNames = cms.vstring("hltEle27noerWPLooseGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle27

process.goodElectronsMeasureHLTEle35 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle35.filterNames = cms.vstring("hltEle35WPLooseGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle35

# double electron
process.goodElectronsMeasureHLTEle17 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle17.filterNames = cms.vstring("hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle17

process.goodElectronsMeasureHLTEle12 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle12.filterNames = cms.vstring("hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle12

process.goodElectronsMeasureHLTEle17Ele12Leg1 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle17Ele12Leg1.filterNames = cms.vstring("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle17Ele12Leg1

#process.goodElectronsMeasureHLTEle17Ele12Leg1L1EG15 = cms.EDProducer("PatElectronL1CandProducer",
#        inputs = cms.InputTag("goodElectronsMeasureHLTEle17Ele12Leg1"),
#        isoObjects = cms.InputTag("l1extraParticles:Isolated"),
#        nonIsoObjects = cms.InputTag("l1extraParticles:NonIsolated"),
#        minET = cms.double(15.),
#        dRmatch = cms.double(.5)
#        )
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle17Ele12Leg1L1EG15

process.goodElectronsMeasureHLTEle17Ele12Leg2 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle17Ele12Leg2.filterNames = cms.vstring("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle17Ele12Leg2

process.goodElectronsMeasureHLTMu17Ele12ELeg = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTMu17Ele12ELeg.filterNames = cms.vstring("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTMu17Ele12ELeg

#################
### SEQUENCES ###
#################

process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag(options['ELECTRON_COLL'])
process.ele_sequence = cms.Sequence(
    process.goodElectrons +
    process.egmGsfElectronIDSequence +
    process.tagElectrons + 
    process.tagElectronsTriggerMatched +
    process.goodElectronsMeasureHLT
    )


#################
### TnP PAIRS ###
#################

process.allTagsAndProbes = cms.Sequence()

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagElectronsTriggerMatched@+ goodElectrons@-"), # charge coniugate states are implied
    cut   = cms.string("40 < mass < 200")
)

process.allTagsAndProbes *= process.tpPairs

##########################
### TREE MAKER OPTIONS ###
##########################
from PhysicsTools.TagAndProbe.treeContent_cfi import *

if (not varOptions.isMC):
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(False)
        )

# remove sc stuff that is broken
del CommonStuffForGsfElectronProbe.variables.probe_sc_energy
del CommonStuffForGsfElectronProbe.variables.probe_sc_et
del CommonStuffForGsfElectronProbe.variables.probe_sc_eta
del CommonStuffForGsfElectronProbe.variables.probe_sc_abseta
del CommonStuffForGsfElectronProbe.tagVariables.sc_energy
del CommonStuffForGsfElectronProbe.tagVariables.sc_et
del CommonStuffForGsfElectronProbe.tagVariables.sc_eta
del CommonStuffForGsfElectronProbe.tagVariables.sc_abseta
del CommonStuffForSuperClusterProbe.tagVariables.sc_energy
del CommonStuffForSuperClusterProbe.tagVariables.sc_et
del CommonStuffForSuperClusterProbe.tagVariables.sc_eta
del CommonStuffForSuperClusterProbe.tagVariables.sc_abseta

process.GsfElectronToTrigger = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                              CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
                                              tagProbePairs = cms.InputTag("tpPairs"),
                                              arbitration   = cms.string("Random2"),
                                              flags         = cms.PSet(
                                                  passingHLTEle22Eta2p1    = cms.InputTag("goodElectronsMeasureHLTEle22Eta2p1"),
                                                  passingHLTEle24Eta2p1    = cms.InputTag("goodElectronsMeasureHLTEle24Eta2p1"),
                                                  passingHLTEle25Eta2p1    = cms.InputTag("goodElectronsMeasureHLTEle25Eta2p1"),
                                                  passingHLTEle23    = cms.InputTag("goodElectronsMeasureHLTEle23"),
                                                  passingHLTEle27    = cms.InputTag("goodElectronsMeasureHLTEle27"),
                                                  passingHLTEle35    = cms.InputTag("goodElectronsMeasureHLTEle35"),
                                                  passingHLTEle17    = cms.InputTag("goodElectronsMeasureHLTEle17"),
                                                  passingHLTEle12    = cms.InputTag("goodElectronsMeasureHLTEle12"),
                                                  passingHLTEle17Ele12Leg1    = cms.InputTag("goodElectronsMeasureHLTEle17Ele12Leg1"),
                                                  #passingHLTEle17Ele12Leg1L1Match    = cms.InputTag("goodElectronsMeasureHLTEle17Ele12Leg1L1EG15"),
                                                  passingHLTEle17Ele12Leg2    = cms.InputTag("goodElectronsMeasureHLTEle17Ele12Leg2"),
                                                  passingHLTMu17Ele12ELeg     = cms.InputTag("goodElectronsMeasureHLTMu17Ele12ELeg"),
                                                                       ),                                               
                                              allProbes     = cms.InputTag("probeElectrons"),
                                              )

if (varOptions.isMC):
    #process.GsfElectronToTrigger.probeMatches  = cms.InputTag("McMatchHLT")
    process.GsfElectronToTrigger.eventWeight   = cms.InputTag("generator")
    process.GsfElectronToTrigger.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")


process.tree_sequence = cms.Sequence()
process.tree_sequence *= process.GsfElectronToTrigger

##########################################################################
## PATHS
##########################################################################

process.out = cms.OutputModule("PoolOutputModule", 
                               fileName = cms.untracked.string(options['OUTPUTEDMFILENAME']),
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
                               )
process.outpath = cms.EndPath(process.out)
if (not options['DEBUG']):
    process.outpath.remove(process.out)

if (varOptions.isMC):
    process.p = cms.Path(
        process.sampleInfo +
        process.hltFilter +
        process.ele_sequence + 
        process.allTagsAndProbes +
        process.pileupReweightingProducer +
        process.eleVarHelper +
        process.tree_sequence
        )
else:
    process.p = cms.Path(
        process.sampleInfo +
        process.hltFilter +
        process.ele_sequence + 
        process.allTagsAndProbes +
        process.eleVarHelper +
        process.tree_sequence
        )

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string(options['OUTPUT_FILE_NAME']),
    closeFileFast = cms.untracked.bool(True)
    )
