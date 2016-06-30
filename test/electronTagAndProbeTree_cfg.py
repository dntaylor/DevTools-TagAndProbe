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
options['ELECTRON_TAG_CUTS']       = "(abs(eta)<=2.5) && !(1.4442<=abs(eta)<=1.566) && pt >= 30.0"
options['SUPERCLUSTER_COLL']       = "reducedEgamma:reducedSuperClusters"
options['SUPERCLUSTER_CUTS']       = "abs(eta)<2.5 && !(1.4442< abs(eta) <1.566) && et>10.0"
options['MAXEVENTS']               = cms.untracked.int32(-1) 
options['useAOD']                  = cms.bool(False)
options['DOTRIGGER']               = cms.bool(True)
options['DORECO']                  = cms.bool(True)
options['DOID']                    = cms.bool(True)
options['OUTPUTEDMFILENAME']       = 'edmFile.root'
options['DEBUG']                   = cms.bool(False)

from PhysicsTools.TagAndProbe.treeMakerOptions_cfi import *

if (varOptions.isMC):
    #options['INPUT_FILE_NAME']     = '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/00000/00F0B3DC-211B-E611-A6A0-001E67248A39.root'
    options['INPUT_FILE_NAME']     = '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/000FF6AC-9F2A-E611-A063-0CC47A4C8EB0.root'
    options['OUTPUT_FILE_NAME']    = "TnPTree_mc_electron.root"
    options['TnPPATHS']            = cms.vstring()#"HLT_Ele23_WPLoose_Gsf_v*")
    options['TnPHLTTagFilters']    = cms.vstring()#"hltEle23WPLooseGsfTrackIsoFilter")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring()
    options['GLOBALTAG']           = 'auto:run2_mc'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()
else:
    options['INPUT_FILE_NAME']     = "/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/06277EC1-181A-E611-870F-02163E0145E5.root"
    options['OUTPUT_FILE_NAME']    = "TnPTree_data_electron.root"
    options['TnPPATHS']            = ["HLT_Ele25_WPTight_Gsf_v*",]
    options['TnPHLTTagFilters']    = ["hltEle25WPTightGsfTrackIsoFilter"]
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("")
    options['GLOBALTAG']           = 'auto:run2_data'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()

###################################################################

setModules(process, options)

# manually fix pileup
from SimGeneral.MixingModule.mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU_cfi import mix
from DevTools.TagAndProbe.pileup_cfi import currentPileup
pu_distribs = { "mc" : mix.input.nbPileupEvents.probValue }
data_pu_distribs = { "golden" : currentPileup }

process.pileupReweightingProducer.PileupMC = cms.vdouble(pu_distribs['mc'])
process.pileupReweightingProducer.PileupData = cms.vdouble(data_pu_distribs["golden"])

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

###################################################################
## ID
###################################################################

from PhysicsTools.TagAndProbe.electronIDModules_cfi import *
setIDs(process, options)

# trigger
#Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
#Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20
#Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1
#Ele25_eta2p1_WPTight_Gsf
#Ele27_WPTight_Gsf
#Ele27_eta2p1_WPLoose_Gsf
#Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1
#Ele27_eta2p1_WPTight_Gsf
#Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1
#Ele32_eta2p1_WPTight_Gsf
#Ele45_WPLoose_Gsf


###############
### Trigger ###
###############
process.goodElectronsMeasureHLT = cms.Sequence()

process.goodElectronsMeasureHLTEle23 = cms.EDProducer("PatElectronTriggerCandProducer",
                                                filterNames = cms.vstring("hltEle23WPLooseGsfTrackIsoFilter"),
                                                inputs      = cms.InputTag("goodElectronsProbeMeasureHLT"),
                                                bits        = cms.InputTag('TriggerResults::HLT'),
                                                objects     = cms.InputTag('selectedPatTrigger'),
                                                dR          = cms.double(0.1),
                                                isAND       = cms.bool(False)
                                                )
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle23

#process.goodElectronsMeasureHLTEle22Eta2p1 = process.goodElectronsMeasureHLTEle23.clone()
#process.goodElectronsMeasureHLTEle22Eta2p1.filterNames = cms.vstring("hltSingleEle22WPLooseGsfTrackIsoFilter")
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle22Eta2p1

#process.goodElectronsMeasureHLTEle24Eta2p1 = process.goodElectronsMeasureHLTEle23.clone()
#process.goodElectronsMeasureHLTEle24Eta2p1.filterNames = cms.vstring("hltSingleEle24WPLooseGsfTrackIsoFilter")
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle24Eta2p1

#process.goodElectronsMeasureHLTEle25Eta2p1 = process.goodElectronsMeasureHLTEle23.clone()
#process.goodElectronsMeasureHLTEle25Eta2p1.filterNames = cms.vstring("hltEle25erWPLooseGsfTrackIsoFilter")
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle25Eta2p1

process.goodElectronsMeasureHLTEle27Eta2p1 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle27Eta2p1.filterNames = cms.vstring("hltEle27erWPLooseGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle27Eta2p1

process.goodElectronsMeasureHLTEle25Eta2p1Tight = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle25Eta2p1Tight.filterNames = cms.vstring("hltEle25erWPTightGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle25Eta2p1Tight

process.goodElectronsMeasureHLTEle27Eta2p1Tight = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle27Eta2p1Tight.filterNames = cms.vstring("hltEle27erWPTightGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle27Eta2p1Tight

process.goodElectronsMeasureHLTEle32Eta2p1Tight = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle32Eta2p1Tight.filterNames = cms.vstring("hltEle32erWPTightGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle32Eta2p1Tight

#process.goodElectronsMeasureHLTEle27 = process.goodElectronsMeasureHLTEle23.clone()
#process.goodElectronsMeasureHLTEle27.filterNames = cms.vstring("hltEle27noerWPLooseGsfTrackIsoFilter")
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle27

#process.goodElectronsMeasureHLTEle35 = process.goodElectronsMeasureHLTEle23.clone()
#process.goodElectronsMeasureHLTEle35.filterNames = cms.vstring("hltEle35WPLooseGsfTrackIsoFilter")
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle35

process.goodElectronsMeasureHLTEle45 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle45.filterNames = cms.vstring("hltEle45WPLooseGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle45

#process.goodElectronsMeasureHLTEle25Tight = process.goodElectronsMeasureHLTEle23.clone()
#process.goodElectronsMeasureHLTEle25Tight.filterNames = cms.vstring("hltEle25WPTightGsfTrackIsoFilter")
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle25Tight

process.goodElectronsMeasureHLTEle27Tight = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle27Tight.filterNames = cms.vstring("hltEle27WPTightGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle27Tight

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


process.goodElectronsMeasureHLTEle17Ele12Leg2 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle17Ele12Leg2.filterNames = cms.vstring("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle17Ele12Leg2

process.goodElectronsMeasureHLTEle17Ele12DZ = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle17Ele12DZ.filterNames = cms.vstring("hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle17Ele12DZ

process.goodElectronsMeasureHLTEle23Ele12Leg1 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle23Ele12Leg1.filterNames = cms.vstring("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle23Ele12Leg1

process.goodElectronsMeasureHLTEle23Ele12Leg2 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle23Ele12Leg2.filterNames = cms.vstring("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle23Ele12Leg2

process.goodElectronsMeasureHLTEle23Ele12DZ = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle23Ele12DZ.filterNames = cms.vstring("hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle23Ele12DZ

# tags
process.tagElectronsMeasureHLTEle17Ele12Leg1 = process.goodElectronsMeasureHLTEle17Ele12Leg1.clone()
process.tagElectronsMeasureHLTEle17Ele12Leg1.inputs = cms.InputTag("goodElectronsTagHLT")
process.goodElectronsMeasureHLT += process.tagElectronsMeasureHLTEle17Ele12Leg1

process.tagElectronsMeasureHLTEle23Ele12Leg1 = process.goodElectronsMeasureHLTEle23Ele12Leg1.clone()
process.tagElectronsMeasureHLTEle23Ele12Leg1.inputs = cms.InputTag("goodElectronsTagHLT")
process.goodElectronsMeasureHLT += process.tagElectronsMeasureHLTEle23Ele12Leg1



## electron muon
#process.goodElectronsMeasureHLTMu17Ele12ELeg = process.goodElectronsMeasureHLTEle23.clone()
#process.goodElectronsMeasureHLTMu17Ele12ELeg.filterNames = cms.vstring("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter")
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTMu17Ele12ELeg
#
#process.goodElectronsMeasureHLTMu8Ele17ELeg = process.goodElectronsMeasureHLTEle23.clone()
#process.goodElectronsMeasureHLTMu8Ele17ELeg.filterNames = cms.vstring("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter")
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTMu8Ele17ELeg
#
#process.goodElectronsMeasureHLTMu8Ele23ELeg = process.goodElectronsMeasureHLTEle23.clone()
#process.goodElectronsMeasureHLTMu8Ele23ELeg.filterNames = cms.vstring("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter")
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTMu8Ele23ELeg
#
#process.goodElectronsMeasureHLTMu23Ele8ELeg = process.goodElectronsMeasureHLTEle23.clone()
#process.goodElectronsMeasureHLTMu23Ele8ELeg.filterNames = cms.vstring("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter")
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTMu23Ele8ELeg
#
#process.goodElectronsMeasureHLTMu23Ele12ELeg = process.goodElectronsMeasureHLTEle23.clone()
#process.goodElectronsMeasureHLTMu23Ele12ELeg.filterNames = cms.vstring("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter")
#process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTMu23Ele12ELeg

# electron tau
process.goodElectronsMeasureHLTEle22Tau20LegSingleL1 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle22Tau20LegSingleL1.filterNames = cms.vstring("hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle22Tau20LegSingleL1

process.goodElectronsMeasureHLTEle24Tau20LegSingleL1 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle24Tau20LegSingleL1.filterNames = cms.vstring("hltEle24WPLooseL1SingleIsoEG22erGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle24Tau20LegSingleL1

process.goodElectronsMeasureHLTEle27Tau20LegSingleL1 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle27Tau20LegSingleL1.filterNames = cms.vstring("hltEle27erWPLooseGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle27Tau20LegSingleL1

process.goodElectronsMeasureHLTEle32Tau20LegSingleL1 = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle32Tau20LegSingleL1.filterNames = cms.vstring("hltEle32erWPLooseGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle32Tau20LegSingleL1

process.goodElectronsMeasureHLTEle24Tau20Leg = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTEle24Tau20Leg.filterNames = cms.vstring("hltEle24WPLooseL1IsoEG22erTau20erGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle24Tau20Leg

# OR

#soup
# ele25eta2p1tight ele27tight ele27eta2p1loose ele45loose
process.goodElectronsMeasureHLTSingleEleSoup = process.goodElectronsMeasureHLTEle23.clone()
process.goodElectronsMeasureHLTSingleEleSoup.isAND = cms.bool(False)
process.goodElectronsMeasureHLTSingleEleSoup.filterNames = cms.vstring("hltEle25erWPTightGsfTrackIsoFilter","hltEle27WPTightGsfTrackIsoFilter","hltEle27erWPLooseGsfTrackIsoFilter","hltEle45WPLooseGsfTrackIsoFilter")
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTSingleEleSoup



#################
### Tag flags ###
#################
TagFlagsToStore = cms.PSet(
    passingEle17Ele12Leg1 = cms.InputTag("tagElectronsMeasureHLTEle17Ele12Leg1"),
    passingEle23Ele12Leg1 = cms.InputTag("tagElectronsMeasureHLTEle23Ele12Leg1"),
)



###################################################################
## SEQUENCES
###################################################################

process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag(options['ELECTRON_COLL'])
process.ele_sequence = cms.Sequence(
    process.goodElectrons +
    process.egmGsfElectronIDSequence +
    process.goodElectronsPROBECutBasedVeto +
    process.goodElectronsPROBECutBasedLoose +
    process.goodElectronsPROBECutBasedMedium +
    process.goodElectronsPROBECutBasedTight +
    process.goodElectronsTAGCutBasedVeto +
    process.goodElectronsTAGCutBasedLoose +
    process.goodElectronsTAGCutBasedMedium +
    process.goodElectronsTAGCutBasedTight +
    process.goodElectronsTagHLT +
    process.goodElectronsProbeHLT +
    process.goodElectronsProbeMeasureHLT +
    process.goodElectronsMeasureHLT
    )

process.sc_sequence = cms.Sequence(process.superClusterCands +
                                   process.goodSuperClusters +
                                   process.goodSuperClustersHLT +
                                   process.GsfMatchedSuperClusterCands
                                   )

###################################################################
## TnP PAIRS
###################################################################

process.allTagsAndProbes = cms.Sequence()

if (options['DOTRIGGER']):
    process.allTagsAndProbes *= process.tagTightHLT

if (options['DORECO']):
    process.allTagsAndProbes *= process.tagTightSC

if (options['DOID']):
    process.allTagsAndProbes *= process.tagTightRECO

process.mc_sequence = cms.Sequence()

#if (varOptions.isMC):
#    process.mc_sequence *= (process.McMatchHLT + process.McMatchTag + process.McMatchSC + process.McMatchRECO)

##########################################################################
## TREE MAKER OPTIONS
##########################################################################
if (not varOptions.isMC):
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(False)
        )

# remove sc stuff that is broken
if hasattr(CommonStuffForGsfElectronProbe.variables,'probe_sc_energy'): del CommonStuffForGsfElectronProbe.variables.probe_sc_energy
if hasattr(CommonStuffForGsfElectronProbe.variables,'probe_sc_et'):     del CommonStuffForGsfElectronProbe.variables.probe_sc_et
if hasattr(CommonStuffForGsfElectronProbe.variables,'probe_sc_eta'):    del CommonStuffForGsfElectronProbe.variables.probe_sc_eta
if hasattr(CommonStuffForGsfElectronProbe.variables,'probe_sc_abseta'): del CommonStuffForGsfElectronProbe.variables.probe_sc_abseta
if hasattr(CommonStuffForGsfElectronProbe.tagVariables,'sc_energy'):    del CommonStuffForGsfElectronProbe.tagVariables.sc_energy
if hasattr(CommonStuffForGsfElectronProbe.tagVariables,'sc_et'):        del CommonStuffForGsfElectronProbe.tagVariables.sc_et
if hasattr(CommonStuffForGsfElectronProbe.tagVariables,'sc_eta'):       del CommonStuffForGsfElectronProbe.tagVariables.sc_eta
if hasattr(CommonStuffForGsfElectronProbe.tagVariables,'sc_abseta'):    del CommonStuffForGsfElectronProbe.tagVariables.sc_abseta
if hasattr(CommonStuffForSuperClusterProbe.tagVariables,'sc_energy'):   del CommonStuffForSuperClusterProbe.tagVariables.sc_energy
if hasattr(CommonStuffForSuperClusterProbe.tagVariables,'sc_et'):       del CommonStuffForSuperClusterProbe.tagVariables.sc_et
if hasattr(CommonStuffForSuperClusterProbe.tagVariables,'sc_eta'):      del CommonStuffForSuperClusterProbe.tagVariables.sc_eta
if hasattr(CommonStuffForSuperClusterProbe.tagVariables,'sc_abseta'):   del CommonStuffForSuperClusterProbe.tagVariables.sc_abseta

process.GsfElectronToTrigger = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                              CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
                                              tagProbePairs = cms.InputTag("tagTightHLT"),
                                              arbitration   = cms.string("Random2"),
                                              flags         = cms.PSet(
                                                  #passingHLTEle22Eta2p1            = cms.InputTag("goodElectronsMeasureHLTEle22Eta2p1"),
                                                  #passingHLTEle24Eta2p1            = cms.InputTag("goodElectronsMeasureHLTEle24Eta2p1"),
                                                  #passingHLTEle25Eta2p1            = cms.InputTag("goodElectronsMeasureHLTEle25Eta2p1"),
                                                  passingHLTEle27Eta2p1            = cms.InputTag("goodElectronsMeasureHLTEle27Eta2p1"),
                                                  passingHLTEle25Eta2p1Tight       = cms.InputTag("goodElectronsMeasureHLTEle25Eta2p1Tight"),
                                                  passingHLTEle27Eta2p1Tight       = cms.InputTag("goodElectronsMeasureHLTEle27Eta2p1Tight"),
                                                  passingHLTEle32Eta2p1Tight       = cms.InputTag("goodElectronsMeasureHLTEle32Eta2p1Tight"),
                                                  #passingHLTEle23                  = cms.InputTag("goodElectronsMeasureHLTEle23"),
                                                  #passingHLTEle27                  = cms.InputTag("goodElectronsMeasureHLTEle27"),
                                                  #passingHLTEle35                  = cms.InputTag("goodElectronsMeasureHLTEle35"),
                                                  passingHLTEle45                  = cms.InputTag("goodElectronsMeasureHLTEle45"),
                                                  #passingHLTEle25Tight             = cms.InputTag("goodElectronsMeasureHLTEle25Tight"),
                                                  passingHLTEle27Tight             = cms.InputTag("goodElectronsMeasureHLTEle27Tight"),
                                                  passingHLTEle17                  = cms.InputTag("goodElectronsMeasureHLTEle17"),
                                                  passingHLTEle12                  = cms.InputTag("goodElectronsMeasureHLTEle12"),
                                                  passingHLTEle17Ele12Leg1         = cms.InputTag("goodElectronsMeasureHLTEle17Ele12Leg1"),
                                                  passingHLTEle17Ele12Leg2         = cms.InputTag("goodElectronsMeasureHLTEle17Ele12Leg2"),
                                                  passingHLTEle17Ele12DZ           = cms.InputTag("goodElectronsMeasureHLTEle17Ele12DZ"),
                                                  passingHLTEle23Ele12Leg1         = cms.InputTag("goodElectronsMeasureHLTEle23Ele12Leg1"),
                                                  passingHLTEle23Ele12Leg2         = cms.InputTag("goodElectronsMeasureHLTEle23Ele12Leg2"),
                                                  passingHLTEle23Ele12DZ           = cms.InputTag("goodElectronsMeasureHLTEle23Ele12DZ"),
                                                  #passingHLTMu8Ele17ELeg           = cms.InputTag("goodElectronsMeasureHLTMu8Ele17ELeg"),
                                                  #passingHLTMu8Ele23ELeg           = cms.InputTag("goodElectronsMeasureHLTMu8Ele23ELeg"),
                                                  #passingHLTMu17Ele12ELeg          = cms.InputTag("goodElectronsMeasureHLTMu17Ele12ELeg"),
                                                  #passingHLTMu23Ele8ELeg           = cms.InputTag("goodElectronsMeasureHLTMu23Ele8ELeg"),
                                                  #passingHLTMu23Ele12ELeg          = cms.InputTag("goodElectronsMeasureHLTMu23Ele12ELeg"),
                                                  passingHLTEle22Tau20LegSingleL1  = cms.InputTag("goodElectronsMeasureHLTEle22Tau20LegSingleL1"),
                                                  passingHLTEle24Tau20LegSingleL1  = cms.InputTag("goodElectronsMeasureHLTEle24Tau20LegSingleL1"),
                                                  passingHLTEle27Tau20LegSingleL1  = cms.InputTag("goodElectronsMeasureHLTEle27Tau20LegSingleL1"),
                                                  passingHLTEle32Tau20LegSingleL1  = cms.InputTag("goodElectronsMeasureHLTEle32Tau20LegSingleL1"),
                                                  passingHLTEle24Tau20Leg          = cms.InputTag("goodElectronsMeasureHLTEle24Tau20Leg"),
                                                  passingHLTSingleEleSoup          = cms.InputTag("goodElectronsMeasureHLTSingleEleSoup"),
                                                                       ),                                               
                                              allProbes     = cms.InputTag("goodElectronsProbeMeasureHLT"),
                                              )

process.GsfElectronToTrigger.tagFlags = cms.PSet(TagFlagsToStore)

if (varOptions.isMC):
    #process.GsfElectronToTrigger.probeMatches  = cms.InputTag("McMatchHLT")
    process.GsfElectronToTrigger.eventWeight   = cms.InputTag("generator")
    process.GsfElectronToTrigger.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")

process.GsfElectronToSC = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                         CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
                                         tagProbePairs = cms.InputTag("tagTightSC"),
                                         arbitration   = cms.string("Random2"),
                                         flags         = cms.PSet(passingRECO   = cms.InputTag("GsfMatchedSuperClusterCands", "superclusters"),         
                                                                  ),                                               
                                         allProbes     = cms.InputTag("goodSuperClustersHLT"),
                                         )

if (varOptions.isMC):
    #process.GsfElectronToSC.probeMatches  = cms.InputTag("McMatchSC")
    process.GsfElectronToSC.eventWeight   = cms.InputTag("generator")
    process.GsfElectronToSC.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")

process.GsfElectronToRECO = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                           mcTruthCommonStuff, CommonStuffForGsfElectronProbe,
                                           tagProbePairs = cms.InputTag("tagTightRECO"),
                                           arbitration   = cms.string("Random2"),
                                           flags         = cms.PSet(passingVeto   = cms.InputTag("goodElectronsPROBECutBasedVeto"),
                                                                    passingLoose  = cms.InputTag("goodElectronsPROBECutBasedLoose"),
                                                                    passingMedium = cms.InputTag("goodElectronsPROBECutBasedMedium"),
                                                                    passingTight  = cms.InputTag("goodElectronsPROBECutBasedTight"),
                                                                    ),                                               
                                           allProbes     = cms.InputTag("goodElectronsProbeHLT"),
                                           )

if (varOptions.isMC):
    #process.GsfElectronToRECO.probeMatches  = cms.InputTag("McMatchRECO")
    process.GsfElectronToRECO.eventWeight   = cms.InputTag("generator")
    process.GsfElectronToRECO.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")

process.tree_sequence = cms.Sequence()
if (options['DOTRIGGER']):
    process.tree_sequence *= process.GsfElectronToTrigger

if (options['DORECO']):
    process.tree_sequence *= process.GsfElectronToSC

if (options['DOID']):
    process.tree_sequence *= process.GsfElectronToRECO

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
        process.sc_sequence +
        process.allTagsAndProbes +
        process.pileupReweightingProducer +
        process.mc_sequence +
        process.eleVarHelper +
        process.tree_sequence
        )
else:
    process.p = cms.Path(
        process.sampleInfo +
        process.hltFilter +
        process.ele_sequence + 
        process.sc_sequence +
        process.allTagsAndProbes +
        process.mc_sequence +
        process.eleVarHelper +
        process.tree_sequence
        )

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string(options['OUTPUT_FILE_NAME']),
    closeFileFast = cms.untracked.bool(True)
    )
