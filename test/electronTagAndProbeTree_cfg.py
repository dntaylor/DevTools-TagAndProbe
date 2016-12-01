import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from collections import OrderedDict
import sys

process = cms.Process("tnp")

###################################################################
## argument line options
###################################################################
varOptions = VarParsing('analysis')
varOptions.register(
    "isMC", True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Compute MC efficiencies"
    )

varOptions.parseArguments()


###################################################################
## Define TnP inputs 
###################################################################

options = dict()
options['useAOD']               = cms.bool(False)

options['HLTProcessName']       = "HLT"

### set input collections
options['ELECTRON_COLL']        = "slimmedElectrons"
options['PHOTON_COLL']          = "slimmedPhotons"
options['SUPERCLUSTER_COLL']    = "reducedEgamma:reducedSuperClusters" ### not used in AOD
if options['useAOD']:
    options['ELECTRON_COLL']        = "gedGsfElectrons"


options['ELECTRON_CUTS']        = "ecalEnergy*sin(superClusterPosition.theta)>5.0 &&  (abs(-log(tan(superClusterPosition.theta/2)))<2.5)"
options['SUPERCLUSTER_CUTS']    = "abs(eta)<2.5 &&  et>5.0"
options['PHOTON_CUTS']          = "(abs(-log(tan(superCluster.position.theta/2)))<=2.5) && pt> 10"

options['ELECTRON_TAG_CUTS']    = "(abs(-log(tan(superCluster.position.theta/2)))<=2.1) && !(1.4442<=abs(-log(tan(superClusterPosition.theta/2)))<=1.566) && pt >= 30.0"

options['MAXEVENTS']            = cms.untracked.int32(varOptions.maxEvents) 
options['DoTrigger']            = cms.bool(True)
options['DoRECO']               = cms.bool(True)
options['DoEleID']              = cms.bool(True)
options['DoPhoID']              = cms.bool(False)

options['OUTPUTEDMFILENAME']    = 'edmFile.root'
options['DEBUG']                = cms.bool(False)


if (varOptions.isMC):
    # TODO add MC triggers
    options['OUTPUT_FILE_NAME']    = "TnPTree_mc_electron.root"
    options['TnPPATHS']            = cms.vstring()
    options['TnPHLTTagFilters']    = cms.vstring()
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("")
    options['GLOBALTAG']           = 'auto:run2_mc'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()
else:
    options['OUTPUT_FILE_NAME']    = "TnPTree_data_electron.root"
    options['TnPPATHS']            = cms.vstring("HLT_Ele27_eta2p1_WPTight_Gsf_v*")
    options['TnPHLTTagFilters']    = cms.vstring("hltEle27erWPTightGsfTrackIsoFilter")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("hltEle27erWPTightGsfTrackIsoFilter")
    options['GLOBALTAG']           = 'auto:run2_data'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()

###################################################################
## Inputs for test
###################################################################
# TODO update MC
filesMC =  cms.untracked.vstring(
    '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/000FF6AC-9F2A-E611-A063-0CC47A4C8EB0.root',
    )

filesData =  cms.untracked.vstring( 
    '/store/data/Run2016H/SingleElectron/MINIAOD/PromptReco-v2/000/283/353/00000/A6165E3D-F696-E611-826F-02163E013542.root',
    )


options['INPUT_FILE_NAME'] = filesData
if varOptions.isMC:
    options['INPUT_FILE_NAME'] = filesMC

###########################
### Trigger definitions ###
###########################
# TODO check all menus
menu = 'v4.2'
trigger_filters = OrderedDict()

#HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*
#HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*
#HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v*
trigger_filters['probeTriggersEle23'] = {
    'filterNames': {
        'v4.2': ["hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter"],
    },
    'isAND': True,
    'inputs': 'goodElectrons',
}
trigger_filters['probeTriggersEle17'] = {
    'filterNames': {
        'v4.2': ["hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter"],
    },
    'isAND': True,
    'inputs': 'goodElectrons',
}
trigger_filters['probeTriggersEle12'] = {
    'filterNames': {
        'v4.2': ["hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter"],
    },
    'isAND': True,
    'inputs': 'goodElectrons',
}

# HLT_Ele27_WPTight_Gsf_v*
# HLT_Ele27_eta2p1_WPTight_Gsf_v*
trigger_filters['probeTriggersEle27WPTight'] = {
    'filterNames': {
        'v4.2': ["hltEle27WPTightGsfTrackIsoFilter"],
    },
    'isAND': True,
    'inputs': 'goodElectrons',
}
trigger_filters['probeTriggersEle27Eta2p1WPTight'] = {
    'filterNames': {
        'v4.2': ["hltEle27erWPTightGsfTrackIsoFilter"],
    },
    'isAND': True,
    'inputs': 'goodElectrons',
}

# HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*
trigger_filters['probeTriggersEle23Leg'] = {
    'filterNames': {
        'v4.2': ["hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter"],
    },
    'isAND': True,
    'inputs': 'goodElectrons',
}
trigger_filters['probeTriggersEle12Leg'] = {
    'filterNames': {
        'v4.2': ["hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter"],
    },
    'isAND': True,
    'inputs': 'goodElectrons',
}
trigger_filters['probeTriggersEle12LegDZ'] = {
    'filterNames': {
        'v4.2': ["hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter"],
    },
    'isAND': True,
    'inputs': 'probeTriggersEle12Leg',
}


###################################################################
## import TnP tree maker pythons and configure for AODs
###################################################################
if options['useAOD']:
    import PhysicsTools.TagAndProbe.treeMakerOptionsAOD_cfi as tnpTreeMaker
else: 
    import PhysicsTools.TagAndProbe.treeMakerOptions_cfi as tnpTreeMaker

tnpTreeMaker.setModules(process,options)

import PhysicsTools.TagAndProbe.treeContent_cfi as tnpVars
if options['useAOD']:
    tnpVars.setupTnPVariablesForAOD()

if not varOptions.isMC:
    tnpVars.mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(False)
        )


###################################################################
## Init and Load
###################################################################
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options['GLOBALTAG'] , '')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            fileNames = options['INPUT_FILE_NAME'],
                            eventsToProcess = options['EVENTSToPROCESS']
                            )

process.maxEvents = cms.untracked.PSet( input = options['MAXEVENTS'])

###################################################################
## ID
###################################################################
import PhysicsTools.TagAndProbe.electronIDModules_cfi as egmEleID
import PhysicsTools.TagAndProbe.photonIDModules_cfi   as egmPhoID
egmEleID.setIDs(process, options)
egmPhoID.setIDs(process, options)
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag(options['ELECTRON_COLL'])
process.egmPhotonIDs.physicsObjectSrc      = cms.InputTag(options['PHOTON_COLL'])


###################################################################
## SEQUENCES
###################################################################
tnpTreeMaker.setSequences(process,options)
process.cand_sequence = cms.Sequence( process.tag_sequence )

if (options['DoEleID']):
    process.cand_sequence += process.ele_sequence
    print "  -- Producing electron SF tree    -- "

if (options['DoPhoID']):
    process.cand_sequence += process.pho_sequence
    print "  -- Producing photon SF tree      -- "

if (options['DoTrigger']):
    process.cand_sequence += process.ele_sequence
    process.cand_sequence += process.hlt_sequence
    print "  -- Producing HLT efficiency tree -- "

if (options['DoRECO']):
    process.cand_sequence += process.sc_sequence
    print "  -- Producing RECO SF tree        -- "



###################################################################
## TnP PAIRS
###################################################################
process.allTagsAndProbes = cms.Sequence()

if (options['DoTrigger']):
    process.allTagsAndProbes *= process.tagTightHLT

if (options['DoRECO']):
    process.allTagsAndProbes *= process.tagTightSC

if (options['DoEleID']):
    process.allTagsAndProbes *= process.tagTightEleID

if (options['DoPhoID']):
    process.allTagsAndProbes *= process.tagTightPhoID


####################
### Add triggers ###
####################
process.goodElectronsMeasureHLT = cms.Sequence()


for trigger,args in trigger_filters.iteritems():
    mod = cms.EDProducer(
        "PatElectronTriggerCandProducer",
        filterNames = cms.vstring(*args['filterNames'][menu]),
        inputs = cms.InputTag(args['inputs']),
        bits        = cms.InputTag('TriggerResults::HLT'),
        objects     = cms.InputTag('selectedPatTrigger'),
        dR          = cms.double(0.1),
        isAND       = cms.bool(args['isAND'])
    )
    setattr(process, trigger, mod)
    process.goodElectronsMeasureHLT += getattr(process,trigger)

# tag trigger
process.tagTriggersEle23Leg = cms.EDProducer("PatElectronTriggerCandProducer",
    filterNames = cms.vstring(*trigger_filters['probeTriggersEle23Leg']['filterNames'][menu]),
    inputs      = cms.InputTag("goodElectronsTagHLT"),
    bits        = cms.InputTag('TriggerResults::HLT'),
    objects     = cms.InputTag('selectedPatTrigger'),
    dR          = cms.double(0.1),
    isAND       = cms.bool(True)
    )
process.goodElectronsMeasureHLT += process.tagTriggersEle23Leg


#################
### Tag flags ###
#################
TagFlagsToStore = cms.PSet(
    passingEle23Ele12Leg1 = cms.InputTag("tagTriggersEle23Leg"),
)



process.mc_sequence = cms.Sequence()


##########################################################################
## TREE MAKER OPTIONS
##########################################################################    
process.GsfElectronToTrigger = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                              tnpVars.CommonStuffForGsfElectronProbe, tnpVars.mcTruthCommonStuff,
                                              tagProbePairs = cms.InputTag("tagTightHLT"),
                                              arbitration   = cms.string("HighestPt"),
                                              flags         = cms.PSet(
                                                                       passingLoose                     = cms.InputTag("goodElectronsPROBECutBasedLoose"),
                                                                       passingMedium                    = cms.InputTag("goodElectronsPROBECutBasedMedium"),
                                                                       passingTight                     = cms.InputTag("goodElectronsPROBECutBasedTight"),
                                                                       ),                                               
                                              allProbes     = cms.InputTag("goodElectrons"),
                                              )
for trigger in trigger_filters:
    setattr(process.GsfElectronToTrigger.flags,trigger.replace('probeTriggers','passing'),cms.InputTag(trigger))

process.GsfElectronToTrigger.tagFlags = cms.PSet(TagFlagsToStore)

process.GsfElectronToSC = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                         tnpVars.mcTruthCommonStuff, tnpVars.CommonStuffForSuperClusterProbe, 
                                         tagProbePairs = cms.InputTag("tagTightSC"),
                                         arbitration   = cms.string("HighestPt"),
                                         flags         = cms.PSet(passingRECO   = cms.InputTag("GsfMatchedSuperClusterCands", "superclusters")
                                                                  ),                                               
                                         allProbes     = cms.InputTag("goodSuperClustersHLT"),
                                         )

process.GsfElectronToEleID = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                            tnpVars.mcTruthCommonStuff, tnpVars.CommonStuffForGsfElectronProbe,
                                            tagProbePairs = cms.InputTag("tagTightEleID"),
                                            arbitration   = cms.string("HighestPt"),
                                            flags         = cms.PSet(passingVeto   = cms.InputTag("goodElectronsPROBECutBasedVeto"),
                                                                     passingLoose  = cms.InputTag("goodElectronsPROBECutBasedLoose"),
                                                                     passingMedium = cms.InputTag("goodElectronsPROBECutBasedMedium"),
                                                                     passingTight  = cms.InputTag("goodElectronsPROBECutBasedTight"),
                                                                     ),                                               
                                            allProbes     = cms.InputTag("goodElectronsProbeHLT"),
                                            )

process.GsfElectronToPhoID = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                            tnpVars.mcTruthCommonStuff, tnpVars.CommonStuffForPhotonProbe,
                                            tagProbePairs = cms.InputTag("tagTightPhoID"),
                                            arbitration   = cms.string("HighestPt"),
                                            flags         = cms.PSet(passingLoose  = cms.InputTag("goodPhotonsPROBECutBasedLoose"),
                                                                     passingMedium = cms.InputTag("goodPhotonsPROBECutBasedMedium"),
                                                                     passingTight  = cms.InputTag("goodPhotonsPROBECutBasedTight"),
                                                                     passingMVA    = cms.InputTag("goodPhotonsPROBEMVA"),
                                                                     ),                                                                                           
                                            allProbes     = cms.InputTag("goodPhotonsProbeHLT"),
                                            )


##############
### Pileup ###
##############
from SimGeneral.MixingModule.mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi import mix
from DevTools.TagAndProbe.pileup_cfi import pileup2016
pu_distribs = { "mc": mix.input.nbPileupEvents.probValue }
data_pu_distribs = { "golden": pileup2016 }

process.pileupReweightingProducer.PileupMC = cms.vdouble(pu_distribs['mc'])
process.pileupReweightingProducer.PileupData = cms.vdouble(data_pu_distribs["golden"])

if (varOptions.isMC):
    process.GsfElectronToTrigger.eventWeight = cms.InputTag("generator")
    process.GsfElectronToEleID.eventWeight   = cms.InputTag("generator")
    process.GsfElectronToPhoID.eventWeight   = cms.InputTag("generator")
    process.GsfElectronToSC.eventWeight      = cms.InputTag("generator")
    process.GsfElectronToTrigger.PUWeightSrc = cms.InputTag("pileupReweightingProducer","pileupWeights")
    process.GsfElectronToEleID.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")
    process.GsfElectronToPhoID.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")
    process.GsfElectronToSC.PUWeightSrc      = cms.InputTag("pileupReweightingProducer","pileupWeights")


process.tree_sequence = cms.Sequence()
if (options['DoTrigger']):
    process.tree_sequence *= process.GsfElectronToTrigger

if (options['DoRECO']):
    process.tree_sequence *= process.GsfElectronToSC

if (options['DoEleID']):
    process.tree_sequence *= process.GsfElectronToEleID

if (options['DoPhoID']):
    process.tree_sequence *= process.GsfElectronToPhoID

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
        process.sampleInfo    +
        process.hltFilter     +
        process.cand_sequence + 
        process.allTagsAndProbes +
        process.pileupReweightingProducer +
        process.mc_sequence +
        process.eleVarHelper +
        process.tree_sequence
        )
else:
    process.p = cms.Path(
        process.sampleInfo    +
        process.hltFilter     +
        process.cand_sequence +
        process.allTagsAndProbes +
        process.mc_sequence  +
        process.eleVarHelper +
        process.tree_sequence
        )

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string(options['OUTPUT_FILE_NAME']),
    closeFileFast = cms.untracked.bool(True)
    )
