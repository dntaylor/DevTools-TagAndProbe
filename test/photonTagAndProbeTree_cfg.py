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
options['DoEleID']              = cms.bool(False)
options['DoPhoID']              = cms.bool(True)

options['OUTPUTEDMFILENAME']    = 'edmFile.root'
options['DEBUG']                = cms.bool(False)

options['UseCalibEn']           = cms.bool(False) # TODO, fix
options['isMC']                 = cms.bool(varOptions.isMC)

if (varOptions.isMC):
    # TODO add MC triggers
    options['OUTPUT_FILE_NAME']    = "TnPTree_mc_photon.root"
    options['TnPPATHS']            = cms.vstring("HLT_Ele27_eta2p1_WPTight_Gsf_v*")
    options['TnPHLTTagFilters']    = cms.vstring("hltEle27erWPTightGsfTrackIsoFilter")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("hltEle27erWPTightGsfTrackIsoFilter")
    options['GLOBALTAG']           = 'auto:run2_mc'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()
else:
    options['OUTPUT_FILE_NAME']    = "TnPTree_data_photon.root"
    options['TnPPATHS']            = cms.vstring("HLT_Ele27_eta2p1_WPTight_Gsf_v*")
    options['TnPHLTTagFilters']    = cms.vstring("hltEle27erWPTightGsfTrackIsoFilter")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("hltEle27erWPTightGsfTrackIsoFilter")
    options['GLOBALTAG']           = 'auto:run2_data'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()

###################################################################
## Inputs for test
###################################################################
filesMC =  cms.untracked.vstring(
    #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/005ED0EB-79F1-E611-B6DA-02163E011C2B.root',
    '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/02A210D6-F5C3-E611-B570-008CFA197BD4.root',
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


# HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*
# 30: hltEG30LRId85ORIso60CaloId15b35eANDHE12R9Id50b80eLegCombLastFilter
# 18: hltEG18R9Id85b90eHE12R9Id50b80eR9UnseededLastFilter or  hltEG18Iso60CaloId15b35eHE12R9Id50b80eTrackIsoUnseededLastFilter
# 90: hltDiEG18R9Id85b90eORIso60CaloId15b35eANDHE12R9Id50b80eMass90CombMassLastFilter

trigger_filters['probeTriggersPho30Leg'] = {
    'filterNames': {
        'v4.2': ["hltEG30LRId85ORIso60CaloId15b35eANDHE12R9Id50b80eLegCombLastFilter"],
    },
    'isAND': True,
    #'inputs': 'goodPhotons',
    'inputs': 'goodElectrons',
}
trigger_filters['probeTriggersPho18Leg'] = {
    'filterNames': {
        'v4.2': ["hltEG18R9Id85b90eHE12R9Id50b80eR9UnseededLastFilter","hltEG18Iso60CaloId15b35eHE12R9Id50b80eTrackIsoUnseededLastFilter"],
    },
    'isAND': False,
    #'inputs': 'goodPhotons',
    'inputs': 'goodElectrons',
}
trigger_filters['probeTriggersPho18LegM90'] = {
    'filterNames': {
        'v4.2': ["hltDiEG18R9Id85b90eORIso60CaloId15b35eANDHE12R9Id50b80eMass90CombMassLastFilter"],
    },
    'isAND': True,
    #'inputs': 'goodPhotons',
    'inputs': 'goodElectrons',
}

# HLT_DoublePhoton60_v*
# 60: hltEG60HEFilter
# 2*60: hltDiEG60HEUnseededFilter

trigger_filters['probeTriggersPho60Leg'] = {
    'filterNames': {
        'v4.2': ["hltEG60HEFilter"],
    },
    'isAND': True,
    #'inputs': 'goodPhotons',
    'inputs': 'goodElectrons',
}

# HLT_Photon175_v*
# 175: hltEG175HEFilter

trigger_filters['probeTriggersPho175'] = {
    'filterNames': {
        'v4.2': ["hltEG175HEFilter"],
    },
    'isAND': True,
    #'inputs': 'goodPhotons',
    'inputs': 'goodElectrons',
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
        dR          = cms.double(0.3),
        isAND       = cms.bool(args['isAND'])
    )
    setattr(process, trigger, mod)
    process.goodElectronsMeasureHLT += getattr(process,trigger)


process.mc_sequence = cms.Sequence()


##########################################################################
## TREE MAKER OPTIONS
##########################################################################    
process.GsfElectronToTrigger = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                              tnpVars.CommonStuffForGsfElectronProbe, tnpVars.mcTruthCommonStuff,
                                              tagProbePairs = cms.InputTag("tagTightHLT"),
                                              arbitration   = cms.string("HighestPt"),
                                              flags         = cms.PSet(
                                                                     #passingLoose  = cms.InputTag("goodPhotonsPROBECutBasedLoose"),
                                                                     #passingMedium = cms.InputTag("goodPhotonsPROBECutBasedMedium"),
                                                                     #passingTight  = cms.InputTag("goodPhotonsPROBECutBasedTight"),
                                                                     #passingMVA    = cms.InputTag("goodPhotonsPROBEMVA"),
                                                                       ),                                               
                                              allProbes     = cms.InputTag("goodElectrons"),
                                              )
for trigger in trigger_filters:
    setattr(process.GsfElectronToTrigger.flags,trigger.replace('probeTriggers','passing'),cms.InputTag(trigger))

#process.GsfElectronToTrigger.tagFlags = cms.PSet(TagFlagsToStore)

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
