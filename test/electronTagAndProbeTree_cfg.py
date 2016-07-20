import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
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
    options['OUTPUT_FILE_NAME']    = "TnPTree_mc_electron.root"
    options['TnPPATHS']            = cms.vstring("HLT*")
    options['TnPHLTTagFilters']    = cms.vstring()
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("")
    options['GLOBALTAG']           = 'auto:run2_mc'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()
else:
    options['OUTPUT_FILE_NAME']    = "TnPTree_data_electron.root"
    options['TnPPATHS']            = cms.vstring("HLT_Ele27_eta2p1_WPLoose_Gsf_v*")
    options['TnPHLTTagFilters']    = cms.vstring("hltEle27erWPLooseGsfTrackIsoFilter")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("hltEle27erWPLooseGsfTrackIsoFilter")
    options['GLOBALTAG']           = 'auto:run2_data'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()

###################################################################
## Inputs for test
###################################################################
filesMC =  cms.untracked.vstring(
    '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/000FF6AC-9F2A-E611-A063-0CC47A4C8EB0.root',
    )

filesData =  cms.untracked.vstring( 
    '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/150/00000/0A6284C7-D719-E611-93E6-02163E01421D.root',
    '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/06277EC1-181A-E611-870F-02163E0145E5.root',
    '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/0A7BD549-131A-E611-8287-02163E0134FC.root',
    )


if options['useAOD']:
    filesMC = cms.untracked.vstring(
        '/store/mc/RunIISpring16DR80/DYToEE_NNPDF30_13TeV-powheg-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/000AB373-AD0C-E611-8261-0002C94CD120.root',
        '/store/mc/RunIISpring16DR80/DYToEE_NNPDF30_13TeV-powheg-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/0020B9BE-D20C-E611-AD17-A0369F310374.root',
        '/store/mc/RunIISpring16DR80/DYToEE_NNPDF30_13TeV-powheg-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/004A2D66-D30C-E611-A844-0CC47A009258.root',
        )
    filesData = cms.untracked.vstring(
       '/store/data/Run2016B/SingleElectron/AOD/PromptReco-v2/000/273/158/00000/100B4FC6-1D1A-E611-AADB-02163E0118D5.root',
       '/store/data/Run2016B/SingleElectron/AOD/PromptReco-v2/000/273/158/00000/1032C8FE-211A-E611-BE71-02163E01387F.root',
       '/store/data/Run2016B/SingleElectron/AOD/PromptReco-v2/000/273/158/00000/12E3D1D6-0A1A-E611-85A9-02163E011BF0.root',
       '/store/data/Run2016B/SingleElectron/AOD/PromptReco-v2/000/273/158/00000/18C2A248-101A-E611-9906-02163E01414F.root',
       '/store/data/Run2016B/SingleElectron/AOD/PromptReco-v2/000/273/158/00000/18CB30C7-181A-E611-BA7E-02163E014573.root',
        )

options['INPUT_FILE_NAME'] = filesData
if varOptions.isMC:
    options['INPUT_FILE_NAME'] = filesMC

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

process.goodElectronsMeasureHLTEle23 = cms.EDProducer("PatElectronTriggerCandProducer",
                                                filterNames = cms.vstring("hltEle23WPLooseGsfTrackIsoFilter"),
                                                #inputs      = cms.InputTag("goodElectronsProbeMeasureHLT"),
                                                inputs      = cms.InputTag("goodElectrons"),
                                                bits        = cms.InputTag('TriggerResults::HLT'),
                                                objects     = cms.InputTag('selectedPatTrigger'),
                                                dR          = cms.double(0.3),
                                                isAND       = cms.bool(False)
                                                )
process.goodElectronsMeasureHLT += process.goodElectronsMeasureHLTEle23

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

process.ele_sequence += process.goodElectronsMeasureHLT


#################
### Tag flags ###
#################
TagFlagsToStore = cms.PSet(
    passingEle17Ele12Leg1 = cms.InputTag("tagElectronsMeasureHLTEle17Ele12Leg1"),
    passingEle23Ele12Leg1 = cms.InputTag("tagElectronsMeasureHLTEle23Ele12Leg1"),
)



process.mc_sequence = cms.Sequence()


##########################################################################
## TREE MAKER OPTIONS
##########################################################################    
process.GsfElectronToTrigger = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                              tnpVars.CommonStuffForGsfElectronProbe, tnpVars.mcTruthCommonStuff,
                                              tagProbePairs = cms.InputTag("tagTightHLT"),
                                              arbitration   = cms.string("HighestPt"),
                                              flags         = cms.PSet(#passingHLT                       = cms.InputTag("goodElectronsMeasureHLT"),
                                                                       passingLoose                     = cms.InputTag("goodElectronsPROBECutBasedLoose"),
                                                                       passingMedium                    = cms.InputTag("goodElectronsPROBECutBasedMedium"),
                                                                       passingTight                     = cms.InputTag("goodElectronsPROBECutBasedTight"),
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
                                              #allProbes     = cms.InputTag("goodElectronsProbeMeasureHLT"),
                                              allProbes     = cms.InputTag("goodElectrons"),
                                              )
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

# manually fix pileup
from SimGeneral.MixingModule.mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU_cfi import mix
from DevTools.TagAndProbe.pileup_cfi import currentPileup
pu_distribs = { "mc" : mix.input.nbPileupEvents.probValue }
data_pu_distribs = { "golden" : currentPileup }

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
