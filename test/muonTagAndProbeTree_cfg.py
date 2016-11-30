import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from collections import OrderedDict
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

from PhysicsTools.TagAndProbe.treeMakerOptions_cfi import *

if (varOptions.isMC):
    # TODO update to new MC
    #options['INPUT_FILE_NAME']     = '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/00000/00F0B3DC-211B-E611-A6A0-001E67248A39.root'
    options['INPUT_FILE_NAME']     = '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/000FF6AC-9F2A-E611-A063-0CC47A4C8EB0.root'
    options['OUTPUT_FILE_NAME']    = "TnPTree_mc_muon.root"
    options['TnPPATHS']            = cms.vstring()#"HLT_IsoTkMu24_v*")
    options['TnPHLTTagFilters']    = cms.vstring()#"hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring()
    options['GLOBALTAG']           = 'auto:run2_mc'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()
else:
    options['INPUT_FILE_NAME']     = "/store/data/Run2016H/SingleMuon/MINIAOD/PromptReco-v3/000/284/036/00000/0E02D50E-989F-E611-A962-FA163EE15C80.root"
    options['OUTPUT_FILE_NAME']    = "TnPTree_data_muon.root"
    options['TnPPATHS']            = ["HLT_IsoTkMu24_v*",]
    options['TnPHLTTagFilters']    = ["hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09"]
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['HLTFILTERTOMEASURE']  = cms.vstring("")
    options['GLOBALTAG']           = 'auto:run2_data'
    options['EVENTSToPROCESS']     = cms.untracked.VEventRange()

###################################################################

###########################
### Trigger definitions ###
###########################
# TODO check all menus
menu = 'v4.2'
trigger_filters = OrderedDict()

# HLT_Mu17_TrkIsoVVL_v*
# HLT_Mu8_TrkIsoVVL_v*
trigger_filters['probeTriggersMu17'] = {
    'filterNames': {
        'v4.2': ["hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4"],
    },
    'isAND': True,
    'inputs': 'probeMuons',
 }
trigger_filters['probeTriggersMu8'] = {
    'filterNames': {
        'v4.2': ["hltL3fL1sMu5L1f0L2f5L3Filtered8TkIsoFiltered0p4"],
    },
    'isAND': True,
    'inputs': 'probeMuons',
}
# HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*
# HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*
trigger_filters['probeTriggersMu17Leg'] = {
    'filterNames': {
        'v4.2': ["hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4","hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17"],
    },
    'isAND': True,
    'inputs': 'probeMuons'
}
trigger_filters['probeTriggersMu8Leg'] = {
    'filterNames': {
        'v4.2': ["hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4","hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8"],
    },
    'isAND': True,
    'inputs': 'probeMuons',
}
trigger_filters['probeTriggersTkMu8Leg'] = {
    'filterNames': {
        'v4.2': ["hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4","hltDiMuonGlbFiltered17TrkFiltered8"],
    },
    'isAND': True,  
    'inputs': 'probeMuons',
}
trigger_filters['probeTriggersMu8ORTkMu8LegPre'] = {
    'filterNames': {
        'v4.2': ["hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8","hltDiMuonGlbFiltered17TrkFiltered8"],
    },
    'isAND': False,
    'inputs': 'probeMuons',
}
trigger_filters['probeTriggersMu8ORTkMu8Leg'] = {
    'filterNames': {
        'v4.2': ["hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4","hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4"],
    },
    'isAND': False, 
    'inputs': 'probeTriggersMu8ORTkMu8LegPre',
}
trigger_filters['probeTriggersMu8LegDZ'] = {
    'filterNames': {
        'v4.2': ["hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2"],
    },
    'isAND': True,
    'inputs': 'probeTriggersMu8Leg',
}
trigger_filters['probeTriggersTkMu8LegDZ'] = {
    'filterNames': {
        'v4.2': ["hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2"],
    },
    'isAND': True,
    'inputs': 'probeTriggersTkMu8Leg',
}
trigger_filters['probeTriggersMu8ORTkMu8LegDZ'] = {
    'filterNames': {
        'v4.2': ["hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2","hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2"],
    },
    'isAND': False,
    'inputs': 'probeTriggersMu8ORTkMu8Leg',
}
# HLT_IsoMu24_v*
# HLT_IsoTkMu24_v*
# HLT_Mu50_v*
trigger_filters['probeTriggersIsoMu24'] = {
    'filterNames': {
        'v4.2': ["hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09"],
    },
    'isAND': True,
    'inputs': 'probeMuons',
}
trigger_filters['probeTriggersIsoTkMu24'] = {
    'filterNames': {
        'v4.2': ["hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09"],
    },
    'isAND': True,
    'inputs': 'probeMuons',
}
trigger_filters['probeTriggersMu50'] = {
    'filterNames': {
        'v4.2': ["hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"],
    },
    'isAND': True,
    'inputs': 'probeMuons',
}
trigger_filters['probeTriggersIsoMu24ORIsoTkMu24'] = {
    'filterNames': {
        'v4.2': ["hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09","hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"],
    },
    'isAND': False,
    'inputs': 'probeMuons',
}
trigger_filters['probeTriggersIsoMu24ORIsoTkMu24OrMu50'] = {
    'filterNames': {
        'v4.2': ["hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09","hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09","hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"],
    },
    'isAND': False,
    'inputs': 'probeMuons',
}



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


######################
### Trigger Probes ###
######################

process.probeTriggerSeq = cms.Sequence()


for trigger,args in trigger_filters.iteritems():
    mod = cms.EDProducer(
        "PatMuonTriggerCandProducer",
        filterNames = cms.vstring(*args['filterNames'][menu]),
        inputs = cms.InputTag(args['inputs']),
        bits        = cms.InputTag('TriggerResults::HLT'),
        objects     = cms.InputTag('selectedPatTrigger'),
        dR          = cms.double(0.1),
        isAND       = cms.bool(args['isAND'])
    )
    setattr(process, trigger, mod)
    process.probeTriggerSeq += getattr(process,trigger)

# tag trigger
process.tagTriggersMu17Leg = cms.EDProducer("PatMuonTriggerCandProducer",
    filterNames = cms.vstring(*trigger_filters['probeTriggersMu17Leg']['filterNames'][menu]),
    inputs      = cms.InputTag("tagMuonsTriggerMatched"),
    bits        = cms.InputTag('TriggerResults::HLT'),
    objects     = cms.InputTag('selectedPatTrigger'),
    dR          = cms.double(0.1),
    isAND       = cms.bool(True)
    )
process.probeTriggerSeq += process.tagTriggersMu17Leg


###################################################################
## TnP PAIRS
###################################################################

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuonsTriggerMatched@+ probeMuons@-"), # charge coniugate states are implied
    cut   = cms.string("40 < mass < 200")
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
from SimGeneral.MixingModule.mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi import mix
from DevTools.TagAndProbe.pileup_cfi import pileup2016
pu_distribs = { "mc": mix.input.nbPileupEvents.probValue }
data_pu_distribs = { "golden": pileup2016 }

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


TagFlagsToStore = cms.PSet(
    passingMu17                 = cms.InputTag("tagTriggersMu17Leg"),
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
        passingHppLooseID                     = cms.string("userInt('isMediumMuonICHEP')==1 && "+trackIsoDef+" < 0.4"), 
        passingHppLooseIso                    = cms.string(isolationDef+" < 0.25"),
        passingHppMediumID                    = cms.string("userInt('isMediumMuonICHEP')==1 && "+trackIsoDef+" < 0.4 && userFloat('dz') < 0.1 && ((pt<20 && userFloat('dB2D') < 0.01) || (pt>=20 && userFloat('dB2D') < 0.02))"), 
        passingHppMediumIso                   = cms.string(isolationDef+" < 0.15"),
    ),
    allProbes     = cms.InputTag("probeMuons"),
    )

for trigger in trigger_filters:
    setattr(process.muonEffs.flags,trigger.replace('probeTriggers','passing'),cms.InputTag(trigger))

process.tpPairSeq = cms.Sequence(
    process.tpPairs
)

if varOptions.isMC :
    process.tpPairSeq += process.muMcMatch
    #process.tpPairSeq += process.tpPairsMCEmbedded
    process.tpPairSeq += process.pileupReweightingProducer
    process.muonEffs.isMC = cms.bool(True)
    process.muonEffs.eventWeight   = cms.InputTag("generator")
    process.muonEffs.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")
    #setattr(process.muonEffs.pairVariables, 'mc_mass', cms.string("userFloat('mc_mass')"))

process.p = cms.Path(
    process.sampleInfo *
    process.idEmbedSequence *
    (process.tagMuons + process.probeMuons) *
    (process.tagMuonsTriggerMatched + process.probeTriggerSeq) *
    process.tpPairSeq *
    process.muonEffs
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

