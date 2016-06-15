import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

process = cms.Process("tnp")

###################################################################
options = VarParsing('analysis')
options.register(
    "isMC",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Compute MC efficiencies"
    )

options.parseArguments()

config = dict()
isolationDef = "(chargedHadronIso+max(photonIso+neutralHadronIso-0.5*puChargedHadronIso,0.0))/pt"
config['HLTProcessName']          = "HLT"
config['MUON_COLL']               = "slimmedMuons"
config['MUON_CUTS']               = "userInt('tightID')==1 && pt > 5 && abs(eta)<2.4 && "+isolationDef+" < 0.15"
config['MUON_TAG_CUTS']           = "userInt('tightID')==1 && pt > 22 && abs(eta) < 2.4 &&"+isolationDef+"<0.1"
#config['MUON_CUTS']               = "isMediumMuon && abs(eta)<2.1 && pt > 5"
#config['MUON_TAG_CUTS']           = "(userInt('tightID')==1 && pt > 25 && abs(eta) < 2.1 && "+isolationDef+" < 0.2)"
config['MAXEVENTS']               = cms.untracked.int32(-1) 
config['INPUT_FILE_NAME']      	  = options.inputFiles
config['OUTPUT_FILE_NAME']        = options.outputFile
config['OUTPUTEDMFILENAME']       = 'debug.root'
config['DEBUG']                   = cms.bool(False)

from PhysicsTools.TagAndProbe.treeMakerOptions_cfi import *

if (options.isMC):
    config['TnPHLTTagFilters']    = cms.vstring() #"hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09")
    config['TnPHLTProbeFilters']  = cms.vstring()
    config['GLOBALTAG']           = 'auto:run2_mc'
else:
    config['TnPHLTTagFilters']    = "hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09"
    config['TnPHLTProbeFilters']  = cms.vstring() #"hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09" #cms.vstring()
    config['GLOBALTAG']           = 'auto:run2_data'
    config['json']                = 'Cert_271036-274421_13TeV_PromptReco_Collisions16_JSON.txt'


###################################################################

#setModules(process, options)
from PhysicsTools.TagAndProbe.treeContent_cfi import *

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, config['GLOBALTAG'], '')


process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(config['INPUT_FILE_NAME']),
			    firstRun  = cms.untracked.uint32(274094)
                            )

process.maxEvents = cms.untracked.PSet( input = config['MAXEVENTS'])

###################################################################
## ID
###################################################################
process.muonID = cms.EDProducer(
    "MiniAODMuonIDEmbedder",
    src = cms.InputTag(config['MUON_COLL']),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
)
process.idEmbedSequence = cms.Sequence(process.muonID)
muonSource = 'muonID'


############
### Tags ###
############
process.tagMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag(muonSource),

    cut = cms.string(config['MUON_TAG_CUTS']),
    filter = cms.bool(True)
)
process.tagMuonsTriggerMatched = cms.EDProducer("PatMuonTriggerCandProducer",
    filterNames = cms.vstring(config['TnPHLTTagFilters']),
    inputs      = cms.InputTag("tagMuons"),
    bits        = cms.InputTag('TriggerResults::HLT'),
    objects     = cms.InputTag('selectedPatTrigger'),
    dR          = cms.double(0.5),
    isAND       = cms.bool(True)
    )

process.probeMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag(muonSource),
    cut = cms.string(config['MUON_CUTS']), 
)

######################
### Trigger Probes ###
######################
process.probeTriggerSeq = cms.Sequence()
# IsoMu
process.probeTriggersMu17Leg = cms.EDProducer("PatMuonTriggerCandProducer",
    filterNames = cms.vstring("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4", "hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17"),
    inputs      = cms.InputTag("probeMuons"),
    bits        = cms.InputTag('TriggerResults::HLT'),
    objects     = cms.InputTag('selectedPatTrigger'),
    dR          = cms.double(0.5),
    isAND       = cms.bool(True)
    )
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


###################################################################
## TnP PAIRS
###################################################################

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuonsTriggerMatched@+ probeMuons@-"), # charge coniugate states are implied
    cut   = cms.string("40 < mass < 200")
)


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
pu_distribs = { "mc" : mix.input.nbPileupEvents.probValue }

data_pu_distribs = { "golden_v2" : [2.2566754964100775, 1950.7101413328778, 17511.456675455047, 27074.20916127319, 31295.38904638094, 53522.20431560162, 35220.460642050515, 90928.4505886674, 158090.51550539696, 244779.69121570754, 475446.24468129413, 1032094.7569565654, 2867278.6011298, 8427431.69028668, 17066351.233420365, 24027307.400081296, 28582023.0191919, 31161312.562541533, 29834230.38765881, 24016276.70437407, 16773532.053158931, 10934908.81338181, 6882987.985089145, 4078652.637076056, 2188705.6729904166, 1038844.4766897545, 433153.3500083798, 160781.35010697396, 56761.72710812498, 23088.051876630037, 13703.115979318482, 11160.425830582903, 10107.562092362476, 9265.663899467123, 8452.643876626442, 7708.148529163162, 7083.028042403028, 6599.59737809236, 6251.554972113813, 6013.330098038255, 5850.7412552672895, 5728.700884743008, 5612.778750900018, 5486.076394019967, 5335.204681008237, 5152.708506844515, 4936.573628932316, 4688.59937381231, 4413.05262576535, 4115.649709073904, 3802.8200587118135, 3481.1800882296766, 3157.154600542835, 2836.702560279824, 2525.1215607577033, 2226.9167983369393, 1945.7261608204267, 1684.2949902672606, 1444.4940213184038, 1227.3732210270105, 1033.2435191116354, 861.7780761805193, 712.1249032724896, 583.0232918504774, 472.9175407256343, 380.0627464092443, 302.61882216537015, 238.730305964887, 186.5908075849736, 144.49205632995265, 110.85839804233966, 84.26823511318716, 63.464311178194166, 47.35493604185581, 35.00826099036062, 25.641591111449614, 18.607501788495206, 13.378250655929659, 9.529677972026128, 6.725494537916399] }

process.pileupReweightingProducer = cms.EDProducer("PileupWeightProducer",
                                                   #hardcodedWeights = cms.untracked.bool(True),
                                                   pileupInfoTag    = cms.InputTag("slimmedAddPileupInfo"),
                                                   PileupMC = cms.vdouble(pu_distribs["mc"]),
                                                   PileupData = cms.vdouble(data_pu_distribs["golden_v2"]),
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
    probe_eta    = cms.string("eta"),
    probe_abseta = cms.string("abs(eta)"),
    probe_pt     = cms.string("pt"),
    probe_et     = cms.string("et"),
    probe_e      = cms.string("energy"),
    probe_q      = cms.string("charge"),
    )

TagVariablesToStore = cms.PSet(
    tag_eta    = cms.string("eta"),
    tag_abseta = cms.string("abs(eta)"),
    tag_pt     = cms.string("pt"),
    tag_et     = cms.string("et"),
    tag_e      = cms.string("energy"),
    tag_q      = cms.string("charge"),
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
    tagFlags       =  cms.PSet(),    
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

if (not options.isMC):
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(False)
        )


process.muonEffs = cms.EDAnalyzer("TagProbeFitTreeProducer",
    CommonStuffForMuonProbe, mcTruthCommonStuff,
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("Random2"),
    flags         = cms.PSet(
	passingHTTID  = cms.string("isMediumMuon && abs(userFloat('dXY')) < 0.045 && abs(userFloat('dZ')) < 0.2"),
        passingMedium = cms.string("isMediumMuon"),
        passingTight  = cms.string("userInt('tightID')==1"), 
        passingHTTIsoLoose = cms.string(isolationDef+" < 0.3"),
        passingHTTIsoTight = cms.string(isolationDef+" < 0.1"),
        passingIsoMu20 = cms.InputTag("probeTriggersIsoMu20"),
        passingIsoMu22 = cms.InputTag("probeTriggersIsoMu22"),
        passingIsoMu22Eta2p1 = cms.InputTag("probeTriggersIsoMu22Eta2p1"),
        passingIsoMu24 = cms.InputTag("probeTriggersIsoMu24"),
        passingIsoMu27 = cms.InputTag("probeTriggersIsoMu27"),
    ),
    allProbes     = cms.InputTag("probeMuons"),
    )

process.tpPairSeq = cms.Sequence(
    process.tpPairs
)

if options.isMC :
    process.tpPairSeq += process.muMcMatch
    #process.tpPairSeq += process.tpPairsMCEmbedded
    process.tpPairSeq += process.pileupReweightingProducer
    process.muonEffs.isMC = cms.bool(True)
    process.muonEffs.eventWeight   = cms.InputTag("generator")
    process.muonEffs.PUWeightSrc   = cms.InputTag("pileupReweightingProducer","pileupWeights")
    #setattr(process.muonEffs.pairVariables, 'mc_mass', cms.string("userFloat('mc_mass')"))
    process.muonEffs.tagProbePairs = cms.InputTag("tpPairs")

if not options.isMC :
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/'+config['json']).getVLuminosityBlockRange()

process.p = cms.Path(
    process.idEmbedSequence *
    (process.tagMuons + process.probeMuons) *
    (process.tagMuonsTriggerMatched + process.probeTriggerSeq) *
    process.tpPairSeq *
    process.muonEffs
    )

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string(config['OUTPUT_FILE_NAME']),
    closeFileFast = cms.untracked.bool(True)
    )

