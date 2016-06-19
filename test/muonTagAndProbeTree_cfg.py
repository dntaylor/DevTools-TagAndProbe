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
options['MUON_TAG_CUTS']           = "(userInt('isTightMuon')==1 && pt > 25 && abs(eta) < 2.1 && "+isolationDef+" < 0.2)"
options['MAXEVENTS']               = cms.untracked.int32(-1) 
options['OUTPUTEDMFILENAME']       = 'edmFile.root'
options['DEBUG']                   = cms.bool(False)

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
    options['TnPPATHS']            = ["HLT_IsoTkMu20_v*",]
    options['TnPHLTTagFilters']    = ["hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09"]
    #options['TnPHLTTagFilters']    = ["hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09"]
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

######################
### Trigger Probes ###
######################
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

#process.probeTriggersMu17LegL1Mu12 = cms.EDProducer("L1MuonMatcher",
#        inputs = cms.InputTag("probeTriggersMu17Leg"),
#        l1extraMuons = cms.InputTag("l1extraParticles"),
#        minET = cms.double(12.),
#        dRmatch = cms.double(.5)
#        )
#process.probeTriggerSeq += process.probeTriggersMu17LegL1Mu12

process.probeTriggersMu8Leg = process.probeTriggersMu17Leg.clone()
process.probeTriggersMu8Leg.filterNames = cms.vstring("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4", "hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8")
process.probeTriggerSeq += process.probeTriggersMu8Leg

process.probeTriggersTkMu8Leg = process.probeTriggersMu17Leg.clone()
process.probeTriggersTkMu8Leg.filterNames = cms.vstring("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4", "hltDiMuonGlbFiltered17TrkFiltered8")
process.probeTriggerSeq += process.probeTriggersTkMu8Leg

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
from SimGeneral.MixingModule.mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU_cfi import mix
pu_distribs = { "mc" : mix.input.nbPileupEvents.probValue }

data_pu_distribs = { "golden" : [2182.0974237815058, 23982.61536263803, 69965.65140406518, 197986.45804659263, 360761.2702934934, 618833.2110524588, 1294653.2625905501, 9002888.851089654, 23288262.214837648, 31117597.39099118, 40944017.230754964, 57733538.73658765, 83682315.57466365, 115990600.07824127, 150941145.77642596, 183785582.24765, 205622340.13201085, 210634879.14830402, 200928607.69355252, 181576740.9714111, 156292041.111647, 127532432.99693091, 97964897.13034214, 70450676.96016556, 47226312.58673038, 29430755.94919289, 17058024.92005088, 9227445.142189559, 4682639.080967543, 2240806.73464195, 1015167.8624996855, 436276.70215337543, 177856.6086062794, 68673.90271215088, 25057.616014526295, 8620.812303806977, 2792.0870938091107, 850.8624973240721, 244.21998405867603, 66.24071137170579, 17.095448894713563, 4.250745652690683, 1.0390089900679278, 0.2565958750195954, 0.06582988128238443, 0.017806844134965572, 0.005044486045194729, 0.0014636359007932143, 0.00042440737017085395, 0.00012072800218017043] }

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

if (not varOptions.isMC):
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(False)
        )


process.muonEffs = cms.EDAnalyzer("TagProbeFitTreeProducer",
    CommonStuffForMuonProbe, mcTruthCommonStuff,
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("Random2"),
    flags         = cms.PSet(
        passingLoose = cms.string("isLooseMuon"),
        passingMedium = cms.string("isMediumMuon"),
        passingTight  = cms.string("userInt('isTightMuon')==1"), 
        passingIsoLoose = cms.string(isolationDef+" < 0.4"),
        passingIsoTight = cms.string(isolationDef+" < 0.15"),
        passingTrackIso = cms.string(trackIsoDef+" < 0.4"),
        passingMu17 = cms.InputTag("probeTriggersMu17Leg"),
        #passingMu17L1Match = cms.InputTag("probeTriggersMu17LegL1Mu12"),
        passingMu8= cms.InputTag("probeTriggersMu8Leg"),
        passingTkMu8 = cms.InputTag("probeTriggersTkMu8Leg"),
        passingIsoMu18 = cms.InputTag("probeTriggersIsoMu18"),
        passingIsoMu20 = cms.InputTag("probeTriggersIsoMu20"),
        passingIsoMu22 = cms.InputTag("probeTriggersIsoMu22"),
        passingIsoMu22Eta2p1 = cms.InputTag("probeTriggersIsoMu22Eta2p1"),
        passingIsoMu24 = cms.InputTag("probeTriggersIsoMu24"),
        passingIsoMu27 = cms.InputTag("probeTriggersIsoMu27"),
        passingIsoTkMu18 = cms.InputTag("probeTriggersIsoTkMu18"),
        passingIsoTkMu20 = cms.InputTag("probeTriggersIsoTkMu20"),
        passingIsoTkMu22 = cms.InputTag("probeTriggersIsoTkMu22"),
        passingIsoTkMu22Eta2p1 = cms.InputTag("probeTriggersIsoTkMu22Eta2p1"),
        passingIsoTkMu24 = cms.InputTag("probeTriggersIsoTkMu24"),
        passingIsoTkMu27 = cms.InputTag("probeTriggersIsoTkMu27"),
    ),
    allProbes     = cms.InputTag("probeMuons"),
    )

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
    process.muonEffs.tagProbePairs = cms.InputTag("tpPairs")

#if not options.isMC :
#    import FWCore.PythonUtilities.LumiList as LumiList
#    process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/'+options['json']).getVLuminosityBlockRange()

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

