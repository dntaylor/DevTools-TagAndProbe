#!/usr/bin/env python
import os
import sys
import logging
import argparse

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
from DevTools.TagAndProbe.PassFailSimulFitter import PassFailSimulFitter
from DevTools.Utilities.utilities import python_mkdir
from DevTools.TagAndProbe.utilities import getBinning

logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

pdfDefinition = []
with open('{0}/src/DevTools/TagAndProbe/data/pdfDefinitions.txt'.format(os.environ['CMSSW_BASE'])) as defFile :
    for line in defFile :
        line = line.strip()
        if len(line) == 0 or line[0] is '#' :
            continue
        pdfDefinition.append(line)

pdfDefinitionAlt = []
with open('{0}/src/DevTools/TagAndProbe/data/pdfDefinitions_alt.txt'.format(os.environ['CMSSW_BASE'])) as defFile :
    for line in defFile :
        line = line.strip()
        if len(line) == 0 or line[0] is '#' :
            continue
        pdfDefinitionAlt.append(line)

def statusInfo(fitResults):
    fitStatus=':'.join(['% d' % fitResults.statusCodeHistory(i) for i in range(fitResults.numStatusHistory())]),
    return fitStatus

def fitBin(name, allProbeCondition, passingProbeCondition, tmc=None, tmcAlt=None, tdata=None, obj=''):
    fitVariable = ROOT.RooRealVar('mass', 'TP Pair Mass', 60, 120, 'GeV')
    fitVariable.setBins(60)

    #mcTruthCondition = ['mcTrue']
    mcTruthCondition = []

    ROOT.gDirectory.mkdir(name).cd()
    fitter = PassFailSimulFitter(name, fitVariable)
    fitter.addDataFromTree(tmc, 'mcData', allProbeCondition+mcTruthCondition, passingProbeCondition, separatePassFail = True)
    fitter.addDataFromTree(tmcAlt, 'mcAltData', allProbeCondition+mcTruthCondition, passingProbeCondition, separatePassFail = True)
    nMCPass = fitter.workspace.data('mcDataPass').sumEntries()
    nMCFail = fitter.workspace.data('mcDataFail').sumEntries()
    if nMCPass!=nMCPass or nMCFail!=nMCFail:
        print 'WARNING NaN: {0}'.format(name)
    mcEff = nMCPass/(nMCPass+nMCFail) if nMCPass+nMCFail else 0.
    mcEffLo = ROOT.TEfficiency.ClopperPearson(int(nMCPass+nMCFail), int(nMCPass), 0.68, False)
    mcEffHi = ROOT.TEfficiency.ClopperPearson(int(nMCPass+nMCFail), int(nMCPass), 0.68, True)
    h=ROOT.TH1F('mc_cutCount', 'Cut & Count', 2, 0, 2)
    h.SetBinContent(1, nMCPass)
    h.SetBinContent(2, nMCPass+nMCFail)

    # All MC templates must be set up by now
    fitter.setPdf(pdfDefinition)

    print '-'*40, 'Central value fit'
    fitter.addDataFromTree(tdata, 'data', allProbeCondition, passingProbeCondition, weightVariable='1')
    res = fitter.fit('simPdf', 'data')
    effValue = res.floatParsFinal().find('efficiency')
    dataEff = effValue.getVal()
    dataEffErrHi = effValue.getErrorHi()
    dataEffErrLo = effValue.getErrorLo()
    scaleFactor = dataEff / mcEff if mcEff else 0.
    maxSf = (dataEff+dataEffErrHi)/mcEffLo if mcEffLo else 0.
    minSf = (dataEff+dataEffErrLo)/mcEffHi if mcEffHi else 0.
    res.SetName('fitresults')
    c = fitter.drawFitCanvas(res)
    c.Write('',ROOT.TObject.kOverwrite)
    h.Write('',ROOT.TObject.kOverwrite)
    res.Write('',ROOT.TObject.kOverwrite)

    print '-'*40, 'Fit with alternate MC template'
    resAlt = fitter.fit('simAltPdf', 'data')
    dataAltEff = resAlt.floatParsFinal().find('efficiency').getVal()
    resAlt.SetName('fitresults_systAltTemplate')
    resAlt.Write('',ROOT.TObject.kOverwrite)

    #print '-'*40, 'Fit with tag pt > 30 (vs. 25)'
    #fitter.addDataFromTree(tdata, 'dataTagPt30', allProbeCondition+['tag_Ele_pt>30'], passingProbeCondition, weightVariable='1')
    #resTagPt30 = fitter.fit('simPdf', 'dataTagPt30')
    #dataTagPt30Eff = resTagPt30.floatParsFinal().find('efficiency').getVal()
    #resTagPt30.Write('',ROOT.TObject.kOverwrite)

    print '-'*40, 'Fit with CMSShape background (vs. Bernstein)'
    resCMSBkg = fitter.fit('simCMSBkgPdf', 'data')
    dataCMSBkgEff = resCMSBkg.floatParsFinal().find('efficiency').getVal()
    resCMSBkg.Write('',ROOT.TObject.kOverwrite)

    fitter.workspace.Write('',ROOT.TObject.kOverwrite)
    print name, ': Data=%.2f, MC=%.2f, Ratio=%.2f' % (dataEff, mcEff, scaleFactor)
    condition = ' && '.join(allProbeCondition+[passingProbeCondition])
    variations = {
            'CENTRAL'  : (scaleFactor, res),
            'STAT_UP'  : (maxSf, res),
            'STAT_DOWN': (minSf, res),
            'SYST_ALT_TEMPL' : (dataAltEff / mcEff if mcEff else 0., resAlt),
            #'SYST_TAG_PT30' : (dataTagPt30Eff / mcEff if mcEff else 0., resTagPt30),
            'SYST_CMSSHAPE' : (dataCMSBkgEff / mcEff if mcEff else 0., resCMSBkg),
            'EFF_DATA' : (dataEff, res),
            'EFF_DATA_ERRSYM' : ((dataEffErrHi-dataEffErrLo)/2, res),
            'EFF_MC' : (mcEff, res),
            'EFF_MC_ERRSYM' : ((mcEffHi-mcEffLo)/2, res),
            }
    cutString = ''
    for varName, value in variations.items() :
        (value, fitResult) = value
        cutString += '    if ( variation == Variation::%s && (%s) ) return %f;\n' % (varName, condition, value)
        print '  Variation {:>15s} : {:.4f}, edm={:f}, status={:s}'.format(varName, value, fitResult.edm(), statusInfo(fitResult))
        if 'STAT' not in varName and 'EFF' not in varName and fitResult.statusCodeHistory(0) < 0 :
            cBad = fitter.drawFitCanvas(fitResult)
            python_mkdir('fits_{0}/badFits/{1}'.format(obj,name))
            cBad.Print('fits_{0}/badFits/{1}/badFit_{1}_{2}.png'.format(obj,name, varName))

    ROOT.TNamed('cutString', cutString).Write('',ROOT.TObject.kOverwrite)
    print
    ROOT.gDirectory.cd('..')

def fitBinAlt(name, allProbeCondition, passingProbeCondition, tdata=None, obj=''):
    fitVariable = ROOT.RooRealVar('mass', 'TP Pair Mass', 60, 120, 'GeV')
    fitVariable.setBins(60)

    ROOT.gDirectory.mkdir(name).cd()
    fitter = PassFailSimulFitter(name, fitVariable)

    fitter.setPdf(pdfDefinitionAlt)

    print '-'*40, 'Central value fit'
    fitter.addDataFromTree(tdata, 'data', allProbeCondition, passingProbeCondition, weightVariable='1')
    effValue = fitter.fit('simPdf', 'data', doCutAndCount=True)
    #effValue = res.floatParsFinal().find('efficiency')
    res = None
    dataEff = effValue.getVal()
    dataEffErrHi = effValue.getErrorHi()
    dataEffErrLo = effValue.getErrorLo()
    #res.SetName('fitresults')
    #c = fitter.drawFitCanvas(res,skip=True)
    #c.Write('',ROOT.TObject.kOverwrite)
    #res.Write('',ROOT.TObject.kOverwrite)

    #print '-'*40, 'Fit with tag pt > 30 (vs. 25)'
    #fitter.addDataFromTree(tdata, 'dataTagPt30', allProbeCondition+['tag_Ele_pt>30'], passingProbeCondition, weightVariable='1')
    #resTagPt30 = fitter.fit('simPdf', 'dataTagPt30')
    #dataTagPt30Eff = resTagPt30.floatParsFinal().find('efficiency').getVal()
    #resTagPt30.Write('',ROOT.TObject.kOverwrite)

    #print '-'*40, 'Fit with CMSShape background (vs. Bernstein)'
    #resCMSBkg = fitter.fit('simCMSBkgPdf', 'data')
    #dataCMSBkgEff = resCMSBkg.floatParsFinal().find('efficiency').getVal()
    #resCMSBkg.Write('',ROOT.TObject.kOverwrite)

    fitter.workspace.Write('',ROOT.TObject.kOverwrite)
    print name, ': Data=%.2f' % (dataEff)
    condition = ' && '.join(allProbeCondition+[passingProbeCondition])
    variations = {
            'CENTRAL'  : (dataEff, res),
            'STAT_UP'  : (dataEff+dataEffErrHi, res),
            'STAT_DOWN': (dataEff-dataEffErrLo, res),
            #'SYST_TAG_PT30' : (dataTagPt30Eff / dataEff if dataEff else 0., resTagPt30),
            #'SYST_CMSSHAPE' : (dataCMSBkgEff / dataEff if dataEff else 0., resCMSBkg),
            'EFF_DATA' : (dataEff, res),
            'EFF_DATA_ERRSYM' : ((dataEffErrHi-dataEffErrLo)/2, res),
            }
    cutString = ''
    for varName, value in variations.items() :
        (value, fitResult) = value
        cutString += '    if ( variation == Variation::%s && (%s) ) return %f;\n' % (varName, condition, value)
        #print '  Variation {:>15s} : {:.4f}, edm={:f}, status={:s}'.format(varName, value, fitResult.edm(), statusInfo(fitResult))
        #if 'STAT' not in varName and 'EFF' not in varName and fitResult.statusCodeHistory(0) < 0 :
        #    cBad = fitter.drawFitCanvas(fitResult,skip=True)
        #    python_mkdir('fits_{0}/badFits/{1}'.format(obj,name))
        #    cBad.Print('fits_{0}/badFits/{1}/badFit_{1}_{2}.png'.format(obj,name, varName))

    ROOT.TNamed('cutString', cutString).Write('',ROOT.TObject.kOverwrite)
    print
    ROOT.gDirectory.cd('..')

def fit(name, allProbeCondition, passingProbeCondition, binningMap, macroVariables, tmc=None, tmcAlt=None, tdata=None, obj='', alt=False):
    ROOT.gDirectory.mkdir(name).cd()
    ROOT.TNamed('variables', ', '.join(macroVariables)).Write('',ROOT.TObject.kOverwrite)
    for binName, cut in sorted(binningMap.items()) :
        if alt:
            fitBinAlt(name+'_'+binName, allProbeCondition+cut, passingProbeCondition, tdata=tdata, obj=obj)
        else:
            fitBin(name+'_'+binName, allProbeCondition+cut, passingProbeCondition, tmc=tmc, tmcAlt=tmcAlt, tdata=tdata, obj=obj)
    ROOT.gDirectory.cd('..')

def runfit(args):
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.Math.MinimizerOptions.SetDefaultTolerance(1.e-2) # default is 1.e-2

    # trees for mc, mc lo (for systematics), and data
    treeNameMap = {
        'electron' : 'GsfElectronToEleID/fitter_tree',
        'electronTrig' : 'GsfElectronToTrigger/fitter_tree',
        'muon'     : 'muonEffs/fitter_tree',
    }

    fmc = ROOT.TFile.Open(args.mcFileName)
    tmc = fmc.Get(treeNameMap[args.object])

    fmcAlt = ROOT.TFile.Open(args.mcLOFileName)
    tmcAlt = fmcAlt.Get(treeNameMap[args.object])

    fdata = ROOT.TFile.Open(args.dataFileName)
    tdata = fdata.Get(treeNameMap[args.object])

    if args.object=='electron':
        tdata_trig = fdata.Get(treeNameMap[args.object+'Trig'])
    if args.object=='muon':
        tdata_trig = tdata


    # binning for the efficiencies
    ptBinMap = {
        'electron' : getBinning('electron','pt'),
        'muon'     : getBinning('muon','pt'),
        'electronTrig' : getBinning('electron','pt',trig=True),
        'muonTrig'     : getBinning('muon','pt',trig=True),
    }

    etaBinMap = {
        'electron' : getBinning('electron','eta'),
        'muon'     : getBinning('muon','eta'),
        'electronTrig' : getBinning('electron','eta',trig=True),
        'muonTrig'     : getBinning('muon','eta',trig=True),
    }

    # single bin
    #ptBinMap = {
    #    'electron' : [10,200],
    #    'muon'     : [10,200],
    #}

    #etaBinMap = {
    #    'electron' : [-2.5,2.5],
    #    'muon'     : [-2.4,2.4],
    #}

    ptVar = {
        'electron'     : 'probe_Ele_pt',
        'electronTrig' : 'probe_Ele_pt',
        'muon'         : 'probe_pt',
        'muonTrig'     : 'probe_pt',
    }

    etaVar = {
        'electron'     : 'probe_Ele_eta',
        'electronTrig' : 'probe_sc_eta',
        'muon'         : 'probe_eta',
        'muonTrig'     : 'probe_eta',
    }

    binning = {}
    for pb in range(len(ptBinMap[args.object][:-1])):
        ptlow = ptBinMap[args.object][pb]
        pthigh = ptBinMap[args.object][pb+1]
        ptname = 'pt{0}to{1}'.format(ptlow,pthigh)
        ptcut = '{0}>={1} && {0}<{2}'.format(ptVar[args.object],ptlow,pthigh)
        for eb in range(len(etaBinMap[args.object][:-1])):
            etalow = etaBinMap[args.object][eb]
            etahigh = etaBinMap[args.object][eb+1]
            etaname = 'eta{0}to{1}'.format(etalow,etahigh)
            etacut = '{0}>={1} && {0}<{2}'.format(etaVar[args.object],etalow,etahigh)
            binning['{0}_{1}'.format(ptname,etaname)] = [ptcut,etacut]

    binning_trig = {}
    for pb in range(len(ptBinMap[args.object+'Trig'][:-1])):
        ptlow = ptBinMap[args.object+'Trig'][pb]
        pthigh = ptBinMap[args.object+'Trig'][pb+1]
        ptname = 'pt{0}to{1}'.format(ptlow,pthigh)
        ptcut_trig = '{0}>={1} && {0}<{2}'.format(ptVar[args.object+'Trig'],ptlow,pthigh)
        for eb in range(len(etaBinMap[args.object+'Trig'][:-1])):
            etalow = etaBinMap[args.object+'Trig'][eb]
            etahigh = etaBinMap[args.object+'Trig'][eb+1]
            etaname = 'eta{0}to{1}'.format(etalow,etahigh)
            etacut_trig = '{0}>={1} && {0}<{2}'.format(etaVar[args.object+'Trig'],etalow,etahigh)
            binning_trig['{0}_{1}'.format(ptname,etaname)] = [ptcut_trig,etacut_trig]

    commonVars = ['float {0}'.format(ptVar[args.object]), 'float {0}'.format(etaVar[args.object])]
    commonVars_trig = ['float {0}'.format(ptVar[args.object+'Trig']), 'float {0}'.format(etaVar[args.object+'Trig'])]

    # run the fits
    fname = args.outputFileName
    fout = ROOT.TFile(fname, 'recreate')
    directory = '{0}Fits'.format(args.object)
    fout.mkdir(directory).cd()

    idArgs = {
        'electron': {
            'CutBasedIDVeto':     {'condition': [], 'variable': 'passingVeto',   'fitVars': ['bool passingVeto']},
            'CutBasedIDLoose':    {'condition': [], 'variable': 'passingLoose',  'fitVars': ['bool passingLoose']},
            'CutBasedIDMedium':   {'condition': [], 'variable': 'passingMedium', 'fitVars': ['bool passingMedium']},
            'CutBasedIDTight':    {'condition': [], 'variable': 'passingTight',  'fitVars': ['bool passingTight']},
        },
        'muon': {
            'HppLooseID':                {'conditions': [],                     'variable': 'passingHppLooseID',  'fitVars': ['bool passingHppLooseID']},
            'HppLooseIsoFromLooseID':    {'conditions': ['passingHppLooseID'],  'variable': 'passingHppLooseIso', 'fitVars': ['bool passingHppLooseID','bool passingHppLooseIso']},
            'HppMediumID':               {'conditions': [],                     'variable': 'passingHppMediumID', 'fitVars': ['bool passingHppMediumID']},
            'HppMediumIsoFromMediumID':  {'conditions': ['passingHppMediumID'], 'variable': 'passingHppMediumIso','fitVars': ['bool passingHppMediumID','bool passingHppMediumIso']},
            'LooseID':                   {'conditions': [],                     'variable': 'passingLoose',       'fitVars': ['bool passingLoose']},
            'LooseIsoFromLooseID':       {'conditions': ['passingLoose'],       'variable': 'probe_isoR04<0.25',  'fitVars': ['bool passingLoose','float probe_isoR04']},
            'MediumID':                  {'conditions': [],                     'variable': 'passingMedium',      'fitVars': ['bool passingMedium']},
            'LooseIsoFromMediumID':      {'conditions': ['passingMedium'],      'variable': 'probe_isoR04<0.25',  'fitVars': ['bool passingMedium','float probe_isoR04']},
            'TightIsoFromMediumID':      {'conditions': ['passingMedium'],      'variable': 'probe_isoR04<0.15',  'fitVars': ['bool passingMedium','float probe_isoR04']},
            'MediumIDICHEP':             {'conditions': [],                     'variable': 'passingMediumICHEP', 'fitVars': ['bool passingMediumICHEP']},
            'LooseIsoFromMediumIDICHEP': {'conditions': ['passingMediumICHEP'], 'variable': 'probe_isoR04<0.25',  'fitVars': ['bool passingMediumICHEP','float probe_isoR04']},
            'TightIsoFromMediumIDICHEP': {'conditions': ['passingMediumICHEP'], 'variable': 'probe_isoR04<0.15',  'fitVars': ['bool passingMediumICHEP','float probe_isoR04']},
            'TightID':                   {'conditions': [],                     'variable': 'passingTight',       'fitVars': ['bool passingTight']},
            'TightIsoFromTightID':       {'conditions': ['passingTight'],       'variable': 'probe_isoR04<0.15',  'fitVars': ['bool passingTight','float probe_isoR04']},

        },
    }

    trigArgs = {
        'electron': {
            'Ele27WPTight':       {'condition': [],                                               'variable': 'passingEle27WPTight',       'fitVars': ['bool passingEle27WPTight']},
            'Ele27Eta2p1WPTight': {'condition': [],                                               'variable': 'passingEle27Eta2p1WPTight', 'fitVars': ['bool passingEle27Eta2p1WPTight']},
            'Ele23Leg':           {'condition': [],                                               'variable': 'passingEle23Leg',           'fitVars': ['bool passingEle23Leg']},
            'Ele12Leg':           {'condition': ['tag_passingEle23Ele12Leg1'],                    'variable': 'passingEle12Leg',           'fitVars': ['bool tag_passingEle23Ele12Leg1','bool passingEle12Leg']},
            'Ele12LegDZ':         {'condition': ['tag_passingEle23Ele12Leg1 && passingEle12Leg'], 'variable': 'passingEle12LegDZ',         'fitVars': ['bool tag_passingEle23Ele12Leg1','bool passingEle12Leg','bool passingEle12LegDZ']},
        },
        'muon': {
            'IsoMu24':                  {'condition': [],                                          'variable': 'passingIsoMu24',                   'fitVars': ['bool passingIsoMu24']},
            'IsoTkMu24':                {'condition': [],                                          'variable': 'passingIsoTkMu24',                 'fitVars': ['bool passingIsoTkMu24']},
            'IsoMu24ORIsoTkMu24':       {'condition': [],                                          'variable': 'passingIsoMu24ORIsoTkMu24',        'fitVars': ['bool passingIsoMu24ORIsoTkMu24']},
            'Mu50':                     {'condition': [],                                          'variable': 'passingMu50',                      'fitVars': ['bool passingMu50']},
            'IsoMu24ORIsoTkMu24ORMu50': {'condition': [],                                          'variable': 'passingIsoMu24ORIsoTkMu24ORMu50',  'fitVars': ['bool passingIsoMu24ORIsoTkMu24ORMu50']},
            'Mu17Leg':                  {'condition': [],                                          'variable': 'passingMu17Leg',                   'fitVars': ['bool passingMu17Leg']},
            'Mu8Leg':                   {'condition': ['tag_passingMu17'],                         'variable': 'passingMu8Leg',                    'fitVars': ['bool tag_passingMu17','bool passingMu8Leg']},
            'TkMu8Leg':                 {'condition': ['tag_passingMu17'],                         'variable': 'passingTkMu8Leg',                  'fitVars': ['bool tag_passingMu17','bool passingTkMu8Leg']},
            'Mu8ORTkMu8Leg':            {'condition': ['tag_passingMu17'],                         'variable': 'passingMu8ORTkMu8Leg',             'fitVars': ['bool tag_passingMu17','bool passingMu8ORTkMu8Leg']},
            'Mu8LegDZ':                 {'condition': ['tag_passingMu17 && passingMu8Leg'],        'variable': 'passingMu8LegDZ',                  'fitVars': ['bool tag_passingMu17','bool passingMu8Leg','bool passingMu8LegDZ']},
            'TkMu8LegDZ':               {'condition': ['tag_passingMu17 && passingTkMu8Leg'],      'variable': 'passingTkMu8LegDZ',                'fitVars': ['bool tag_passingMu17','bool passingTkMu8Leg','bool passingTkMu8LegDZ']},
            'Mu8ORTkMu8LegDZ':          {'condition': ['tag_passingMu17 && passingMu8ORTkMu8Leg'], 'variable': 'passingMu8ORTkMu8LegDZ',           'fitVars': ['bool tag_passingMu17','bool passingMu8ORTkMu8Leg','bool passingMu8ORTkMu8LegDZ']},
        },
    }

    if not args.fitsToRun:
        args.fitsToRun = idArgs[args.object].keys()+trigArgs[args.object].keys()

    for trigArg,vals in trigArgs[args.object].iteritems():
        if trigArg in args.fitsToRun: fit(trigArg, vals['condition'], vals['variable'], binning_trig, commonVars_trig+vals['fitVars'], tmc=tmc, tmcAlt=tmcAlt, tdata=tdata_trig, obj=args.object)
    for idArg,vals in idArgs[args.object].iteritems():
        if idArg in args.fitsToRun: fit(idArg, vals['condition'], vals['variable'], binning, commonVars+vals['fitVars'], tmc=tmc, tmcAlt=tmcAlt, tdata=tdata, obj=args.object)



def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='TagAndProbe Fitter')

    parser.add_argument('object', type=str, choices=['electron','muon'], help='Physics object')
    parser.add_argument('--mcFileName', '-mc', type=str, default='TnPTree_mc.root', help='Filename for MC TagAndProbe Tree')
    parser.add_argument('--mcLOFileName', '-mcLO', type=str, default='TnPTree_mcLO.root', help='Filename for MC LO TagAndProbe Tree')
    parser.add_argument('--dataFileName', '-data', type=str, default='TnPTree_data.root', help='Filename for Data TagAndProbe Tree')
    parser.add_argument('--fitsToRun', type=str, nargs='*', help='Run selected fits')
    parser.add_argument('--outputFileName', type=str, default='fits.root', help='Filename for output')

    return parser.parse_args(argv)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    runfit(args)

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)

