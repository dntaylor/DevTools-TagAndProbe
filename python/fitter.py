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

def statusInfo(fitResults):
    fitStatus=':'.join(['% d' % fitResults.statusCodeHistory(i) for i in range(fitResults.numStatusHistory())]),
    return fitStatus

def fitBin(name, allProbeCondition, passingProbeCondition, tmc=None, tmcAlt=None, tdata=None, obj=''):
    fitVariable = ROOT.RooRealVar('mass', 'TP Pair Mass', 60, 120, 'GeV')
    fitVariable.setBins(60)

    mcTruthCondition = ['mcTrue']
    mcTruthCondition = []

    ROOT.gDirectory.mkdir(name).cd()
    fitter = PassFailSimulFitter(name, fitVariable)
    fitter.addDataFromTree(tmc, 'mcData', allProbeCondition+mcTruthCondition, passingProbeCondition, separatePassFail = True)
    fitter.addDataFromTree(tmcAlt, 'mcAltData', allProbeCondition+mcTruthCondition, passingProbeCondition, separatePassFail = True)
    nMCPass = fitter.workspace.data('mcDataPass').sumEntries()
    nMCFail = fitter.workspace.data('mcDataFail').sumEntries()
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
    c.Write()
    h.Write()
    res.Write()

    print '-'*40, 'Fit with alternate MC template'
    resAlt = fitter.fit('simAltPdf', 'data')
    dataAltEff = resAlt.floatParsFinal().find('efficiency').getVal()
    resAlt.SetName('fitresults_systAltTemplate')
    resAlt.Write()

    print '-'*40, 'Fit with tag pt > 30 (vs. 25)'
    fitter.addDataFromTree(tdata, 'dataTagPt30', allProbeCondition+['tag_Ele_pt>30'], passingProbeCondition, weightVariable='1')
    resTagPt30 = fitter.fit('simPdf', 'dataTagPt30')
    dataTagPt30Eff = resTagPt30.floatParsFinal().find('efficiency').getVal()
    resTagPt30.Write()

    print '-'*40, 'Fit with CMSShape background (vs. Bernstein)'
    resCMSBkg = fitter.fit('simCMSBkgPdf', 'data')
    dataCMSBkgEff = resCMSBkg.floatParsFinal().find('efficiency').getVal()
    resCMSBkg.Write()

    fitter.workspace.Write()
    print name, ': Data=%.2f, MC=%.2f, Ratio=%.2f' % (dataEff, mcEff, scaleFactor)
    condition = ' && '.join(allProbeCondition+[passingProbeCondition])
    variations = {
            'CENTRAL'  : (scaleFactor, res),
            'STAT_UP'  : (maxSf, res),
            'STAT_DOWN': (minSf, res),
            'SYST_ALT_TEMPL' : (dataAltEff / mcEff if mcEff else 0., resAlt),
            'SYST_TAG_PT30' : (dataTagPt30Eff / mcEff if mcEff else 0., resTagPt30),
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

    ROOT.TNamed('cutString', cutString).Write()
    print
    ROOT.gDirectory.cd('..')

def fit(name, allProbeCondition, passingProbeCondition, binningMap, macroVariables, tmc=None, tmcAlt=None, tdata=None, obj=''):
    ROOT.gDirectory.mkdir(name).cd()
    ROOT.TNamed('variables', ', '.join(macroVariables)).Write()
    for binName, cut in sorted(binningMap.items()) :
        fitBin(name+'_'+binName, allProbeCondition+cut, passingProbeCondition, tmc=tmc, tmcAlt=tmcAlt, tdata=tdata, obj=obj)
    ROOT.gDirectory.cd('..')

def runfit(args):
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.Math.MinimizerOptions.SetDefaultTolerance(1.e-2) # default is 1.e-2

    # trees for mc, mc lo (for systematics), and data
    treeNameMap = {
        'electron' : 'GsfElectronToRECO/fitter_tree',
        'muon'     : 'muonEffs/fitter_tree',
    }

    fmc = ROOT.TFile.Open(args.mcFileName)
    tmc = fmc.Get(treeNameMap[args.object])

    fmcAlt = ROOT.TFile.Open(args.mcLOFileName)
    tmcAlt = fmcAlt.Get(treeNameMap[args.object])

    fdata = ROOT.TFile.Open(args.dataFileName)
    tdata = fdata.Get(treeNameMap[args.object])

    # binning for the efficiencies
    ptBinMap = {
        'electron' : getBinning('electron','pt'),
        'muon'     : getBinning('muon','pt'),
    }

    etaBinMap = {
        'electron' : getBinning('electron','eta'),
        'muon'     : getBinning('muon','eta'),
    }

    ptVar = {
        'electron' : 'probe_Ele_pt',
        'muon'     : 'probe_pt',
    }

    etaVar = {
        'electron' : 'probe_Ele_eta',
        'muon'     : 'probe_eta',
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

    commonVars = ['float {0}'.format(ptVar[args.object]), 'float {0}'.format(etaVar[args.object])]

    # run the fits
    fout = ROOT.TFile('fits_{0}.root'.format(args.object), 'recreate')
    fout.mkdir('{0}Fits'.format(args.object)).cd()

    if args.object=='electron':
        fit('CutBasedIDVeto',   [], 'passingVeto',   binning, commonVars+['bool passingVeto'],   tmc=tmc, tmcAlt=tmcAlt, tdata=tdata, obj=args.object)
        fit('CutBasedIDLoose',  [], 'passingLoose',  binning, commonVars+['bool passingLoose'],  tmc=tmc, tmcAlt=tmcAlt, tdata=tdata, obj=args.object)
        fit('CutBasedIDMedium', [], 'passingMedium', binning, commonVars+['bool passingMedium'], tmc=tmc, tmcAlt=tmcAlt, tdata=tdata, obj=args.object)
        fit('CutBasedIDTight',  [], 'passingTight',  binning, commonVars+['bool passingTight'],  tmc=tmc, tmcAlt=tmcAlt, tdata=tdata, obj=args.object)
    if args.object=='muon':
        fit('MediumID',   [], 'passingMedium',   binning, commonVars+['bool passingMedium'],   tmc=tmc, tmcAlt=tmcAlt, tdata=tdata, obj=args.object)
        fit('TightID',    [], 'passingTight',    binning, commonVars+['bool passingTight'],    tmc=tmc, tmcAlt=tmcAlt, tdata=tdata, obj=args.object)


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='TagAndProbe Fitter')

    parser.add_argument('object', type=str, choices=['electron','muon'], help='Physics object')
    parser.add_argument('--mcFileName', '-mc', type=str, default='TnPTree_mc.root', help='Filename for MC TagAndProbe Tree')
    parser.add_argument('--mcLOFileName', '-mcLO', type=str, default='TnPTree_mcLO.root', help='Filename for MC LO TagAndProbe Tree')
    parser.add_argument('--dataFileName', '-data', type=str, default='TnPTree_data.root', help='Filename for Data TagAndProbe Tree')

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

