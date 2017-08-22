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

def fitBin(name, allProbeCondition, passingProbeCondition, tmc=None, tmcAlt=None, tdata=None, obj='photon'):
    fitVariable = ROOT.RooRealVar('eg_mass', 'TP Pair Mass', 60, 120, 'GeV')
    fitVariable.setBins(60)

    #mcTruthCondition = ['mcTrue']
    mcTruthCondition = []

    ROOT.gDirectory.mkdir(name).cd()
    fitter = PassFailSimulFitter(name, fitVariable)
    fitter.addDataFromTree(tmc, 'mcData', allProbeCondition+mcTruthCondition, passingProbeCondition, separatePassFail = True, weightVariable='genWeight*pileupWeight*triggerEfficiency*e_mediumScale')
    fitter.addDataFromTree(tmcAlt, 'mcAltData', allProbeCondition+mcTruthCondition, passingProbeCondition, separatePassFail = True, weightVariable='genWeight*pileupWeight*triggerEfficiency*e_mediumScale')
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
        #cAll = fitter.drawFitCanvas(fitResult)
        #python_mkdir('fits_{0}/allFits/{1}'.format(obj,name))
        #cAll.Print('fits_{0}/allFits/{1}/{1}_{2}.png'.format(obj,name, varName))

    ROOT.TNamed('cutString', cutString).Write('',ROOT.TObject.kOverwrite)
    print
    ROOT.gDirectory.cd('..')

def fitBinAlt(name, allProbeCondition, passingProbeCondition, tdata=None, obj='photon'):
    fitVariable = ROOT.RooRealVar('eg_mass', 'TP Pair Mass', 60, 120, 'GeV')
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

def fit(name, allProbeCondition, passingProbeCondition, binningMap, macroVariables, tmc=None, tmcAlt=None, tdata=None, alt=False):
    ROOT.gDirectory.mkdir(name).cd()
    ROOT.TNamed('variables', ', '.join(macroVariables)).Write('',ROOT.TObject.kOverwrite)
    for binName, cut in sorted(binningMap.items()) :
        if alt:
            fitBinAlt(name+'_'+binName, allProbeCondition+cut, passingProbeCondition, tdata=tdata)
        else:
            fitBin(name+'_'+binName, allProbeCondition+cut, passingProbeCondition, tmc=tmc, tmcAlt=tmcAlt, tdata=tdata)
    ROOT.gDirectory.cd('..')

def runfit(args):
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.Math.MinimizerOptions.SetDefaultTolerance(1.e-2) # default is 1.e-2

    # trees for mc, mc lo (for systematics), and data
    treeName = 'EGTree'

    fmc = ROOT.TFile.Open(args.mcFileName)
    tmc = fmc.Get(treeName)

    fmcAlt = ROOT.TFile.Open(args.mcLOFileName)
    tmcAlt = fmcAlt.Get(treeName)

    fdata = ROOT.TFile.Open(args.dataFileName)
    tdata = fdata.Get(treeName)

    # binning for the efficiencies
    ptBin = getBinning('photon','pt')
    etaBin = getBinning('photon','eta')

    ptVar = 'g_pt'
    etaVar = 'g_eta'

    binning = {}
    for pb in range(len(ptBin[:-1])):
        ptlow = ptBin[pb]
        pthigh = ptBin[pb+1]
        ptname = 'pt{0}to{1}'.format(ptlow,pthigh)
        ptcut = '{0}>={1} && {0}<{2}'.format(ptVar,ptlow,pthigh)
        for eb in range(len(etaBin[:-1])):
            etalow = etaBin[eb]
            etahigh = etaBin[eb+1]
            etaname = 'eta{0}to{1}'.format(etalow,etahigh)
            etacut = '{0}>={1} && {0}<{2}'.format(etaVar,etalow,etahigh)
            binning['{0}_{1}'.format(ptname,etaname)] = [ptcut,etacut]

    commonVars = ['float {0}'.format(ptVar), 'float {0}'.format(etaVar)]

    # run the fits
    fname = args.outputFileName
    fout = ROOT.TFile(fname, 'recreate')
    directory = 'PhotonFits'
    fout.mkdir(directory).cd()

    idArgs = {
        'Preselection' : {'condition': ['g_passElectronVeto<0.5'], 'variable': 'g_passPreselectionNoElectronVeto>0.5',                           'fitVars': ['float g_mvaNonTrigValues', 'bool g_passPreselectionNoElectronVeto','bool g_passElectronVeto']},
        'MVA0p0Pre'    : {'condition': ['g_passElectronVeto<0.5'], 'variable': 'g_mvaNonTrigValues>0.0 && g_passPreselectionNoElectronVeto>0.5', 'fitVars': ['float g_mvaNonTrigValues', 'bool g_passPreselectionNoElectronVeto','bool g_passElectronVeto']},
        'MVA0p0PreFail': {'condition': ['g_passElectronVeto<0.5'], 'variable': 'g_mvaNonTrigValues<0.0 && g_passPreselectionNoElectronVeto>0.5', 'fitVars': ['float g_mvaNonTrigValues', 'bool g_passPreselectionNoElectronVeto','bool g_passElectronVeto']},
        'MVA0p0'       : {'condition': ['g_passElectronVeto<0.5'], 'variable': 'g_mvaNonTrigValues>0.0',                                         'fitVars': ['float g_mvaNonTrigValues', 'bool g_passPreselectionNoElectronVeto','bool g_passElectronVeto']},
        'MVA0p0Fail'   : {'condition': ['g_passElectronVeto<0.5'], 'variable': 'g_mvaNonTrigValues<0.0',                                         'fitVars': ['float g_mvaNonTrigValues', 'bool g_passPreselectionNoElectronVeto','bool g_passElectronVeto']},
    }

    for idArg,vals in idArgs.iteritems():
        fit(idArg, vals['condition'], vals['variable'], binning, commonVars+vals['fitVars'], tmc=tmc, tmcAlt=tmcAlt, tdata=tdata)



def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='TagAndProbe Fitter')

    parser.add_argument('--mcFileName', '-mc', type=str, default='TnPTree_mc.root', help='Filename for MC TagAndProbe Tree')
    parser.add_argument('--mcLOFileName', '-mcLO', type=str, default='TnPTree_mcLO.root', help='Filename for MC LO TagAndProbe Tree')
    parser.add_argument('--dataFileName', '-data', type=str, default='TnPTree_data.root', help='Filename for Data TagAndProbe Tree')
    parser.add_argument('--outputFileName', type=str, default='fits_photon.root', help='Filename for output')

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

