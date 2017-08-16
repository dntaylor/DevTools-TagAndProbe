#!/usr/bin/env python
import array
import sys
import os
import subprocess
import argparse

from DevTools.TagAndProbe.utilities import getBinning
from DevTools.Plotter.utilities import getLumi
from DevTools.Utilities.utilities import python_mkdir

__gitversion__ = subprocess.check_output(["git", "describe", "--always"]).strip()

import ROOT
import DevTools.Plotter.CMS_lumi as CMS_lumi
import DevTools.Plotter.tdrstyle as tdrstyle
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPalette(1)

def setStyle(pad,position=11,preliminary=True):
    '''Set style for plots based on the CMS TDR style guidelines.
       https://twiki.cern.ch/twiki/bin/view/CMS/Internal/PubGuidelines#Figures_and_tables
       https://ghm.web.cern.ch/ghm/plots/'''
    # set period (used in CMS_lumi)
    # period : sqrts
    # 1 : 7, 2 : 8, 3 : 7+8, 4 : 13, ... 7 : 7+8+13
    period_int = 4
    # set position
    # 11: upper left, 33 upper right
    CMS_lumi.wrtieExtraText = preliminary
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (float(getLumi())/1000.)
    if getLumi() < 1000:
        CMS_lumi.lumi_13TeV = "%0.1f pb^{-1}" % (float(getLumi()))
    CMS_lumi.CMS_lumi(pad,period_int,position)


def save2D(eff,savename,outdir):
    canvas = ROOT.TCanvas(savename,savename,50,50,600,600)
    canvas.SetLogx(True)
    canvas.SetRightMargin(0.18) # for Z axis
    eff.GetXaxis().SetTitleOffset(1.4)
    eff.GetYaxis().SetTitleOffset(1.4)
    eff.GetXaxis().SetTitleSize(0.045)
    eff.GetYaxis().SetTitleSize(0.045)
    eff.Draw('colz')
    setStyle(canvas,position=0)
    for t in ['png','pdf']:
        name = '{0}/{2}/{1}.{2}'.format(outdir,savename,t)
        python_mkdir(os.path.dirname(name))
        canvas.Print(name)


def plot(args):

    trigArgs = {
        'electron': {
            'Ele27WPTight':       25,
            'Ele27Eta2p1WPTight': 25,
            'Ele23Leg':           25,
            'Ele12Leg':           10,
            'Ele12LegDZ':         10,
        },
        'muon': {
            'IsoMu24':                  25,
            'IsoTkMu24':                25,
            'IsoMu24ORIsoTkMu24':       25,
            'Mu50':                     50,
            'IsoMu24ORIsoTkMu24ORMu50': 25,
            'Mu17Leg':                  15,
            'Mu8Leg':                   10,
            'TkMu8Leg':                 10,
            'Mu8ORTkMu8Leg':            10,
            'Mu8LegDZ':                 10,
            'TkMu8LegDZ':               10,
            'Mu8ORTkMu8LegDZ':          10,
        },
    }

    if args.trig:
        ptbins = getBinning(args.object,'pt',trig=args.trig,threshold=trigArgs[args.object][args.idName])
    else:
        ptbins = getBinning(args.object,'pt',trig=args.trig)
    etabins = getBinning(args.object,'eta',trig=args.trig)
    colors = [
        ROOT.kRed, ROOT.kGreen, ROOT.kBlue, ROOT.kBlack, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kGreen+2, ROOT.kRed-3, ROOT.kCyan+1, ROOT.kMagenta-3, ROOT.kViolet-1, ROOT.kSpring+10,
        ROOT.kRed, ROOT.kGreen, ROOT.kBlue, ROOT.kBlack, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kGreen+2, ROOT.kRed-3, ROOT.kCyan+1, ROOT.kMagenta-3, ROOT.kViolet-1, ROOT.kSpring+10,
    ]
    
    ROOT.gROOT.ProcessLine('.L {output}/{id}/{id}.C+'.format(output=args.output,id=args.idName))
    
    variations = [
        ROOT.STAT_UP,
        ROOT.STAT_DOWN,
        ROOT.SYST_ALT_TEMPL,
        ROOT.SYST_TAG_PT30,
        ROOT.SYST_CMSSHAPE
        ]
    
    if args.trig:
        if 'DZ' in args.idName:
            eff = lambda pt, eta, var : getattr(ROOT, args.idName)(pt, eta, True, True, True, var)
        elif 'Mu8Leg' in args.idName or 'Ele12Leg' in args.idName:
            eff = lambda pt, eta, var : getattr(ROOT, args.idName)(pt, eta, True, True, var)
        else:
            eff = lambda pt, eta, var : getattr(ROOT, args.idName)(pt, eta, True, var)
    else:
        if args.object=='photon':
            eff = lambda pt, eta, var : getattr(ROOT, args.idName)(pt, eta, -1.0 if 'Fail' in args.idName else 1., True, False, var) # args: pt, eta, mva, pre, el-veto
        elif 'Iso' in args.idName:
            eff = lambda pt, eta, var : getattr(ROOT, args.idName)(pt, eta, True, True if 'Hpp' in args.idName else 0., var)
        else:
            eff = lambda pt, eta, var : getattr(ROOT, args.idName)(pt, eta, True, var)
    
    effCentral = lambda pt, eta : eff(pt, eta, ROOT.CENTRAL)
    effMax = lambda pt, eta : max(map(lambda v : eff(pt,eta,v), variations))-effCentral(pt,eta)
    effMin = lambda pt, eta : effCentral(pt,eta)-min(map(lambda v : eff(pt,eta,v), variations))
    
    effData = lambda pt, eta : eff(pt, eta, ROOT.EFF_DATA)
    effDataErr = lambda pt, eta : eff(pt, eta, ROOT.EFF_DATA_ERRSYM)
    
    effMC = lambda pt, eta : eff(pt, eta, ROOT.EFF_MC)
    effMCErr = lambda pt, eta : eff(pt, eta, ROOT.EFF_MC_ERRSYM)
    
    # pt bins
    xbins = array.array('d', [0.5*sum(ptbins[i:i+2]) for i in range(len(ptbins)-1)])
    xlo  = lambda bins: array.array('d', map(lambda (a,b): a-b, zip(bins, ptbins)))
    xhi  = lambda bins: array.array('d', map(lambda (a,b): a-b, zip(ptbins[1:], bins)))
    def y_eta(eta)   : return array.array('d', map(lambda pt : effCentral(pt, eta), xbins))
    def eyl_eta(eta) : return array.array('d', map(lambda pt : effMin(pt, eta), xbins))
    def eyh_eta(eta) : return array.array('d', map(lambda pt : effMax(pt, eta), xbins))

    def yData_eta(eta)   : return array.array('d', map(lambda pt : effData(pt, eta), xbins))
    def eylData_eta(eta) : return array.array('d', map(lambda pt : effData(pt, eta)+effDataErr(pt, eta), xbins))
    def eyhData_eta(eta) : return array.array('d', map(lambda pt : effData(pt, eta)-effDataErr(pt, eta), xbins))

    def yMC_eta(eta)   : return array.array('d', map(lambda pt : effMC(pt, eta), xbins))
    def eylMC_eta(eta) : return array.array('d', map(lambda pt : effMC(pt, eta)+effMCErr(pt, eta), xbins))
    def eyhMC_eta(eta) : return array.array('d', map(lambda pt : effMC(pt, eta)-effMCErr(pt, eta), xbins))

    canvas = ROOT.TCanvas()
    mg = ROOT.TMultiGraph('alletaBins', ';Probe p_{T};Scale Factor')
    mgData = ROOT.TMultiGraph('alletaBinsData', ';Probe p_{T};Data Efficiency')
    mgMC = ROOT.TMultiGraph('alletaBinsMC', ';Probe p_{T};MC Efficiency')
    
    for i in range(len(etabins)-1) :
        eta = .5*sum(etabins[i:i+2])
        bins2 = array.array('d', [b for b in xbins])
        graph = ROOT.TGraphAsymmErrors(len(xbins), bins2, y_eta(eta), xlo(bins2), xhi(bins2), eyl_eta(eta), eyh_eta(eta))
        graph.SetName('eff_etaBin%d'%i)
        graph.SetTitle('%.1f #leq #eta #leq %.1f' % tuple(etabins[i:i+2]))
        graph.SetMarkerColor(colors[i])
        graph.SetLineColor(colors[i])
        mg.Add(graph, 'p')
    
        graphData = ROOT.TGraphAsymmErrors(len(xbins), bins2, yData_eta(eta), xlo(bins2), xhi(bins2), eylData_eta(eta), eyhData_eta(eta))
        graphData.SetName('effData_etaBin%d'%i)
        graphData.SetTitle('%.1f #leq #eta #leq %.1f' % tuple(etabins[i:i+2]))
        graphData.SetMarkerColor(colors[i])
        graphData.SetLineColor(colors[i])
        mgData.Add(graphData, 'p')
    
        graphMC = ROOT.TGraphAsymmErrors(len(xbins), bins2, yMC_eta(eta), xlo(bins2), xhi(bins2), eylMC_eta(eta), eyhMC_eta(eta))
        graphMC.SetName('effMC_etaBin%d'%i)
        graphMC.SetTitle('%.1f #leq #eta #leq %.1f' % tuple(etabins[i:i+2]))
        graphMC.SetMarkerColor(colors[i])
        graphMC.SetLineColor(colors[i])
        mgMC.Add(graphMC, 'p')
    
    mg.SetMinimum(0.8)
    mg.SetMaximum(1.2)
    mg.Draw('a')
    leg = canvas.BuildLegend(.5,.2,.9,.4)
    for entry in leg.GetListOfPrimitives() :
        entry.SetOption('lp')
    leg.SetHeader(args.idNameNice)
    
    canvas.Print('{0}/{1}/scaleFactor_vs_pt.png'.format(args.output,args.idName))
    canvas.Print('{0}/{1}/scaleFactor_vs_pt.pdf'.format(args.output,args.idName))
    canvas.Print('{0}/{1}/scaleFactor_vs_pt.root'.format(args.output,args.idName))

    mgData.SetMinimum(0.0)
    mgData.SetMaximum(1.2)
    mgData.Draw('a')
    leg = canvas.BuildLegend(.5,.2,.9,.4)
    for entry in leg.GetListOfPrimitives() :
        entry.SetOption('lp')
    leg.SetHeader(args.idNameNice)
    
    canvas.Print('{0}/{1}/effData_vs_pt.png'.format(args.output,args.idName))
    canvas.Print('{0}/{1}/effData_vs_pt.pdf'.format(args.output,args.idName))
    canvas.Print('{0}/{1}/effData_vs_pt.root'.format(args.output,args.idName))

    mgMC.SetMinimum(0.0)
    mgMC.SetMaximum(1.2)
    mgMC.Draw('a')
    leg = canvas.BuildLegend(.5,.2,.9,.4)
    for entry in leg.GetListOfPrimitives() :
        entry.SetOption('lp')
    leg.SetHeader(args.idNameNice)
    
    canvas.Print('{0}/{1}/effMC_vs_pt.png'.format(args.output,args.idName))
    canvas.Print('{0}/{1}/effMC_vs_pt.pdf'.format(args.output,args.idName))
    canvas.Print('{0}/{1}/effMC_vs_pt.root'.format(args.output,args.idName))

    # eta bins
    xbins = array.array('d', [0.5*sum(etabins[i:i+2]) for i in range(len(etabins)-1)])
    xlo  = lambda bins: array.array('d', map(lambda (a,b): a-b, zip(bins, etabins)))
    xhi  = lambda bins: array.array('d', map(lambda (a,b): a-b, zip(etabins[1:], bins)))
    def y_pt(pt) : return array.array('d', map(lambda eta : effCentral(pt, eta), xbins))
    def eyl_pt(pt) : return array.array('d', map(lambda eta : effMin(pt, eta), xbins))
    def eyh_pt(pt) : return array.array('d', map(lambda eta : effMax(pt, eta), xbins))
    
    def yData_pt(pt)   : return array.array('d', map(lambda eta : effData(pt, eta), xbins))
    def eylData_pt(pt) : return array.array('d', map(lambda eta : effData(pt, eta)+effDataErr(pt, eta), xbins))
    def eyhData_pt(pt) : return array.array('d', map(lambda eta : effData(pt, eta)-effDataErr(pt, eta), xbins))

    def yMC_pt(pt)   : return array.array('d', map(lambda eta : effMC(pt, eta), xbins))
    def eylMC_pt(pt) : return array.array('d', map(lambda eta : effMC(pt, eta)+effMCErr(pt, eta), xbins))
    def eyhMC_pt(pt) : return array.array('d', map(lambda eta : effMC(pt, eta)-effMCErr(pt, eta), xbins))

    canvas = ROOT.TCanvas()
    mg = ROOT.TMultiGraph('allptBins', ';Probe #eta;Scale Factor')
    mgData = ROOT.TMultiGraph('allptBinsData', ';Probe #eta;Data Efficiency')
    mgMC = ROOT.TMultiGraph('allptBinsMC', ';Probe #eta;MC Efficiency')
    
    for i in range(len(ptbins)-1) :
        pt = .5*sum(ptbins[i:i+2])
        #bins2 = array.array('d', [b-1.5+i for b in xbins])
        bins2 = array.array('d', [b for b in xbins])
        graph = ROOT.TGraphAsymmErrors(len(xbins), bins2, y_pt(pt), xlo(bins2), xhi(bins2), eyl_pt(pt), eyh_pt(pt))
        graph.SetName('eff_ptBin%d'%i)
        graph.SetTitle('%.1f #leq p_{T} #leq %.1f' % tuple(ptbins[i:i+2]))
        graph.SetMarkerColor(colors[i])
        graph.SetLineColor(colors[i])
        mg.Add(graph, 'p')
    
        graphData = ROOT.TGraphAsymmErrors(len(xbins), bins2, yData_eta(pt), xlo(bins2), xhi(bins2), eylData_eta(pt), eyhData_eta(pt))
        graphData.SetName('effData_ptBin%d'%i)
        graphData.SetTitle('%.1f #leq p_{T} #leq %.1f' % tuple(ptbins[i:i+2]))
        graphData.SetMarkerColor(colors[i])
        graphData.SetLineColor(colors[i])
        mgData.Add(graphData, 'p')
    
        graphMC = ROOT.TGraphAsymmErrors(len(xbins), bins2, yMC_eta(pt), xlo(bins2), xhi(bins2), eylMC_eta(pt), eyhMC_eta(pt))
        graphMC.SetName('effMC_ptBin%d'%i)
        graphMC.SetTitle('%.1f #leq p_{T} #leq %.1f' % tuple(ptbins[i:i+2]))
        graphMC.SetMarkerColor(colors[i])
        graphMC.SetLineColor(colors[i])
        mgMC.Add(graphMC, 'p')

    mg.SetMinimum(0.8)
    mg.SetMaximum(1.2)
    mg.Draw('a')
    leg = canvas.BuildLegend(.5,.2,.9,.4)
    for entry in leg.GetListOfPrimitives() :
        entry.SetOption('lp')
    leg.SetHeader(args.idNameNice)
    
    canvas.Print('{0}/{1}/scaleFactor_vs_eta.png'.format(args.output,args.idName))
    canvas.Print('{0}/{1}/scaleFactor_vs_eta.pdf'.format(args.output,args.idName))
    canvas.Print('{0}/{1}/scaleFactor_vs_eta.root'.format(args.output,args.idName))

    mgData.SetMinimum(0.0)
    mgData.SetMaximum(1.2)
    mgData.Draw('a')
    leg = canvas.BuildLegend(.5,.2,.9,.4)
    for entry in leg.GetListOfPrimitives() :
        entry.SetOption('lp')
    leg.SetHeader(args.idNameNice)
    
    canvas.Print('{0}/{1}/effData_vs_eta.png'.format(args.output,args.idName))
    canvas.Print('{0}/{1}/effData_vs_eta.pdf'.format(args.output,args.idName))
    canvas.Print('{0}/{1}/effData_vs_eta.root'.format(args.output,args.idName))

    mgMC.SetMinimum(0.0)
    mgMC.SetMaximum(1.2)
    mgMC.Draw('a')
    leg = canvas.BuildLegend(.5,.2,.9,.4)
    for entry in leg.GetListOfPrimitives() :
        entry.SetOption('lp')
    leg.SetHeader(args.idNameNice)
    
    canvas.Print('{0}/{1}/effMC_vs_eta.png'.format(args.output,args.idName))
    canvas.Print('{0}/{1}/effMC_vs_eta.pdf'.format(args.output,args.idName))
    canvas.Print('{0}/{1}/effMC_vs_eta.root'.format(args.output,args.idName))

    # pt eta bins
    tfile = ROOT.TFile('{0}/{1}_scalefactor.root'.format(args.output,args.idName),'recreate')

    pteta = ROOT.TH2F('scalefactor','scalefactor;Probe p_{T};Probe #eta',len(ptbins)-1,array.array('d',ptbins),len(etabins)-1,array.array('d',etabins))
    for p in range(len(ptbins)-1):
        pt = .5*sum(ptbins[p:p+2])
        for e in range(len(etabins)-1):
            eta = .5*sum(etabins[e:e+2])
            sf = effCentral(pt,eta)
            hi = effMax(pt,eta)
            lo = effMin(pt,eta)
            err = abs(hi-sf + sf-lo)/2.
            #print pt,eta,sf,hi,lo,err
            pteta.SetBinContent(pteta.FindBin(pt,eta),sf)
            pteta.SetBinError(pteta.FindBin(pt,eta),err)

    save2D(pteta.Clone(),'scaleFactor','{0}/{1}'.format(args.output,args.idName))
    pteta.Write()

    ptetaData = ROOT.TH2F('effData','effData;Probe p_{T};Probe #eta',len(ptbins)-1,array.array('d',ptbins),len(etabins)-1,array.array('d',etabins))
    for p in range(len(ptbins)-1):
        pt = .5*sum(ptbins[p:p+2])
        for e in range(len(etabins)-1):
            eta = .5*sum(etabins[e:e+2])
            sf = effData(pt,eta)
            err = effDataErr(pt,eta)
            #print pt,eta,sf,hi,lo,err
            ptetaData.SetBinContent(ptetaData.FindBin(pt,eta),sf)
            ptetaData.SetBinError(ptetaData.FindBin(pt,eta),err)

    save2D(ptetaData.Clone(),'effData','{0}/{1}'.format(args.output,args.idName))
    ptetaData.Write()
    
    ptetaMC = ROOT.TH2F('effMC','effMC;Probe p_{T};Probe #eta',len(ptbins)-1,array.array('d',ptbins),len(etabins)-1,array.array('d',etabins))
    for p in range(len(ptbins)-1):
        pt = .5*sum(ptbins[p:p+2])
        for e in range(len(etabins)-1):
            eta = .5*sum(etabins[e:e+2])
            sf = effMC(pt,eta)
            err = effMCErr(pt,eta)
            #print pt,eta,sf,hi,lo,err
            ptetaMC.SetBinContent(ptetaMC.FindBin(pt,eta),sf)
            ptetaMC.SetBinError(ptetaMC.FindBin(pt,eta),err)

    save2D(ptetaMC.Clone(),'effMC','{0}/{1}'.format(args.output,args.idName))
    ptetaMC.Write()
    
    # ------- Latex
    def formatValue(var, varerr, etabin, ptbin) :
        eta = etabins[etabin-1]+.001
        pt = ptbins[ptbin-1]+.001
        value = eff(pt, eta, var)
        if var is ROOT.CENTRAL :
            return '$%1.4f^{+%1.4f}_{-%1.4f}$' % (value, effMax(pt,eta), effMin(pt,eta))
        else :
            err = eff(ptbins[ptbin-1]+.001, etabins[etabin-1]+.001, varerr)
            return '$%1.4f \\pm %1.4f$' % (value, err)
    
    output = ''
    
    neta = len(etabins)-1
    output += '''\\begin{table}[htbp]
    \\centering
    \\begin{tabular}{%s}
    \hline
    ''' % 'c'.join(['|']*(neta+3))
    
    etaLabels = ['$p_T$', '-']
    for etabin in xrange(1, neta+1) :
        etalo = etabins[etabin-1]
        etahi = etabins[etabin]
        etaLabels.append('$%1.1f < |\\eta| < %1.1f$' % (etalo, etahi))
    
    output += '        ' + ' & '.join(etaLabels) + '\\\\ \hline\n'
    
    npt = len(ptbins)-1
    for ptbin in xrange(1, npt+1) :
        ptlo = ptbins[ptbin-1]
        pthi = ptbins[ptbin]
        ptLabel = '$%3.0f - %3.0f$' % (ptlo, pthi)
    
        output += '      \\multirow{3}{*}{%s} \n' % ptLabel
    
        dataLine, mcLine, ratioLine = [['', name] for name in ['Data', 'MC', 'Ratio']]
        for etabin in xrange(1, neta+1) :
            dataLine.append(formatValue(ROOT.EFF_DATA, ROOT.EFF_DATA_ERRSYM, etabin, ptbin))
            mcLine.append(formatValue(ROOT.EFF_MC, ROOT.EFF_MC_ERRSYM, etabin, ptbin))
            ratioLine.append(formatValue(ROOT.CENTRAL, None, etabin, ptbin))
    
    
        output += '        %s \\\\ \n' % ' & '.join(dataLine)
        output += '        %s \\\\ \n' % ' & '.join(mcLine)
        output += '        %s \\\\ \\hline\n' % ' & '.join(ratioLine)
    
    output += '''   \\end{tabular}
    \\caption{Efficiency table for %s}
    \\end{table}
    %% Generated with DevTools/TagAndProbe version %s
    ''' % (args.idName, __gitversion__)
    with open('{0}/{1}/table.tex'.format(args.output,args.idName), 'w') as fout :
        fout.write(output)


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='TagAndProbe Plotter')

    parser.add_argument('object', type=str, choices=['electron','muon','photon'], help='Physics object')
    parser.add_argument('output', type=str, help='Directory for output')
    parser.add_argument('idName', type=str, help='Name of ID')
    parser.add_argument('idNameNice', type=str, help='Name of ID for printing')
    parser.add_argument('--trig', action='store_true', help='Is a trigger')

    return parser.parse_args(argv)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    plot(args)

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)

