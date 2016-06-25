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


def plot(obj,idName,idNameNice):

    ptbins = getBinning(obj,'pt')
    etabins = getBinning(obj,'eta')
    colors = [ROOT.kRed, ROOT.kGreen, ROOT.kBlue, ROOT.kBlack, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kGreen+2, ROOT.kRed-3, ROOT.kCyan+1, ROOT.kMagenta-3, ROOT.kViolet-1, ROOT.kSpring+10]
    
    ROOT.gROOT.ProcessLine('.L fits_{obj}/{id}/{id}.C+'.format(obj=obj,id=idName))
    
    variations = [
        ROOT.STAT_UP,
        ROOT.STAT_DOWN,
        ROOT.SYST_ALT_TEMPL,
        ROOT.SYST_TAG_PT30,
        ROOT.SYST_CMSSHAPE
        ]
    
    if 'Iso' in idName:
        eff = lambda pt, eta, var : getattr(ROOT, idName)(pt, eta, True, 0., var)
    else:
        eff = lambda pt, eta, var : getattr(ROOT, idName)(pt, eta, True, var)
    
    effCentral = lambda pt, eta : eff(pt, eta, ROOT.CENTRAL)
    effMax = lambda pt, eta : max(map(lambda v : eff(pt,eta,v), variations))-effCentral(pt,eta)
    effMin = lambda pt, eta : effCentral(pt,eta)-min(map(lambda v : eff(pt,eta,v), variations))
    
    
    # pt bins
    xbins = array.array('d', [0.5*sum(ptbins[i:i+2]) for i in range(len(ptbins)-1)])
    xlo  = lambda bins: array.array('d', map(lambda (a,b): a-b, zip(bins, ptbins)))
    xhi  = lambda bins: array.array('d', map(lambda (a,b): a-b, zip(ptbins[1:], bins)))
    def y_eta(eta) : return array.array('d', map(lambda pt : effCentral(pt, eta), xbins))
    def eyl_eta(eta) : return array.array('d', map(lambda pt : effMin(pt, eta), xbins))
    def eyh_eta(eta) : return array.array('d', map(lambda pt : effMax(pt, eta), xbins))

    canvas = ROOT.TCanvas()
    mg = ROOT.TMultiGraph('alletaBins', ';Probe p_{T};Scale Factor')
    
    for i in range(len(etabins)-1) :
        eta = .5*sum(etabins[i:i+2])
        #bins2 = array.array('d', [b-1.5+i for b in xbins])
        bins2 = array.array('d', [b for b in xbins])
        graph = ROOT.TGraphAsymmErrors(len(xbins), bins2, y_eta(eta), xlo(bins2), xhi(bins2), eyl_eta(eta), eyh_eta(eta))
        graph.SetName('eff_etaBin%d'%i)
        graph.SetTitle('%.1f #leq #eta #leq %.1f' % tuple(etabins[i:i+2]))
        graph.SetMarkerColor(colors[i])
        graph.SetLineColor(colors[i])
        mg.Add(graph, 'p')
    
    mg.SetMinimum(0.8)
    mg.SetMaximum(1.05)
    mg.Draw('a')
    leg = canvas.BuildLegend(.5,.2,.9,.4)
    for entry in leg.GetListOfPrimitives() :
        entry.SetOption('lp')
    leg.SetHeader(idNameNice)
    
    canvas.Print('fits_{0}/{1}/scaleFactor_vs_pt.png'.format(obj,idName))
    canvas.Print('fits_{0}/{1}/scaleFactor_vs_pt.pdf'.format(obj,idName))
    canvas.Print('fits_{0}/{1}/scaleFactor_vs_pt.root'.format(obj,idName))

    # eta bins
    xbins = array.array('d', [0.5*sum(etabins[i:i+2]) for i in range(len(etabins)-1)])
    xlo  = lambda bins: array.array('d', map(lambda (a,b): a-b, zip(bins, etabins)))
    xhi  = lambda bins: array.array('d', map(lambda (a,b): a-b, zip(etabins[1:], bins)))
    def y_pt(pt) : return array.array('d', map(lambda eta : effCentral(pt, eta), xbins))
    def eyl_pt(pt) : return array.array('d', map(lambda eta : effMin(pt, eta), xbins))
    def eyh_pt(pt) : return array.array('d', map(lambda eta : effMax(pt, eta), xbins))
    
    canvas = ROOT.TCanvas()
    mg = ROOT.TMultiGraph('allptBins', ';Probe #eta;Scale Factor')
    
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
    
    mg.SetMinimum(0.8)
    mg.SetMaximum(1.05)
    mg.Draw('a')
    leg = canvas.BuildLegend(.5,.2,.9,.4)
    for entry in leg.GetListOfPrimitives() :
        entry.SetOption('lp')
    leg.SetHeader(idNameNice)
    
    canvas.Print('fits_{0}/{1}/scaleFactor_vs_eta.png'.format(obj,idName))
    canvas.Print('fits_{0}/{1}/scaleFactor_vs_eta.pdf'.format(obj,idName))
    canvas.Print('fits_{0}/{1}/scaleFactor_vs_eta.root'.format(obj,idName))

    # pt eta bins
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

    save2D(pteta.Clone(),'scaleFactor','fits_{0}/{1}'.format(obj,idName))
            
    
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
    ''' % (idName, __gitversion__)
    with open('fits_{0}/{1}/table.tex'.format(obj,idName), 'w') as fout :
        fout.write(output)


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='TagAndProbe Plotter')

    parser.add_argument('object', type=str, choices=['electron','muon'], help='Physics object')
    parser.add_argument('idName', type=str, help='Name of ID')
    parser.add_argument('idNameNice', type=str, help='Name of ID for printing')

    return parser.parse_args(argv)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    plot(args.object,args.idName,args.idNameNice)

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)

