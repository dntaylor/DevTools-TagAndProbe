#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
import array
import sys
import subprocess
import argparse

from DevTools.TagAndProbe.utilities import getBinning

__gitversion__ = subprocess.check_output(["git", "describe", "--always"]).strip()

def plot(args):
    ptbins = getBinning(args.object,'pt')
    etabins = getBinning(args.object,'eta')
    colors = [ROOT.kRed, ROOT.kGreen, ROOT.kBlue, ROOT.kBlack, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kGreen+2, ROOT.kRed-3, ROOT.kCyan+1, ROOT.kMagenta-3, ROOT.kViolet-1, ROOT.kSpring+10]
    
    ROOT.gROOT.ProcessLine('.L fits_{obj}/{id}/{id}.C+'.format(obj=args.object,id=args.idName))
    
    variations = [
        ROOT.STAT_UP,
        ROOT.STAT_DOWN,
        ROOT.SYST_ALT_TEMPL,
        ROOT.SYST_TAG_PT30,
        ROOT.SYST_CMSSHAPE
        ]
    
    if 'Iso' in args.idName:
        eff = lambda pt, eta, var : getattr(ROOT, args.idName)(pt, eta, True, 0., var)
    else:
        eff = lambda pt, eta, var : getattr(ROOT, args.idName)(pt, eta, True, var)
    
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
    leg.SetHeader(args.idNameNice)
    
    canvas.Print('fits_{0}/{1}/scaleFactor_vs_pt.png'.format(args.object,args.idName))
    canvas.Print('fits_{0}/{1}/scaleFactor_vs_pt.root'.format(args.object,args.idName))

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
    leg.SetHeader(args.idNameNice)
    
    canvas.Print('fits_{0}/{1}/scaleFactor_vs_eta.png'.format(args.object,args.idName))
    canvas.Print('fits_{0}/{1}/scaleFactor_vs_eta.root'.format(args.object,args.idName))
    
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
    %% Generated with DiBosonTP version %s
    ''' % (args.idName, __gitversion__)
    with open('fits_{0}/{1}/table.tex'.format(args.object,args.idName), 'w') as fout :
        fout.write(output)


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='TagAndProbe Plotter')

    parser.add_argument('object', type=str, choices=['electron','muon'], help='Physics object')
    parser.add_argument('idName', type=str, help='ID name to plot')
    parser.add_argument('idNameNice', type=str, help='Pretty name for ID')

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

