import os
import sys
import glob

import ROOT

import DevTools.Plotter.CMS_lumi as CMS_lumi
import DevTools.Plotter.tdrstyle as tdrstyle
from DevTools.TagAndProbe.utilities import getBinning
from DevTools.Plotter.utilities import getLumi
from DevTools.Utilities.utilities import python_mkdir

ROOT.gROOT.SetBatch(ROOT.kTRUE)
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
    eff.GetXaxis().SetTitleOffset(1.5)
    eff.GetYaxis().SetTitleOffset(1.5)
    eff.GetZaxis().SetTitleOffset(1.5)
    eff.Draw('colz')
    setStyle(canvas,position=0)
    for t in ['png','pdf']:
        name = '{0}/{2}/{1}.{2}'.format(outdir,savename,t)
        python_mkdir(os.path.dirname(name))
        canvas.Print(name)

def save1D(eff,var,savename,outdir):
    colors = [ROOT.kRed, ROOT.kGreen, ROOT.kBlue, ROOT.kBlack, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kGreen+2, ROOT.kRed-3, ROOT.kCyan+1, ROOT.kMagenta-3, ROOT.kViolet-1, ROOT.kSpring+10]
    canvas = ROOT.TCanvas(savename,savename,50,50,600,600)
    canvas.SetLogx(True if var=='pt' else False)
    canvas.SetRightMargin(0.04)
    nbins = eff.GetNbinsX() if var=='eta' else eff.GetNbinsY()

    for b in range(nbins):
       args = ['eff_{0}_{1}'.format(var,b),b+1,b+1]
       proj = eff.ProjectionX(*args) if var=='pt' else eff.ProjectionY(*args)
       #proj.SetLineColor(colors[b])
       if b==0:
           proj.GetXaxis().SetTitleOffset(1.5)
           proj.GetYaxis().SetTitleOffset(1.5)
           proj.GetYaxis().SetTitle(eff.GetZaxis().GetTitle())
           proj.GetYaxis().SetTitleSize(0.04)
           proj.Draw()
           proj.SetMinimum(0.)
           proj.SetMaximum(1.)
       else:
           proj.Draw('same')

    setStyle(canvas,position=0)
    for t in ['png','pdf']:
        name = '{0}/{2}/{1}.{2}'.format(outdir,savename,t)
        python_mkdir(os.path.dirname(name))
        canvas.Print(name)

def saveEff(obj,trigname,outfile,outdir):
    effdir = {
        'electron' : 'GsfElectronToTrigger',
        'photon' : 'GsfElectronToTrigger',
        'muon'     : 'muonEffs',
    }
    
    ptvar = {
        'electron' : 'probe_Ele_pt',
        'photon' : 'probe_Ele_pt',
        'muon'     : 'probe_pt',
    }
    
    etavar = {
        'electron' : 'probe_Ele_eta',
        'photon' : 'probe_Ele_eta',
        'muon'     : 'probe_eta',
    }
    
    # get efficiency from file
    fname = 'efficiency-data-{0}.root'.format(trigname)
    tfile = ROOT.TFile(fname)
    #effname = '{0}/{1}/cnt_eff_plots/{2}_{3}_PLOT'.format(effdir[obj],trigname,ptvar[obj],etavar[obj])
    directory = tfile.Get(effdir[obj]).GetListOfKeys()[0].GetName()
    effname = [x.GetName() for x in tfile.Get('{0}/{1}/cnt_eff_plots'.format(effdir[obj],directory)).GetListOfKeys() if x.GetName().startswith('{0}_{1}'.format(ptvar[obj],etavar[obj]))][0]
    canvas = tfile.Get('{0}/{1}/cnt_eff_plots/{2}'.format(effdir[obj],directory,effname))
    eff = canvas.GetPrimitive(effname)

    mcfname = 'efficiency-mc-{0}.root'.format(trigname)
    mctfile = ROOT.TFile(mcfname)
    mcdirectory = mctfile.Get(effdir[obj]).GetListOfKeys()[0].GetName()
    mceffname = [x.GetName() for x in mctfile.Get('{0}/{1}/cnt_eff_plots'.format(effdir[obj],mcdirectory)).GetListOfKeys() if x.GetName().startswith('{0}_{1}'.format(ptvar[obj],etavar[obj]))][0]
    mccanvas = mctfile.Get('{0}/{1}/cnt_eff_plots/{2}'.format(effdir[obj],mcdirectory,mceffname))
    mceff = mccanvas.GetPrimitive(mceffname)

    outfile.cd()
    outfile.mkdir(trigname).cd()
    eff.Write()
    mceff.Write()
    
    # plot 2d and 1d efficiencies
    ptbins = getBinning(obj,'pt',trig=True)
    etabins = getBinning(obj,'eta',trig=True)
    
    save2D(eff.Clone(),trigname,outdir)
    save1D(eff,'pt',trigname+'_pt',outdir)
    save1D(eff,'eta',trigname+'_eta',outdir)
    save2D(mceff.Clone(),trigname+'_mc',outdir)
    save1D(mceff,'pt',trigname+'_mc_pt',outdir)
    save1D(mceff,'eta',trigname+'_mc_eta',outdir)

def doAllEfficiencies(obj):
    outdir = 'fits_{0}_trigger'.format(obj)
    python_mkdir(outdir)
    outfile = ROOT.TFile('{0}/efficiencies_{1}_trigger.root'.format(outdir,obj),'recreate')

    effList = {
        'electron': [
            'passingHLTEle23Ele12DZ',
            'passingHLTEle23Ele12Leg1',
            'passingHLTEle23Ele12Leg2',
            'passingHLTEle27Tight',
            'passingHLTEle45',
        ],
        'muon': [
            'passingIsoMu24',
            'passingIsoMu24ORIsoTkMu24',
            'passingIsoTkMu24',
            'passingMu17Leg',
            'passingMu8Leg',
            'passingMu8LegDZ',
            'passingMu8ORTkMu8Leg',
            'passingMu8ORTkMu8LegDZ',
            'passingTkMu8Leg',
            'passingTkMu8LegDZ',
            'passingMu50',
        ],
        'photon': [
            'passingPho18Leg',
            'passingPho30Leg',
            'passingPho60Leg',
            'passingPho175',
        ],
    }

    for eff in effList[obj]:
        print 'processing {0} {1}'.format(obj,eff)
        saveEff(obj,eff,outfile,outdir)



def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    #doAllEfficiencies('muon')
    #doAllEfficiencies('electron')
    doAllEfficiencies('photon')


if __name__ == "__main__":
    status = main()
    sys.exit(status)




