#!/usr/bin/env python

import os
import sys
import glob
import ROOT

def combineScaleFactors():
    inputDir1 = 'fits_charge_z1'
    inputDir2 = 'fits_charge_z2'
    outputFileName = 'scalefactors_charge.root'
    outfile = ROOT.TFile(outputFileName,'recreate')
    for f1, f2 in zip(sorted(glob.glob('{0}/*_scalefactor.root'.format(inputDir1))), sorted(glob.glob('{0}/*_scalefactor.root'.format(inputDir2)))):
        idname = os.path.basename(f1).split('_')[0]
        idfile1 = ROOT.TFile(f1)
        sf = idfile1.Get('scalefactor').Clone()
        em = idfile1.Get('effMC').Clone()
        ed = idfile1.Get('effData').Clone()
        idfile2 = ROOT.TFile(f2)
        sf2 = idfile1.Get('scalefactor').Clone()
        em2 = idfile2.Get('effMC').Clone()
        ed2 = idfile2.Get('effData').Clone()
        outfile.cd()
        for p in range(sf.GetNbinsX()):
            for e in range(sf.GetNbinsY()):
                pt = sf.GetXaxis().GetBinCenter(p+1)
                eta = sf.GetYaxis().GetBinCenter(e+1)
                if pt<20:
                    b = sf.FindBin(pt,eta)
                    sf.SetBinContent(b,sf2.GetBinContent(b))
                    sf.SetBinError(b,sf2.GetBinError(b))
                    b = em.FindBin(pt,eta)
                    em.SetBinContent(b,em2.GetBinContent(b))
                    em.SetBinError(b,em2.GetBinError(b))
                    b = ed.FindBin(pt,eta)
                    ed.SetBinContent(b,ed2.GetBinContent(b))
                    ed.SetBinError(b,ed2.GetBinError(b))
        sf.SetName(idname)
        sf.SetTitle(idname)
        sf.Write()
        em.SetName(idname+'effMC')
        em.SetTitle(idname+'effMC')
        em.Write()
        ed.SetName(idname+'effData')
        ed.SetTitle(idname+'effData')
        ed.Write()


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    combineScaleFactors()


if __name__ == "__main__":
    status = main()
    sys.exit(status)

