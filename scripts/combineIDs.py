#!/usr/bin/env python

import os
import sys
import glob
import ROOT

def combineScaleFactors(obj):
    inputDir = 'fits_{0}'.format(obj)
    outputFileName = '{0}/scalefactors_{1}.root'.format(inputDir,obj)
    outfile = ROOT.TFile(outputFileName,'recreate')
    for fname in glob.glob('{0}/*_scalefactor.root'.format(inputDir)):
        idname = os.path.basename(fname).split('_')[0]
        idfile = ROOT.TFile(fname)
        sf = idfile.Get('scalefactor').Clone()
        em = idfile.Get('effMC').Clone()
        ed = idfile.Get('effMC').Clone()
        outfile.cd()
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

    combineScaleFactors('muon')
    combineScaleFactors('electron')
    combineScaleFactors('photon')


if __name__ == "__main__":
    status = main()
    sys.exit(status)

