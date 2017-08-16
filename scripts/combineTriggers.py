#!/usr/bin/env python

import os
import sys
import glob
import ROOT

def combineScaleFactors(obj):
    inputDir = 'fits_{0}_trigger'.format(obj)
    outputFileName = '{0}/scalefactors_{1}_trigger.root'.format(inputDir,obj)
    outfile = ROOT.TFile(outputFileName,'recreate')
    for fname in glob.glob('{0}/*_scalefactor.root'.format(inputDir)):
        idname = os.path.basename(fname).split('_')[0]
        idfile = ROOT.TFile(fname)
        sf = idfile.Get('scalefactor').Clone()
        ed = idfile.Get('effData').Clone()
        em = idfile.Get('effMC').Clone()
        outfile.cd()
        sf.SetName(idname)
        sf.SetTitle(idname)
        sf.Write()
        ed.SetName(idname+'_effData')
        ed.SetTitle(idname+'_effData')
        ed.Write()
        em.SetName(idname+'_effMC')
        em.SetTitle(idname+'_effMC')
        em.Write()


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    combineScaleFactors('muon')
    combineScaleFactors('electron')


if __name__ == "__main__":
    status = main()
    sys.exit(status)

