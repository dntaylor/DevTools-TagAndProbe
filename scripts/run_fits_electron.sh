python DevTools/TagAndProbe/python/fitter.py electron -mc tagAndProbe/electron/dy_nlo.root -mcLO tagAndProbe/electron/dy_lo.root -data tagAndProbe/electron/data.root
python DevTools/TagAndProbe/python/dumpTagProbeTreeHTML.py --data fits_electron.root -i electronFits -o fits_electron/
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py electron CutBasedIDVeto "Cut Based ID Veto"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py electron CutBasedIDLoose "Cut Based ID Loose"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py electron CutBasedIDMedium "Cut Based ID Medium"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py electron CutBasedIDTight "Cut Based ID Tight"
