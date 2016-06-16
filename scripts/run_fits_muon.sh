python DevTools/TagAndProbe/python/fitter.py muon -mc tagAndProbe/muon/dy_nlo.root -mcLO tagAndProbe/muon/dy_lo.root -data tagAndProbe/muon/data.root
python DevTools/TagAndProbe/python/dumpTagProbeTreeHTML.py --data fits_muon.root -i muonFits -o fits_muon/
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon MediumID "Medium ID"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon TightID "Tight ID"
