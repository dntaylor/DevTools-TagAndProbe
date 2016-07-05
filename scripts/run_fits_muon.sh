python DevTools/TagAndProbe/python/fitter.py muon -mc tagAndProbe/muon/dy_nlo.root -mcLO tagAndProbe/muon/dy_lo.root -data tagAndProbe/muon/data.root
python DevTools/TagAndProbe/python/dumpTagProbeTreeHTML.py --data fits_muon.root -i muonFits -o fits_muon/
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon LooseID "Loose ID"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon LooseIsoFromLooseID "Loose Iso from Loose ID"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon MediumID "Medium ID"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon LooseIsoFromMediumID "Loose Iso from Medium ID"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon TightIsoFromMediumID "Tight Iso from Medium ID"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon MediumIDICHEP "Medium ID (ICHEP)"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon LooseIsoFromMediumIDICHEP "Loose Iso from Medium ID (ICHEP)"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon TightIsoFromMediumIDICHEP "Tight Iso from Medium ID (ICHEP)"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon TightID "Tight ID"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon TightIsoFromTightID "Tight Iso from Tight ID"
