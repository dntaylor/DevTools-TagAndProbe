baseDir=/hdfs/store/user/dntaylor/tagAndProbe/photon
python DevTools/TagAndProbe/python/photonFitter.py --mcFileName $baseDir/dy_nlo.root --mcLOFileName $baseDir/dy_lo.root --dataFileName $baseDir/data.root
python DevTools/TagAndProbe/python/dumpTagProbeTreeHTML.py --data fits_photon.root -i PhotonFits -o fits_photon
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py photon fits_photon Preselection "Photon Preselection"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py photon fits_photon MVA0p0Pre "MVA>0.0 w/ Preselection"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py photon fits_photon MVA0p0 "MVA>0.0"
#python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py photon fits_photon MVA0p0PreFail "MVA<0.0 w/ Preselection"
#python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py photon fits_photon MVA0p0Fail "MVA<0.0"
