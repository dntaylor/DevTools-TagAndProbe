baseDir=/nfs_scratch/dntaylor/data/photon
#python DevTools/TagAndProbe/python/photonFitter.py --mcFileName $baseDir/mcNLO.root --mcLOFileName $baseDir/mcLO.root --dataFileName $baseDir/data.root
python DevTools/TagAndProbe/python/dumpTagProbeTreeHTML.py --data fits_photon.root -i PhotonFits -o fits_photon
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py photon fits_photon Preselection "Photon Preselection"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py photon fits_photon MVA0p0Pre "MVA>0.0 w/ Preselection"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py photon fits_photon MVA0p0 "MVA>0.0"
