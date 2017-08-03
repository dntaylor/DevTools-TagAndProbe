baseDir=/nfs_scratch/dntaylor/data/photon
python photonFitter.py --mcFileName $baseDir/mcNLO.root --mcLOFileName $baseDir/mcLO.root --dataFileName $baseDir/data.root
python DevTools/TagAndProbe/python/dumpTagProbeTreeHTML.py --data fits_photon.root -i PhotonFits -o fits_photon
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py photon fits_photon Preselection "Photon Preselection"
