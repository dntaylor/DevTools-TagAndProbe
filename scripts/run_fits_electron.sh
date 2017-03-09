baseDir=/hdfs/store/user/dntaylor/tagAndProbe/electron
python DevTools/TagAndProbe/python/fitter.py electron -mc $baseDir/dy_nlo.root -mcLO $baseDir/dy_lo.root -data $baseDir/data.root
#python DevTools/TagAndProbe/python/dumpTagProbeTreeHTML.py --data fits_electron.root -i electronFits -o fits_electron/
#python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py electron CutBasedIDVeto "Cut Based ID Veto"
#python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py electron CutBasedIDLoose "Cut Based ID Loose"
#python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py electron CutBasedIDMedium "Cut Based ID Medium"
#python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py electron CutBasedIDTight "Cut Based ID Tight"
