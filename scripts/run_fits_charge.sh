python DevTools/TagAndProbe/python/chargeFitter.py -mc /hdfs/store/user/dntaylor/tagAndProbe/charge/mcNLO.root -mcLO /hdfs/store/user/dntaylor/tagAndProbe/charge/mcLO.root -data /hdfs/store/user/dntaylor/tagAndProbe/charge/data.root
python DevTools/TagAndProbe/python/dumpTagProbeTreeHTML.py --data fits_charge.root -i ChargeFits -o fits_charge/
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py charge fits_charge/ OppositeSign "Opposite Sign"
python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py charge fits_charge/ SameSign "Same Sign"
