obj=photon
baseDir=/hdfs/store/user/dntaylor/tagAndProbe/$obj
for id in Pho30Leg Pho18Leg Pho60Leg Pho175; do
    nohup python DevTools/TagAndProbe/python/fitter.py $obj -mc $baseDir/dy_nlo.root -mcLO $baseDir/dy_lo.root -data $baseDir/data.root --fitsToRun $id --outputFileName fits_"$id".root > fits_"$id".log 2>&1 &
done

#python DevTools/TagAndProbe/python/fitter.py $obj -mc $baseDir/dy_nlo.root -mcLO $baseDir/dy_lo.root -data $baseDir/data.root --fitsToRun $id --outputFileName fits_photon_trigger.root > fits_photon_trigger.log 2>&1
