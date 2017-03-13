obj=muon
for id in IsoMu24 IsoTkMu24 IsoMu24ORIsoTkMu24 Mu50 IsoMu24ORIsoTkMu24ORMu50 Mu17Leg Mu8Leg TkMu8Leg Mu8ORTkMu8Leg Mu8LegDZ TkMu8LegDZ Mu8ORTkMu8LegDZ; do
    baseDir=/hdfs/store/user/dntaylor/tagAndProbe/$obj
    nohup python DevTools/TagAndProbe/python/fitter.py $obj -mc $baseDir/dy_nlo.root -mcLO $baseDir/dy_lo.root -data $baseDir/data.root --fitsToRun $id --outputFileName fits_"$id".root > fits_"$id".log 2>&1 &
done

obj=electron
for id in Ele27WPTight Ele27Eta2p1WPTight Ele23Leg Ele12Leg Ele12LegDZ; do
    baseDir=/hdfs/store/user/dntaylor/tagAndProbe/$obj
    nohup python DevTools/TagAndProbe/python/fitter.py $obj -mc $baseDir/dy_nlo.root -mcLO $baseDir/dy_lo.root -data $baseDir/data.root --fitsToRun $id --outputFileName fits_"$id".root > fits_"$id".log 2>&1 &
done
