#python DevTools/TagAndProbe/python/dumpTagProbeTreeHTML.py fits_I*.root fits_M*.root fits_T*.root -i muonFits -o fits_muon_trigger/
#python DevTools/TagAndProbe/python/dumpTagProbeTreeHTML.py fits_E*.root -i electronFits -o fits_electron_trigger/

#for path in IsoMu24 IsoMu24ORIsoTkMu24 IsoMu24ORIsoTkMu24ORMu50 IsoTkMu24 Mu8Leg Mu8LegDZ Mu8ORTkMu8Leg Mu8ORTkMu8LegDZ Mu17Leg Mu50 TkMu8Leg TkMu8LegDZ; do
#    python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py muon fits_muon_trigger "$path" "$path" --trig
#done
#
#for path in Ele12Leg Ele12LegDZ Ele23Leg Ele27Eta2p1WPTight Ele27WPTight; do
#    python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py electron fits_electron_trigger "$path" "$path" --trig
#done

for path in Pho18Leg Pho30Leg Pho60Leg Pho175; do
    python DevTools/TagAndProbe/python/plotTagProbeEfficiencies.py photon fits_photon_trigger "$path" "$path" --trig
done

