#source me
rm -rf fits/*
python fitter.py muon -mc TnPTree_mc_muon.root -data TnPTree_data_muon.root  -mcLO TnPTree_mcLO_muon.root
mkdir -p fits/badFits
mv fits/badFit_* fits/badFits
python dumpTagProbeTreeHTML.py --data fits.root -i muonFits -o fits/
python plot.py HTTID "HTT ID"
python plot.py HTTISO "HTT Iso"
cp -r fits/ ~/www/ztt/2016/MuTau/
