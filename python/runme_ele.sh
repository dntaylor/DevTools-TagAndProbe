#source me
rm -rf fits/*
python fitter.py electron -mc TnPTree_mc_electron.root -data TnPTree_data_electron.root  -mcLO TnPTree_mc_electron.root
mkdir -p fits/badFits
mv fits/badFit_* fits/badFits
python dumpTagProbeTreeHTML.py --data fits.root -i electronFits -o fits/
#python plot.py MVAID80 "MVAID 80"
#python plot.py MVAID90 "MVAID 90"
python plot.py HTTID "HTT ID"
python plot.py HTTISOwrtID "HTT ISO"
cp -r fits/ ~/www/ztt/2016/ETau/
