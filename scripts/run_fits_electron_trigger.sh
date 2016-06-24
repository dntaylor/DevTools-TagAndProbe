# single electron
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=electron inputFileName=tagAndProbe/electron/data.root idName=passingHLTEle25Eta2p1
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=electron inputFileName=tagAndProbe/electron/data.root idName=passingHLTEle25Tight
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=electron inputFileName=tagAndProbe/electron/data.root idName=passingHLTEle35
# double electron
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=electron inputFileName=tagAndProbe/electron/data.root idName=passingHLTEle23Ele12Leg1
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=electron inputFileName=tagAndProbe/electron/data.root idName=passingHLTEle23Ele12Leg2
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=electron inputFileName=tagAndProbe/electron/data.root conditions=passingHLTEle23Ele12Leg2 idName=passingHLTEle23Ele12DZ
# electron electron
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=electron inputFileName=tagAndProbe/electron/data.root idName=passingHLTMu23Ele12ELeg
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=electron inputFileName=tagAndProbe/electron/data.root idName=passingHLTMu23Ele8ELeg
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=electron inputFileName=tagAndProbe/electron/data.root idName=passingHLTMu8Ele17ELeg
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=electron inputFileName=tagAndProbe/electron/data.root idName=passingHLTMu8Ele23ELeg
# tau electron
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=electron inputFileName=tagAndProbe/electron/data.root idName=passingHLTEle24Tau20LegSingleL1
