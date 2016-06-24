# single muon
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingIsoMu20
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingIsoTkMu20
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingIsoMu20ORIsoTkMu20
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingIsoMu22
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingIsoTkMu22
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingIsoMu22ORIsoTkMu22
# double muon
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingMu17
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingMu8
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingTkMu8
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingMu8ORTkMu8
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root conditions=passingMu8 idName=passingMu8DZ
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root conditions=passingTkMu8 idName=passingTkMu8DZ
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root conditions=passingMu8ORTkMu8 idName=passingMu8ORTkMu8DZ
# muon eg
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingMu17Ele12MLeg
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingMu23Ele12MLeg
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingMu23Ele8MLeg
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingMu8Ele17MLeg
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingMu8Ele23MLeg
# muon tau
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingMu17Tau20MLegSingleL1
cmsRun DevTools/TagAndProbe/python/oldfitter.py object=muon inputFileName=tagAndProbe/muon/data.root idName=passingMu19Tau20MLegSingleL1
