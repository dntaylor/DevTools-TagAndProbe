#!/bin/bash
DATE=`date +%Y-%m-%d`
if [ "$1" == "" ]; then
    NAME="v1"
else
    NAME="$1"
fi
submit_job.py crabSubmit --sampleList $CMSSW_BASE/src/DevTools/TagAndProbe/data/datasetList_TagAndProbe_MC.txt --filesPerJob 10 "$DATE"_DevTools_TagAndProbe_Photon_80X_"$NAME" DevTools/TagAndProbe/test/photonTagAndProbeTree_cfg.py isMC=1
submit_job.py crabSubmit --sampleList $CMSSW_BASE/src/DevTools/TagAndProbe/data/datasetList_TagAndProbe_Data_Electron.txt --lumisPerJob 600 --applyLumiMask "Collisions16" "$DATE"_DevTools_TagAndProbe_Photon_80X_"$NAME" DevTools/TagAndProbe/test/photonTagAndProbeTree_cfg.py isMC=0
