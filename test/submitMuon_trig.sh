#!/bin/bash

# Submit efficiency ntuple jobs on the mu+jet skim
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
fi

farmoutAnalysisJobs $1-triggerEffi \
  --infer-cmssw-path \
  --input-dbs-path=/SingleMuon/Run2016B-PromptReco-v2/MINIAOD \
  --assume-input-files-exist \
  ./muonTriggerTreeFarmout_cfg.py  \
  isMC=0 \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'
