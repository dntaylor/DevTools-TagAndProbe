#!/bin/bash

# Submit efficiency ntuple jobs on the mu+jet skim
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
fi

#farmoutAnalysisJobs $1-DATA \
#  --infer-cmssw-path \
#  --input-dbs-path=/SingleMuon/Run2016B-PromptReco-v2/MINIAOD \
#  --assume-input-files-exist \
#  ./muonTagAndProbeTreeFarmout_cfg.py  \
#  isMC=0 \
#  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'

farmoutAnalysisJobs $1-MC \
  --infer-cmssw-path \
  --input-dbs-path=/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM \
  --assume-input-files-exist \
  ./muonTagAndProbeTreeFarmout_cfg.py  \
  isMC=1  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'

