#!/usr/bin/env bash

# CMSSW packages
pushd $CMSSW_BASE/src
# currently brings in a lot of junk
#git cms-merge-topic fcouderc:tnp_egm_80X_Moriond17_v1.0
git cms-merge-topic dntaylor:tnp_moriond_8026p1
popd
