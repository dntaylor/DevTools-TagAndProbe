#!/usr/bin/env bash

# CMSSW packages
pushd $CMSSW_BASE/src

# TODO: update to 92X
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/ElectronTagAndProbe
# currently brings in a lot of junk
#git cms-merge-topic fcouderc:tnp_egm_80X_Moriond17_v1.0
# rebased on 8026p1
#git cms-merge-topic dntaylor:tnp_moriond_8026p1

popd
