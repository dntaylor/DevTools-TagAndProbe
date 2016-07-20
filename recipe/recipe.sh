#!/usr/bin/env bash

# CMSSW packages
pushd $CMSSW_BASE/src
git cms-merge-topic fcouderc:tnp_egm_80X
popd
