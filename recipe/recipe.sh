#!/usr/bin/env bash

# CMSSW packages
pushd $CMSSW_BASE/src
git cms-merge-topic -u matteosan1:egm_tnp_80X
popd
