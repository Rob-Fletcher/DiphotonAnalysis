#!/bin/bash

START=`pwd`

# Change to the compilation directory
cd packages/

# Setup ATLAS local ROOT base if necessary
[ "`compgen -a | grep localSetupROOT`x" == "x" ] \
   && echo "Going to set up ATLAS local ROOT base from cvmfs" \
   && export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase \
   && source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet 1>/dev/null 2>&1

# Setup FAX if not in PATH
[ "${ROOTSYS}x" == "x" ] \
   && echo "Going to set up FAX from cvmfs" \
   && export RUCIO_ACCOUNT="$USER" \
   && localSetupFAX --quiet 1>/dev/null 2>&1

# RootCore configs
export CERN_USER="$USER"
export ROOTCORE_NCPUS=4
export ROOTCORE_AUTHOR="HGAM_USER_TAG"
 
base="Base,HGAM_BASE_VERSION"
[ -d "RootCoreBin" ] \
  && base=""
 
rcSetup $base

# A few useful aliases after things are setup
alias p='localSetupPandaClient --noAthenaCheck'
alias v='voms-proxy-init -voms atlas'

cd $START
