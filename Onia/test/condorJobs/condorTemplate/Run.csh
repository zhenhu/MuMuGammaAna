#!/bin/csh
source /cvmfs/cms.cern.ch/cmsset_default.csh
cd CMSSWDIR
cmsenv
cd ${_CONDOR_SCRATCH_DIR}
root -b -q -l FILENAME+
