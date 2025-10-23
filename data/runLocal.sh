#!/bin/bash

# Setup (optional, depending on environment)
# export FASTJET='/gpfs/mnt/gpfs01/star/pwg/svomich/install/x64/fastjet-install'
# export FASTJET_CONTRIB='/gpfs/mnt/gpfs01/star/pwg/svomich/install/x64/fjcontrib-1.046'
# export PATH=$FASTJET/bin:$PATH
# export LD_LIBRARY_PATH=$FASTJET/lib:$LD_LIBRARY_PATH
# export CPLUS_INCLUDE_PATH=$FASTJET/include:$CPLUS_INCLUDE_PATH

setup 64b
cons

mode="$1"  # "embedding" or "data" or empty for both

# Embedding test file
# embedding_input="testPythia6picoDsts_pt50_-1.list"
 embedding_input="/star/data105/embedding/Run14_AuAu200_Pythia6_PicoDst/pt50_-1_0/out/st_physics_adc_15097061_raw_0000000.picoDst.root"

# Real data test files
real_input="filelists/test_data.list"
# real_input="root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/AuAu_200_production_2014/ReversedFullField/P18ih.SL20d/2014/080/15080057/st_physics_15080057_raw_3000008.picoDst.root"


if [[ "$mode" == "embedding" || -z "$mode" ]]; then
  echo "Running test on embedding..."
  root4star -l -b -q "StRoot/macros/runPicoHFJetMaker.C(\"$embedding_input\", \"test_embedding\", 0, \"picoDst\", true)"
fi

if [[ "$mode" == "data" || -z "$mode" ]]; then
  echo "Running tests on data..."
  root4star -l -b -q "StRoot/macros/runPicoHFJetMaker.C(\"$real_input\", \"test_data\", 0, \"picoDst\", false)"
fi
