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
embedding_input="/star/data105/embedding/Run14_AuAu200_Pythia6_PicoDst/pt50_-1_0/out/st_physics_adc_15097061_raw_0000000.picoDst.root"

# Real data test files
real_inputs=(
  "/star/data40/reco/AuAu_200_production_low_2014/ReversedFullField/P18ih/2014/095/15095020/st_physics_15095020_raw_4500043.picoDst.root"
  "/star/data40/reco/AuAu_200_production_low_2014/ReversedFullField/P18ih/2014/166/15166044/st_physics_15166044_raw_3000019.picoDst.root"
  "/star/data46/reco/AuAu_200_production_low_2014/ReversedFullField/P18ih/2014/108/15108023/st_physics_15108023_raw_2000053.picoDst.root"
  "/star/data46/reco/AuAu_200_production_low_2014/ReversedFullField/P18ih/2014/154/15154021/st_physics_15154021_raw_5500048.picoDst.root"
)

if [[ "$mode" == "embedding" || -z "$mode" ]]; then
  echo "Running test on embedding..."
  root4star -l -b -q "StRoot/macros/runPicoHFJetMaker.C(\"$embedding_input\", \"test_embedding\", 0, \"picoDst\", true)"
fi

if [[ "$mode" == "data" || -z "$mode" ]]; then
  echo "Running tests on data..."
  for input_file in "${real_inputs[@]}"; do
    tag=$(basename "$input_file" .picoDst.root)
  root4star -l -b -q "StRoot/macros/runPicoHFJetMaker.C(\"$input_file\", \"test_data_${tag}\", 0, \"picoDst\", false)"
  done
fi
