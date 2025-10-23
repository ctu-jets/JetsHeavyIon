#!/bin/bash

setup 64b
cons

run_mode="$1"  # "embedding", "data", or empty for both

embedding_lists=(
    "filelists/embedding/pt3_5.list"
    "filelists/embedding/pt5_7.list"
    "filelists/embedding/pt7_9.list"
    "filelists/embedding/pt9_11.list"
    "filelists/embedding/pt11_15.list"
    "filelists/embedding/pt15_20.list"
    "filelists/embedding/pt20_25.list"
    "filelists/embedding/pt25_30.list"
    "filelists/embedding/pt30_40.list"
    "filelists/embedding/pt40_50.list"
    "filelists/embedding/pt50_-1.list"
)

real_data_list="filelists/pico_low_14.list"

if [[ "$run_mode" == "embedding" || -z "$run_mode" ]]; then
  for list_file in "${embedding_lists[@]}"; do
    ./submit/submit.sh "$list_file" embedding
    if [[ $? -ne 0 ]]; then
        echo "Error processing $list_file"
        exit 1
    fi
  done
fi

if [[ "$run_mode" == "data" || -z "$run_mode" ]]; then
  ./submit/submit.sh "$real_data_list" data
  if [[ $? -ne 0 ]]; then
      echo "Error processing $real_data_list"
      exit 1
  fi
fi

