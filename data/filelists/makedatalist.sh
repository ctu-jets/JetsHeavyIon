#!/bin/bash

rm -f test_data.list
get_file_list.pl -keys path,filename -cond trgsetupname=AuAu_200_production_2014,library=SL20d,production=P18ih,filetype=daq_reco_picoDst,storage!=HPSS -limit 0 -delim "/" | \
sed 's|^/home/starlib|root://xrdstar.rcf.bnl.gov:1095//home/starlib|g' > pico_low_14.list

