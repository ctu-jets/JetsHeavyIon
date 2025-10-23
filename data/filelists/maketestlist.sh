#!bin/bash

nfiles=${5:-10}

rm -f test_data.list
get_file_list.pl -keys path,filename -cond trgsetupname=AuAu_200_production_2014,library=SL20d,production=P18ih,filetype=daq_reco_picoDst,storage!=HPSS -limit 105 -delim "/" | sed 's|^/home/starlib|root://xrdstar.rcf.bnl.gov:1095//home/starlib|g' > tmp.list
tail -n $nfiles tmp.list > test_data.list
rm -f tmp.list 
