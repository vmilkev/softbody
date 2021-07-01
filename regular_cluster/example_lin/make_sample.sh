#!/usr/bin/bash
echo RUNNING 'SAMPLE' PROJECT: soil sample fabrication
../bin/clusterdem.exe param_data/param.dat
#echo sample fabrication is completed

perl -e 'print "Checking the result: "'
perl check_tests.pl sample_runtime.log "SAMPLE PROJECT COMPLETED"
perl -e 'print "Sample file is created succesfully!\n" if -f "sample/sample.smpl"'
