#!/usr/bin/bash
echo RUNNING 'SIMULATION' PROJECT: soil compression simulation
../bin/clusterdem.exe param_data/param.dat sample/sample.smpl
#echo simulation study is completed

perl -e 'print "Checking the result: "'
perl check_tests.pl simulation_runtime.log "SAMPLE PROJECT COMPLETED"
perl -e 'print "Sample simulation file is created succesfully!\n" if -f "simulation/simulation.smpl"'
