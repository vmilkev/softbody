#!/usr/bin/bash
echo CLEAN 'SAMPLE' FOLDER
rm -R sample
echo sample is removed
echo CLEAN 'SIMULATION' FOLDER
rm -R simulation
echo simulation is removed
rm sample_runtime.log
rm simulation_runtime.log
