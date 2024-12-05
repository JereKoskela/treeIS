#!/bin/bash

date "+%H:%M:%S"
../huw infinitesites-55.dat 5.0 100000 100000 1 0.0 > ../Results/cost-55-huw.txt
date "+%H:%M:%S"
../huw infinitesites-550.dat 5.0 100000 100000 1 0.0 > ../Results/cost-550-huw.txt
date "+%H:%M:%S"
../huw infinitesites-5500.dat 5.0 100000 100000 1 0.0 > ../Results/cost-5500-huw.txt
date "+%H:%M:%S"

../sd infinitesites-55.dat 5.0 100000 100000 1 0.0 > ../Results/cost-55-sd.txt
date "+%H:%M:%S"
../sd infinitesites-550.dat 5.0 100000 100000 1 0.0 > ../Results/cost-550-sd.txt
date "+%H:%M:%S"
../sd infinitesites-5500.dat 5.0 100000 100000 1 0.0 > ../Results/cost-5500-sd.txt
date "+%H:%M:%S"

../gt infinitesites-55.dat 5.0 100000 100000 1 0.0 > ../Results/cost-55-gt.txt
date "+%H:%M:%S"
../gt infinitesites-550.dat 5.0 100000 100000 1 0.0 > ../Results/cost-550-gt.txt
date "+%H:%M:%S"
../gt infinitesites-5500.dat 5.0 100000 100000 1 0.0 > ../Results/cost-5500-gt.txt
date "+%H:%M:%S"


