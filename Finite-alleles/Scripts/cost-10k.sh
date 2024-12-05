#!/bin/bash

date "+%H:%M:%S"
../sd 50 20 0.5 10000 10000 1 0.0 > ../Results/cost-50-sd.txt
date "+%H:%M:%S"
../sd 500 20 0.5 10000 10000 1 0.0 > ../Results/cost-500-sd.txt
date "+%H:%M:%S"
../sd 5000 20 0.5 10000 10000 1 0.0 > ../Results/cost-5000-sd.txt
date "+%H:%M:%S"

../gt 50 10000 20 0.5 1000 > ../Results/cost-50-gt.txt
date "+%H:%M:%S"
../gt 500 10000 20 0.5 1000 > ../Results/cost-500-gt.txt
date "+%H:%M:%S"
../gt 5000 10000 20 0.5 1000 > ../Results/cost-5000-gt.txt
date "+%H:%M:%S"
