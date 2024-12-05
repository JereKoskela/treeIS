#!/bin/bash

date "+%H:%M:%S"

../sd 50 20 0.1 10000 10000 0 0.1 > ../Results/likelihood-10k-results-50-re.txt
../sd 50 20 0.2 10000 10000 0 0.1 >> ../Results/likelihood-10k-results-50-re.txt
../sd 50 20 0.3 10000 10000 0 0.1 >> ../Results/likelihood-10k-results-50-re.txt
../sd 50 20 0.4 10000 10000 0 0.1 >> ../Results/likelihood-10k-results-50-re.txt
../sd 50 20 0.5 10000 10000 0 0.1 >> ../Results/likelihood-10k-results-50-re.txt
../sd 50 20 0.6 10000 10000 0 0.1 >> ../Results/likelihood-10k-results-50-re.txt
../sd 50 20 0.7 10000 10000 0 0.1 >> ../Results/likelihood-10k-results-50-re.txt
../sd 50 20 0.8 10000 10000 0 0.1 >> ../Results/likelihood-10k-results-50-re.txt
../sd 50 20 0.9 10000 10000 0 0.1 >> ../Results/likelihood-10k-results-50-re.txt

date "+%H:%M:%S"

../sd 50 20 0.1 300 300 0 0.1 > ../Results/likelihood-mean-results-50-re.txt
../sd 50 20 0.2 343 343 0 0.1 >> ../Results/likelihood-mean-results-50-re.txt
../sd 50 20 0.3 443 443 0 0.1 >> ../Results/likelihood-mean-results-50-re.txt
../sd 50 20 0.4 586 586 0 0.1 >> ../Results/likelihood-mean-results-50-re.txt
../sd 50 20 0.5 760 760 0 0.1 >> ../Results/likelihood-mean-results-50-re.txt
../sd 50 20 0.6 958 958 0 0.1 >> ../Results/likelihood-mean-results-50-re.txt
../sd 50 20 0.7 1170 1170 0 0.1 >> ../Results/likelihood-mean-results-50-re.txt
../sd 50 20 0.8 1390 1390 0 0.1 >> ../Results/likelihood-mean-results-50-re.txt
../sd 50 20 0.9 1613 1613 0 0.1 >> ../Results/likelihood-mean-results-50-re.txt

date "+%H:%M:%S"

../sd 50 20 0.1 100 10000 0 0.1 > ../Results/likelihood-variable-results-50-re.txt
../sd 50 20 0.2 100 10000 0 0.1 >> ../Results/likelihood-variable-results-50-re.txt
../sd 50 20 0.3 100 10000 0 0.1 >> ../Results/likelihood-variable-results-50-re.txt
../sd 50 20 0.4 100 10000 0 0.1 >> ../Results/likelihood-variable-results-50-re.txt
../sd 50 20 0.5 100 10000 0 0.1 >> ../Results/likelihood-variable-results-50-re.txt
../sd 50 20 0.6 100 10000 0 0.1 >> ../Results/likelihood-variable-results-50-re.txt
../sd 50 20 0.7 100 10000 0 0.1 >> ../Results/likelihood-variable-results-50-re.txt
../sd 50 20 0.8 100 10000 0 0.1 >> ../Results/likelihood-variable-results-50-re.txt
../sd 50 20 0.9 100 10000 0 0.1 >> ../Results/likelihood-variable-results-50-re.txt

date "+%H:%M:%S"

../sd 50 20 0.1 100 100 0 0.1 > ../Results/likelihood-100-results-50-re.txt
../sd 50 20 0.2 100 100 0 0.1 >> ../Results/likelihood-100-results-50-re.txt
../sd 50 20 0.3 100 100 0 0.1 >> ../Results/likelihood-100-results-50-re.txt
../sd 50 20 0.4 100 100 0 0.1 >> ../Results/likelihood-100-results-50-re.txt
../sd 50 20 0.5 100 100 0 0.1 >> ../Results/likelihood-100-results-50-re.txt
../sd 50 20 0.6 100 100 0 0.1 >> ../Results/likelihood-100-results-50-re.txt
../sd 50 20 0.7 100 100 0 0.1 >> ../Results/likelihood-100-results-50-re.txt
../sd 50 20 0.8 100 100 0 0.1 >> ../Results/likelihood-100-results-50-re.txt
../sd 50 20 0.9 100 100 0 0.1 >> ../Results/likelihood-100-results-50-re.txt

date "+%H:%M:%S"
