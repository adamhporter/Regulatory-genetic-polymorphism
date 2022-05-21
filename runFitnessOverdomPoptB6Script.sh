#!/bin/bash
# file is called runFitnessOverdomB6Script.sh
# runs in the loop can have very diferent runtimes, up to 20 days for the omega=0.2 case.
#  So, set this to 72 hrs to get the ones that run at fast to moderate rates, and then run the others individually by splitting
#   the maximization into concurrently running pieces 
lowRange=0      #these are all integers
highRange=1000
stepSize=5
PoptValsPerRun=1
bitstringLen=6
Ntfsat=100
omega=0.05
totalSteps=1000
highPopt=0
currentPoptValue=$lowRange

while [ $currentPoptValue -le $highRange ]; do  # -le is <=
        bsub -q long -n 16 -W 70:00 -R span[hosts=1] /home/ap73a/simulations/fitnessOverdomOptGtypeByPopt/fitnessOverdomOptGtypeByPopt "$bitstringLen" "$Ntfsat" "$omega" "$currentPoptValue" "$currentPoptValue" "$totalSteps" "$stepSize"
        echo "currentPoptVal="$currentPoptValue
        echo "highPopt="$highPopt
        let currentPoptValue=currentPoptValue+stepSize
    done

