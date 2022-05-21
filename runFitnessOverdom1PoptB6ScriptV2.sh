#!/bin/bash
# This file is called runFitnessOverdom1PoptB6Script.sh
# Some individual 6-bit Popt runs can take weeks, so this script take breaks it into pieces for a subset of values for Tf0’s dosage

#  To determine what Popt runs require this script, first use runFitnessOverdomB6Script.sh, with the time cap set to 72 hrs or so
#  Runs that time out will have file sizes of 0, and those will be the ones to re-run here

bitstringLen=6
Ntfsat=100
omega=0.05

Popt=405			       #integer value
PoptTotalSteps=1000      #integer value, usually 1000 since that’s what’s run in runFitnessOverdomB6Script.sh
lowRange=0      #these are all integers
stepSize=5
PoptValsPerRun=1

tf0valuesPerOutputFile=2		# = a number divisible by maxTFvalue, including 1
#maxTfValue=64              # = 2^bitstringLen

initialTf0value=0
finalTf0value=63                                          # up to maximum = 2^bitstringLen - 1
currentTfValue=0
let currentTFvalue=initialTf0value
nextTFvalue=0
let nextTFvalue=$initialTf0value+$tf0valuesPerOutputFile-1


for((tfVal=0;tfVal<$finalTf0value;tfVal+=$tf0valuesPerOutputFile));do
        let nextTFvalue=$tfVal+1
        echo "tfVal="$tfVal   " nextTFvalue="$nextTFvalue
	bsub -q long -n 16 -W 300:00 -R span[hosts=1] -R rusage[mem=16]  /home/ap73a/simulations/fitnessOverdomOptGtypeByPopt/fitnessOverdomOptGtypeByPopt "$bitstringLen" "$Ntfsat" "$omega" "$Popt" "$Popt" "$PoptTotalSteps" "$stepSize" "$tfVal" "$nextTFvalue"
   
	done
 
