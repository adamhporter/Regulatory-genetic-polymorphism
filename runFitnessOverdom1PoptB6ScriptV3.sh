#!/bin/bash
# This file is called runFitnessOverdom1PoptB6ScriptV3.sh
# It re-runs some Popt by TF combinations that timed out in runFitnessOverdom1PoptB6ScriptV2.sh

bitstringLen=6
Ntfsat=100
omega=0.05

Popt=105                 #arbitrary integer value from 0 to PoptTotalSteps; the program will divide this by PoptTotalSteps to get the actual Popt
PoptTotalSteps=1000      #integer value, usually 1000 since that’s what’s run in runFitnessOverdomB6Script.sh
lowRange=0      #these are all integers and ignored by the program; they're only here for continuity with the other batch scripts
stepSize=5
PoptValsPerRun=1




# this generates the file fitnessOverdomSummaryTable_b6_Ntf100_Popt0.4_tf6.txt

let Popt=400
bsub -q long -n 16 -W 720:00 -R span[hosts=1] -R rusage[mem=16]  /home/ap73a/simulations/fitnessOverdomOptGtypeByPopt/fitnessOverdomOptGtypeByPopt "$bitstringLen" "$Ntfsat" "$omega" "$Popt" "$Popt" "$PoptTotalSteps" "$stepSize" 6 6

# this generates the file fitnessOverdomSummaryTable_b6_Ntf100_Popt0.4_tf7.txt
bsub -q long -n 16 -W 720:00 -R span[hosts=1] -R rusage[mem=16]  /home/ap73a/simulations/fitnessOverdomOptGtypeByPopt/fitnessOverdomOptGtypeByPopt "$bitstringLen" "$Ntfsat" "$omega" "$Popt" "$Popt" "$PoptTotalSteps" "$stepSize" 7 7

