#!/bin/bash
# this script concatenates output produced by runFitnessOverdom1PoptB6ScriptV2.sh under settings that specify TF dosage
# which itself runs instances of the c++ program runFitnessOverdomOptGtypeByPopt

# first, create a new file called fitnessOverdomSummaryTable_b[x]6_Ntf[y]_Popt[z].txt, where x=bitstring length,
# y = Ntf value and z = Popt value

# In this example x=6,y=100, z=0.570.
#
# this will hold the concatenated files
fileForOutput="fitnessOverdomSummaryTable_b6_Ntf100_Popt0.570.txt"
echo -n "" > $fileForOutput

find . -name "fitnessOverdomSummaryTable_b6_Ntf100_Popt0.57_tf*" -print0 | xargs -0 cat > $fileForOutput

# run a separate script to delete the concatenated files, after confirming that the concatenation worked
#
