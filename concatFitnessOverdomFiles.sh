#!/bin/bash
# this script concatenates output produced by runFitnessOverdomScript.sh
# which itself runs c++ program runFitnessOverdomOptGtypeByPopt for specific Popt and Ntf values.
# in these runs, the output files are named "fitnessOverdomSummaryTable_b6_Ntf100_[Popt].txt",
# [Popt] is the popt value for the specific run (the names of the output files are created by the C++ code
# based on the input Ntf value)

# this will hold the concatenated files
fileForOutput="fitnessOverdomSummaryTable_b6_o0.05_Ntf100.txt"
echo -n "" > $fileForOutput

find . -name "fitnessOverdomSummaryTable_b6_Ntf100_*" -print0 | xargs -0 cat > $fileForOutput

# run a separate script to delete the concatenated files, after confirming that the concatenation worked
#
