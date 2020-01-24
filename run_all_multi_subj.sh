#!/bin/bash
# Run multi_subj.sh to run full FEAT analysis 

codedir=`dirname $0` #assume that other scripts are in the same directory

#Print out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0  <sourcedata folder> e.g. ./run_all_multi_subj.sh /Users/colette/sourcedata/ "
fi

src=$1

echo "Running run_multi_subj.sh for ... ~~~default~~~"
${codedir}/run_multi_subj.sh ${src} default
echo "Running run_multi_subj.sh for ... ~~~paramA~~~"
${codedir}/run_multi_subj.sh ${src} paramA
echo "Running run_multi_subj.sh for ... ~~~paramB~~~"
${codedir}/run_multi_subj.sh ${src} paramB

echo "~~~ End of script run_all_multi_subj.sh ~~~"