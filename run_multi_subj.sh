#!/bin/bash
# Run multi_subj.sh to run full FEAT analysis 

codedir=`dirname $0` #assume that other scripts are in the same directory

#Print out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0  <sourcedata folder> <MRIParam> e.g. ./run_multi_subj.sh /Users/colette/sourcedata/ paramA"
fi

#src=$1 #input 1 is the directory with the scripts in e.g.
#/Users/colette/sourcedata/code/StatsSkip/OneSubjectFullAnalysis
# e.g. srcin=/Users/colette/sourcedata/ #path to sourcedata
src=$1 
MRIParam=$2
#cd ${src}

#echo $PWD

#e.g. ./one_subj.sh /Users/colette/sourcedata/ sub-01 default
for n in {1..6}
do

	subj_id=sub-`printf "%02d\n" $n`
	echo "Subject ID:"$subj_id
	
	echo "Running multi_subj.sh for: ${src} ${subj_id} ${MRIParam}"
	${codedir}/multi_subj.sh ${src} ${subj_id} ${MRIParam}
	echo "Running run_multi_subj_stats.sh for: ${src} ${subj_id} ${MRIParam}"
	${codedir}/run_multi_subj_stats.sh ${src} ${subj_id} ${MRIParam}
	
#chmod u+w ./*
#chmod a+x ./*
done