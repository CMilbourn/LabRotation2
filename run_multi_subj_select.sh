#!/bin/bash
# Nicholas Blockley 2019
# Run one_subj.sh to run full FEAT analysis 

codedir=`dirname $0` #assume that other scripts are in the same directory

#Print out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0  <sourcedata folder> <fromsubj> <tosubj> e.g. ./run_one_subj.sh /Users/colette/sourcedata/ 2 3"
fi

#src=$1 #input 1 is the directory with the scripts in e.g.
#/Users/colette/sourcedata/code/StatsSkip/OneSubjectFullAnalysis
# e.g. srcin=/Users/colette/sourcedata/ #path to sourcedata
src=$1
fromsubj=$2
tosubj=$3
#cd ${src}

#echo $PWD

#e.g. ./one_subj.sh /Users/colette/sourcedata/ sub-01 default
for n in {${fromsubj}..${tosubj}}
do

	subj_id=sub-`printf "%02d\n" $n`
	echo $subj_id
	#echo
	${codedir}/one_subj.sh ${src} ${subj_id} default
	${codedir}/run_one_subj_stats.sh ${src} ${subj_id}
#chmod u+w ./*
#chmod a+x ./*
done