#!/bin/bash
# Nicholas Blockley 2019
# Run one_subj.sh to run full FEAT analysis 

codedir=`dirname $0` #assume that other scripts are in the same directory

#Print out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0 <sourcedata folder> <subjectnumber> e.g. ./run_one_subj_stats.sh Users/colette/sourcedata sub-01"
fi

#src=$1 #Input 1 is the directory with the scripts in e.g.
# e.g. srcin=/Users/colette/sourcedata/ #path to sourcedata
src=$1 
subj=$2 #Input 2 is the subject ID in the format sub-01
echo "Source Data folder is:" $1
echo "Subject Number is:" $2
#cd ${src}
#echo $PWD
#e.g. ./one_subj.sh /Users/colette/sourcedata/ sub-01 default
	
for n in $(seq -w 0 4 120);
do
	#((SkipNumber= n))
	echo $n
	${codedir}/one_subj_stats.sh ${src} ${subj} default ${n}
#done	
#chmod u+w ./*
#chmod a+x ./*
done

