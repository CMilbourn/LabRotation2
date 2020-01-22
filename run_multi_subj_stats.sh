#!/bin/bash
#OLD VERSIPON Nicholas Blockley 2019
# Run one_subj.sh to run full FEAT analysis 

codedir=`dirname $0` #assume that other scripts are in the same directory

#Print out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0 <sourcedata folder> <subjectnumber> e.g. ./run_one_subj_stats_test.sh /Users/colette/sourcedata sub-04"

fi

#src=$1 #Input 1 is the directory with the scripts in e.g.
# e.g. srcin=/Users/colette/sourcedata/ #path to sourcedata
src=$1 
echo "Sourcedata Folder:" ${1}
subj=$2 #Input 2 is the subject ID in the format sub-01
echo "Subject selected:" ${2}
echo "Number of arguments passed:" $# 
#MRIparam=$3
#cd ${src}


#echo $PWD
	
for n in $(seq -w 0 4 120);
do
	echo $n
	${codedir}/one_subj_stats.sh ${src} ${subj} default ${n}	
	#${codedir}/one_subj_stats.sh ${src} ${subj} ${MRIparam} ${n}
done
 
srcout=${src}/derivatives/StatsOutput/${subj}
cd ${srcout}/Zstat2_${subj}
echo "Changing to folder:" ${srcout}/Zstat2_${subj}
chmod u+w ./*
#echo "current directory is:"${PWD};

echo "Merging zstat2 Files..."
if ! [ -d "${srcout}/Zstat2_${subj}" ]; then #check if the output directory exists, if not make it
echo "File" ${srcout}/Zstat2_${subj} "NOT FOUND!"
fi
fslmerge -t zstat2_concat_${subj} $zstat2_${subj}*    

 


