#!/bin/bash
#OLD VERSIPON Nicholas Blockley 2019
# Run one_subj.sh to run full FEAT analysis 

codedir=`dirname $0` #assume that other scripts are in the same directory

#Print out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0 <sourcedata folder> <subjectnumber> <MRIParam> Output: Feat Stats folders for one subject, skip numbers 000 to 120, 4sec intervals e.g. ./run_multi_subj_stats.sh /Users/colette/sourcedata sub-02 paramA"

fi
echo "Running run_multi_subj_stats.sh for..."
#src=$1 #Input 1 is the directory with the scripts in e.g.
# e.g. srcin=/Users/colette/sourcedata/ #path to sourcedata
src=$1 
echo "Sourcedata Folder:" ${1}
subj=$2 #Input 2 is the subject ID in the format sub-01
MRIParam=$3 #Input 3 is the MRI parameter default, paramA or paramB
echo "Subject Number:" ${2}
#echo "Number of arguments passed:" $# 

#echo $PWD
#Loop through skip number 000 in steps of 4 seconds through to 120, n=Skip Number	
for n in $(seq -w 0 4 120);
do
	echo $n
	${codedir}/multi_subj_stats.sh ${src} ${subj} ${MRIParam} ${n}
done
 
srcout=${src}/derivatives/StatsOutput/${subj}/${subj}_${MRIParam}
srcout=${src}/derivatives/StatsOutput/${subj}/${subj}_${MRIParam}
cd ${srcout}/Zstat2_${subj}_${MRIParam}
echo "Changing to folder:" ${srcout}/Zstat2_${subj}_${MRIParam}
chmod u+w ./*
#echo "current directory is:"${PWD};

echo "Merging zstat2 Files..."
if ! [ -d "${srcout}/Zstat2_${subj}_${MRIParam}" ]; then #check if the output directory exists, if not make it
echo "File" ${srcout}/Zstat2_${subj}_${MRIParam} "NOT FOUND!"
fi
fslmerge -t ${srcout}/Zstat2_${subj}_${MRIParam}/zstat2_concat_${subj}_${MRIParam} ${srcout}/Zstat2_${subj}_${MRIParam}/zstat2_${subj}*    

echo "~ ~ ~ End of run_multi_subj_stats.sh for... Subject:" $2 
