#!/bin/bash
# Run multi_subj.sh to run full FEAT analysis 

codedir=`dirname $0` #assume that other scripts are in the same directory

#Print out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0  <sourcedata folder> e.g. time ./run_pipeline.sh /Users/colette/sourcedata/ "
fi

src=$1
#Check which subjects to run through
#while true;
#do
#    read -r -p "Have you updated the subject numbers to run through in multi_subj.sh Yes or no? " response   
#    if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
#    then
#        echo "You chose yes, continuing to next question"
 #   else
    	
#        exit 0
#    fi
#done

#while true;
#do
#    read -r -p "Have you updated the subject numbers to run in this script Yes or no? " response   
#    if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
#    then
#        echo "You chose yes, continuing to now run the scripts"
#    else
#        exit 0
#    fi
#done
echo "testing script now"

echo "~ Running run_multi_subj.sh"
#echo "Running run_multi_subj.sh for ... ~~~default~~~"
${codedir}/run_multi_subj.sh ${src} default
echo "Running run_multi_subj.sh for ... ~~~paramA~~~"
${codedir}/run_multi_subj.sh ${src} paramA
#/Users/colette/sourcedata/code/StatsSkip/Scripts/multi_subj/run_multi_subj.sh ${src} paramA
echo "Running run_multi_subj.sh for ... ~~~paramB~~~"
${codedir}/run_multi_subj.sh ${src} paramB
#/Users/colette/sourcedata/code/StatsSkip/Scripts/multi_subj /run_multi_subj.sh ${src} paramB

#echo "~~~ End of script run_all_multi_subj.sh ~~~"

#cd ${src}

# Run Do_analysis for a range of subjects below

#Check which subjects to run through
echo "Running do_analysis ..."
for n in {10..14}
do

	subj_id=sub-`printf "%02d\n" $n`
	echo "Subject ID:"$subj_id
	
	echo "Running do_analysis_edit_3.sh for: ${src} ${subj_id} "
	${codedir}/do_analysis_edit_3.sh ${src} ${subj_id}

done

echo "~~~ End of run_pipeline.sh ~~~"
echo "Now run cvr_processing.m and zstat2_processing.m"