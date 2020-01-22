#!/bin/bash

#Runs Full Analysis for one subject
#Inputs are: 1.Subject Number <01> 2. MRIParameter <default> (default,paramA,paramB)
#Skip:0, On:120, off:180
#NOTE that all the 'exits' have been commented out for now


#Print out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0 <sourcedata folder> <SubjectNumber> <MRIParameter> e.g [default] (default,paramA,paramB) [options] Output: .Feat Full analysis folder. e.g. ./one_subj.sh /Users/colette/sourcedata/ sub-01 default"; 
 #exit 1 ;
fi

#INPUTS 
src=$1 #Input 1 is the source data - the folder with the datafiles in the format Users/name/sourcedata
subj=$2 #Input 2 is the subject ID in the format sub-01
MRIParam=$3 #Input 3 is the MRI parameters e.g. default

#######################################################
#HOUSEKEEPING
#######################################################

srcin=${src}/${subj} #source in is the sourcedata/subj-0x
srcout=${src}/derivatives/${subj} #source out is the sourcedata/derivatives/subj-0x

#Check if directories exist
if ! [ -d "${src}/derivatives/" ]; then
mkdir ${src}/derivatives/ #if the directory does not exist, then make it
fi

if ! [ -d "${srcout}" ]; then #check if the output directory exists, if not make it
mkdir ${srcout}
fi


#File path to the directory with your data in it:
#PWD=/Users/colette/sourcedata/sub-${subj}/func
#echo $PWD

#Change directory to that PWD
#cd $PWD
#echo current folder to double check
#echo "File Path is:" `pwd`

#Print the Number of arguments
#echo "Number of arguments passed: ${#}"
 
#Print the Subject Number
#echo "Subject Number:" ${subj}

#Print the MRI Parameter
#echo "MRI Parameter is:" ${MRIParam}

#INSERT CODE TO MODIFY TEMPLATE
#Copy the template .fsf file and put one called designDefault in the source out folder
cp template_one_subj2.fsf ${srcout}/designDefault.fsf
inputdata=${srcin}/func/${subj}_task-hyper_acq-default_asl #input data is sourcedata/func/subj-0x_task-hyper_acq-default_asl
sed -i -e "s|INPUTDATA|${inputdata}|" ${srcout}/designDefault.fsf #in the new designDefault.fsf file, find and replace 'INPUTDATA' with the inputdata file above
outputdata=${srcout}/${subj}_FullDefault.feat #output data is named sourcedata/subj-0x_FUllDefault.feat
sed -i -e "s|OUTPUTDIR|${outputdata}|" ${srcout}/designDefault.fsf #in the designDefault.fsf file find and replace 'OUTPUTDIR' with the outputdata file above

#Run full FEAT analysis following design.Default.fsf
echo "Running FEAT..."
feat ${srcout}/designDefault.fsf 

