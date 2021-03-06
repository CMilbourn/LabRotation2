#!/bin/bash

#Runs Full Analysis for one subject
#
#Skip:0, On:120, off:180
#NOTE that all the 'exits' have been commented out for now
 
#Print out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0 <sourcedata folder> <SubjectNumber> <MRIParameter> e.g [default] (default,paramA,paramB) [options] Output: .Feat Full analysis folder. e.g. ./one_subj_stats.sh /Users/colette/sourcedata/ sub-02 default 060"; 
 #exit 1 ;
fi

#INPUTS 
src=$1 #Input 1 is the source data - the folder with the datafiles in the format Users/name/sourcedata
subj=$2 #Input 2 is the subject ID in the format sub-01
MRIParam=$3 #Input 3 is the MRI parameters e.g. default
SkipNo=$4 #Input 4 is a specific SKIPNUMBER
#######################################################
#HOUSEKEEPING
#######################################################
#Print the Subject Number
echo "Subject Number:" ${subj}

#Print the MRI Parameter
#echo "MRI Parameter is:" ${MRIParam}

#Print Skip Number
echo "Skip Number:" ${SkipNo}

srcin=${src}/derivatives/${subj}/ #source in is the sourcedata/subj-0x/sub-0x_FullDefault.feat
srcout=${src}/derivatives/StatsOutput/${subj} #source out is the sourcedata/Derivatives/StatsOutput/subj-0x

#Check if directories exist
if ! [ -d "${src}/derivatives/StatsOutput" ]; then
mkdir ${src}/derivatives/StatsOutput #if the directory does not exist, then make it
fi
 
if ! [ -d "${srcout}" ]; then #check if the output directory exists, if not make it
mkdir ${srcout}
fi

#if there is an existing stats folder - rename it statsOLD
StatsFolder=${srcin}/${subj}_FullDefault.feat/stats
#if [ -d "${StatsFolder}" ]; then
#echo "Stats Folder Exists, re-naming to statsOLD"
#cp -r ${StatsFolder} ${srcin}/${subj}_FullDefault.feat/statsOLD #if the stats directory exists, then copy it and rename it statsOLD
#fi

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
echo "Subject Number:" ${subj}

#Print the MRI Parameter
#echo "MRI Parameter is:" ${MRIParam}

#Print Skip Number
echo "Skip Number:" ${SkipNo}
#INSERT CODE TO MODIFY TEMPLATE
#Copy the template .fsf file and put one called designDefault in the source out folder
cp template_one_subj_stats2.fsf ${srcout}/designStats.fsf
#Input feat directory
inputfeatdir=${srcin}/${subj}_FullDefault.feat
sed -i -e "s|INPUTFEATDIR|${inputfeatdir}|" ${srcout}/designStats.fsf #in the new designDefault.fsf file, find and replace 'INPUTDATA' with the inputdata file above

#Change Skip Number
skipno=${SkipNo}
sed -i -e "s|SKIPNUMBER|${skipno}|" ${srcout}/designStats.fsf #in the designDefault.fsf file change the Skip Number


#Run full FEAT analysis following design.Default.fsf
echo "Running Feat Stats ..."
feat ${srcout}/designStats.fsf 

#Check the stats folder exists
#outputdirectory

#Copy the stats folder & rename it
echo "Creating Stats Folder ..."
if [ -d "${srcout}/stats_${subj}_${SkipNo}" ]; then
	echo "Deleting existing stats folder"
	rm -r ${srcout}/stats_${subj}_${SkipNo}
fi
cp -r ${StatsFolder} ${srcout}/stats_${subj}_${SkipNo}

#re-name the files in that folder to end in _${subj}_${SkipNo}
cd ${srcout}/stats_${subj}_${SkipNo}
pwd 
ls
chmod u+w stats_${subj}_${SkipNo}
#chmod a+x ./*
#if ! [ -f "cope2_${subj}_${SkipNo}.nii.gz" ]; then
for file in *.nii.gz; do mv "${file}" "${file/.nii.gz/_${subj}_${SkipNo}.nii.gz}"; done
#fi

chmod u+w *.nii.gz
#create the CVRmap
echo "Creating CVRmap..."
echo "fslmaths cope2_${subj}_${SkipNo} -div cope6_${subj}_${SkipNo} CVRmap_${subj}_${SkipNo}"
fslmaths cope2_${subj}_${SkipNo} -div cope6_${subj}_${SkipNo} CVRmap_${subj}_${SkipNo}

#check if CVRmap_sub-0x exists, make one if not
if ! [ -d "${srcout}/CVRmaps_${subj}" ]; then #check if the output directory exists, if not make it
mkdir ${srcout}/CVRmaps_${subj}
fi
chmod u+w *.nii.gz
#Move copy of CVRMap_sub-0x_SkipNo to CVRmaps_sub-0x
echo "Moving CVRmap to subfolder..."
cp ${srcout}/stats_${subj}_${SkipNo}/CVRmap_${subj}_${SkipNo}.nii.gz ${srcout}/CVRmaps_${subj}

#Check if Zstats folder exists, make if not
if ! [ -d "${srcout}/Zstat2_${subj}" ]; then #check if the output directory exists, if not make it
mkdir ${srcout}/Zstat2_${subj}
fi
chmod u+w *.nii.gz
#Move copy of zstat to Zstat2_sub-0x
echo "Moving copy of zstat2 to subfolder..."
cp ${srcout}/stats_${subj}_${SkipNo}/zstat2_${subj}_${SkipNo}.nii.gz ${srcout}/Zstat2_${subj}

#sed -i -e "s|OUTPUTDIR|${srcout}|" ${srcout}/designStats.fsf #in the designDefault.fsf file find and replace 'OUTPUTDIR' with the outputdata file above

