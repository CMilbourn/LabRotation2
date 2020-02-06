# !/bin/bash
# do_analysis_applyxfm4d.sh sourcedata subj_no
#Prints out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0 <sourcedata folder> <SubjectNumber> e.g. time ./do_analysis_zstat2.sh /Users/colette/sourcedata/ sub-01"; 
 #exit 1 ;
fi
# 
src=$1
echo "sourcefolder in is $1"
subj=$2 
echo "Subject Number is $2"
PVthresh=0.5
srcout=/Users/colette/sourcedata/derivatives/
#MRIParam='paramA' # 'default paramA paramB'

#${srcout}/${subj}/${subj}_${MRIParam}_analysis
#${srcout}/StatsOutput/Zstat2maps_all/Zstat2maps_concat_${subj}_${MRIParam}.nii.gz -mul 100 ${srcout}/StatsOutput/CVRmaps_percentage/CVRmaps_concat_${subj}_${MRIParam}_percentage.nii.gz

# Make CVR maps ##
## Make directory for Percentage CVRmaps ##
if ! [ -d "${srcout}/StatsOutput/Zstat2_percentage" ]; then #check if the output directory exists, if not make it
mkdir ${srcout}/StatsOutput/Zstat2_percentage
echo "Making Zstats_percentage folder:${srcout}/StatsOutput/Zstat2_percentage"
fi 

echo "Making Zstat2 percentages for $2..."

#${FSLDIR}/bin/fslmaths ${srcout}/StatsOutput/Zstat2_all/zstat2_concat_sub-01_default.nii.gz -mul 100 ${srcout}/StatsOutput/Zstat2_percentage/Zstat2_concat_sub-01_default_percentage.nii.gz


## Loop to make CVRmaps percentages
MRIParam='default paramA paramB'
for MRIParam in ${MRIParam}
do 
${FSLDIR}/bin/fslmaths ${srcout}/StatsOutput/Zstats2_all/zstat2_concat_${subj}_${MRIParam}.nii.gz -mul 100 ${srcout}/StatsOutput/Zstat2_percentage/Zstat2_concat_${subj}_${MRIParam}_percentage.nii.gz
#/Users/colette/sourcedata/derivatives/StatsOutput/Zstats2_all/zstat2_concat_sub-01_default.nii.gz
echo "... Zstat2_concat_${subj}_${MRIParam}_percentage.nii.gz"
done 

echo "Applying Transforms "

MRIParam='default paramA paramB'
for MRIParam in ${MRIParam}
do 
${FSLDIR}/bin/flirt -in ${srcout}/StatsOutput/Zstat2_percentage/Zstat2_concat_${subj}_${MRIParam}_percentage.nii.gz -ref ${srcout}/${subj}/${subj}_MZeroScan_brain.nii.gz -out ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_${MRIParam}_Zstat2_MZeroScan -init ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_${MRIParam}_2MZeroScan.mat -applyxfm -v
echo "Zstat2_concat_${subj}_${MRIParam}_percentage.nii.gz + ${subj}_${MRIParam}_2MZeroScan ==> ${subj}_${MRIParam}_Zstat2_MZeroScan..."
done 
 

#applyxfm4d Zstat2_concat_sub-01_default_percentage.nii.gz ../../sub-01/sub-01_MZeroScan_brain.nii.gz CVRmaps_concat_sub-01_default_MZeroScan.nii.gz ../../sub-01/sub-01_default_analysis/sub-01_default_2MZeroScan.mat -singlematrix
MRIParam='default paramA paramB'
for MRIParam in ${MRIParam}
do 
echo "Running applyxfm4D for: ${subj} ${MRIParam} ..."
applyxfm4d ${srcout}/StatsOutput/Zstat2_percentage/Zstat2_concat_${subj}_${MRIParam}_percentage.nii.gz ${srcout}/${subj}/${subj}_MZeroScan_brain.nii.gz ${srcout}/${subj}/${subj}_${MRIParam}_analysis/Zstat2_concat_${subj}_${MRIParam}_MZeroScan.nii.gz ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_${MRIParam}_2MZeroScan.mat -singlematrix
#echo "applyxfm4d ${srcout}/StatsOutput/Zstat2_percentage/Zstat2_concat_${subj}_${MRIParam}_percentage.nii.gz ${srcout}/${subj}/${subj}_MZeroScan_brain.nii.gz ${srcout}/${subj}/${subj}_${MRIParam}_analysis/Zstat2maps_concat_${subj}_${MRIParam}_MZeroScan.nii.gz ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_${MRIParam}_2MZeroScan.mat -singlematrix
done


echo "~ End of do_analysis_applyxfm4d.sh ~"#${srcout}/StatsOutput/Zstat2maps_percentage/Zstat2maps_concat_${subj}_${MRIParam}_percentage.nii.gz
#${srcout}/${subj}/${subj}_MZeroScan

