# !/bin/bash
# do_analysis_applyxfm4d.sh sourcedata subj_no
#Prints out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0 <sourcedata folder> <SubjectNumber> e.g. ./do_analysis_applyxfm4d.sh /Users/colette/sourcedata/ sub-01"; 
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
#${srcout}/StatsOutput/CVRmaps_all/CVRmaps_concat_${subj}_${MRIParam}.nii.gz -mul 100 ${srcout}/StatsOutput/CVRmaps_percentage/CVRmaps_concact_${subj}_${MRIParam}_percentage.nii.gz


#applyxfm4d CVRmaps_concact_sub-01_default_percentage.nii.gz ../../sub-01/sub-01_MZeroScan_brain.nii.gz CVRmaps_concact_sub-01_default_MZeroScan.nii.gz ../../sub-01/sub-01_default_analysis/sub-01_default_2MZeroScan.mat -singlematrix
MRIParam='default paramA paramB'
for MRIParam in ${MRIParam}
do 
echo "Running applyxfm4D for: ${subj} ${MRIParam} ..."
applyxfm4d ${srcout}/StatsOutput/CVRmaps_percentage/CVRmaps_concact_${subj}_${MRIParam}_percentage.nii.gz ${srcout}/${subj}/${subj}_MZeroScan_brain.nii.gz ${srcout}/${subj}/${subj}_${MRIParam}_analysis/CVRmaps_concact_${subj}_${MRIParam}_MZeroScan.nii.gz ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_${MRIParam}_2MZeroScan.mat -singlematrix
#echo "applyxfm4d ${srcout}/StatsOutput/CVRmaps_percentage/CVRmaps_concact_${subj}_${MRIParam}_percentage.nii.gz ${srcout}/${subj}/${subj}_MZeroScan_brain.nii.gz ${srcout}/${subj}/${subj}_${MRIParam}_analysis/CVRmaps_concact_${subj}_${MRIParam}_MZeroScan.nii.gz ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_${MRIParam}_2MZeroScan.mat -singlematrix
done


echo "~ End of do_analysis_applyxfm4d.sh ~"#${srcout}/StatsOutput/CVRmaps_percentage/CVRmaps_concact_${subj}_${MRIParam}_percentage.nii.gz
#${srcout}/${subj}/${subj}_MZeroScan

