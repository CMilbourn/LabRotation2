# !/bin/bash
# do_analysis.sh sourcedata subj_no
#Prints out usage of the code and input parameters
if [ $# -lt 1 ]; then
  echo "Usage: $0 <sourcedata folder> <SubjectNumber> e.g. ./do_analysis_edit_2.sh /Users/colette/sourcedata/ sub-01"; 
 #exit 1 ;
fi
# 
src=$1
subj=$2
PVthresh=0.5
srcout=/Users/colette/sourcedata/derivatives/
MRIParam='default paramA paramB'
## Make CVR maps ##
## Make directory for Percentage CVRmaps ##
if ! [ -d "${srcout}/StatsOutput/CVRmaps_percentage" ]; then #check if the output directory exists, if not make it
mkdir ${srcout}/StatsOutput/CVRmaps_percentage
echo "Making CVRmaps_all folder:${srcout}/StatsOutput/CVRmaps_percentage"
fi 

echo "Making CVRmaps percentatges for $2..."

## Loop to make CVRmaps percentages
MRIParam='default paramA paramB'
for MRIParam in ${MRIParam}
do 
${FSLDIR}/bin/fslmaths ${srcout}/StatsOutput/CVRmaps_all/CVRmaps_concat_${subj}_${MRIParam}.nii.gz -mul 100 ${srcout}/StatsOutput/CVRmaps_percentage/CVRmaps_concact_${subj}_${MRIParam}_percentage.nii.gz
echo "... CVRmaps_concact_${subj}_${MRIParam}_percentage.nii.gz"
done 
 
## Anatomical analysis ##
echo "Running anatomical analysis for $2..."
echo "Running BET..."
${FSLDIR}/bin/bet ${src}/${subj}/anat/${subj}_T1w ${srcout}/${subj}/${subj}_T1w_brain -f 0.3 -g 0

echo "Running Fast for $2...only run once then comment out"
${FSLDIR}/bin/fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o ${srcout}/${subj}/${subj}_T1w_brain ${srcout}/${subj}/${subj}_T1w_brain


## Make reference images for registrationfunc ##
echo "Making reference images"

echo "Making reference images for $2 Mzeroscan..." 
${FSLDIR}/bin/mcflirt -in ${src}/${subj}/func/${subj}_MZeroScan -out ${srcout}/${subj}/${subj}_MZeroScan_mc 
${FSLDIR}/bin/fslmaths ${srcout}/${subj}/${subj}_MZeroScan_mc -Tmean ${srcout}/${subj}/${subj}_MZeroScan_ref
${FSLDIR}/bin/bet ${srcout}/${subj}/${subj}_MZeroScan_ref ${srcout}/${subj}/${subj}_MZeroScan_brain -Z

MRIParam='default paramA paramB'
for MRIParam in ${MRIParam}
do 

if ! [ -d "${srcout}/${subj}/${subj}_${MRIParam}_analysis" ]; then #check if the output directory exists, if not make it
mkdir ${srcout}/${subj}/${subj}_${MRIParam}_analysis
echo "Making folder: ${srcout}/${subj}/${subj}_${MRIParam}_analysis"
fi 

echo "MCFLIRT for $2 ${MRIParam}..."
${FSLDIR}/bin/mcflirt -in ${src}/${subj}/func/${subj}_task-hyper_acq-${MRIParam}_asl -out ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_task-hyper_acq-${MRIParam}_asl_mc
echo "FSLmaths Tmean for $2 ${MRIParam}..."
${FSLDIR}/bin/fslmaths ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_task-hyper_acq-${MRIParam}_asl_mc -Tmean ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_task-hyper_acq-${MRIParam}_asl_ref
echo "BET for $2 ${MRIParam}..."
${FSLDIR}/bin/bet ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_task-hyper_acq-${MRIParam}_asl_ref ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_task-hyper_acq-${MRIParam}_asl_brain -Z
done  
 
## Registration of functional and anatomical images ##
echo "Registering functional and anatomical images..."

## T1w to MNI ##
echo "T1w to MNI..."
${FSLDIR}/bin/flirt -in ${srcout}/${subj}/${subj}_T1w_brain -ref ${FSLDIR}/data/standard/MNI152_T1_1mm_brain -out ${srcout}/${subj}/${subj}_T1w_MNI -omat ${srcout}/${subj}/${subj}_T1w2MNI.mat -v


## MZeroScan to T1w ##
echo "MZeroScan to T1w..."
${FSLDIR}/bin/flirt -in ${srcout}/${subj}/${subj}_MZeroScan_brain -ref ${srcout}/${subj}/${subj}_T1w_brain -out ${srcout}/${subj}/${subj}_MZeroScan_T1w -omat ${srcout}/${subj}/${subj}_MZeroScan2T1w.mat -dof 6 -v

## PCASL to MZeroScan ##
echo "PCASL to MZeroScan..."

MRIParam='default paramA paramB'
for MRIParam in ${MRIParam}
do 
# old version - ${FSLDIR}/bin/flirt -in ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_task-hyper_acq-${MRIParam}_asl_brain -ref ${srcout}/${subj}/${subj}_MZeroScan_brain -out ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_task-hyper_acq-${MRIParam}_asl_MZeroScan -dof 6 -omat ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_${MRIParam}_2MZeroScan.mat -v 
# from terminal /usr/local/fsl/bin/flirt -in /Users/colette/sourcedata/derivatives/sub-01/sub-01_paramA_analysis/sub-01_task-hyper_acq-paramA_asl_brain.nii.gz -ref /Users/colette/sourcedata/derivatives/sub-01/sub-01_MZeroScan_brain.nii.gz -out /Users/colette/sourcedata/derivatives/sub-01/sub-01_paramA_analysis/test.nii.gz -omat /Users/colette/sourcedata/derivatives/sub-01/sub-01_paramA_analysis/test.mat -bins 256 -cost corratio -searchrx 0 0 -searchry 0 0 -searchrz 0 0 -dof 6  -interp trilinear

${FSLDIR}/bin/flirt -in ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_task-hyper_acq-${MRIParam}_asl_brain -ref ${srcout}/${subj}/${subj}_MZeroScan_brain -out ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_task-hyper_acq-${MRIParam}_asl_MZeroScan -omat ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_${MRIParam}_2MZeroScan.mat -bins 256 -cost corratio -searchrx 0 0 -searchry 0 0 -searchrz 0 0 -dof 6  -interp trilinear

#echo "Creating ... ${subj}_task-hyper_acq-${MRIParam}_asl_brain, ${subj}_MZeroScan_brain, ${subj}_task-hyper_acq-${MRIParam}_asl_MZeroScan, ${subj}_${MRIParam}_2MZeroScan.mat"

done 

## Create transform from MNI space to MZeroScan space ##
echo "Create transform from MNI space to MZeroScan space..."
${FSLDIR}/bin/convert_xfm -omat ${srcout}/${subj}/${subj}_T1w2MZeroScan.mat -inverse ${srcout}/${subj}//${subj}_MZeroScan2T1w.mat

echo "T1w2MNI ==> MNI2T1w..."
${FSLDIR}/bin/convert_xfm -omat ${srcout}/${subj}/${subj}_MNI2T1w.mat -inverse ${srcout}/${subj}/${subj}_T1w2MNI.mat

echo "T1w2MZeroScan + MNI2T1w ==> MNI2MZeroScan..."
${FSLDIR}/bin/convert_xfm -omat ${srcout}/${subj}/${subj}_MNI2MZeroScan.mat -concat ${srcout}/${subj}/${subj}_T1w2MZeroScan.mat ${srcout}/${subj}/${subj}_MNI2T1w.mat


## Apply transformations ##
echo "Apply Transformations..."
echo "Creating ${subj}_GM_MZeroScan..."
${FSLDIR}/bin/applywarp --in=${srcout}/${subj}/${subj}_T1w_brain_pve_1 --ref=${srcout}/${subj}/${subj}_MZeroScan_brain --out=${srcout}/${subj}/${subj}_GM_MZeroScan --premat=${srcout}/${subj}/${subj}_T1w2MZeroScan.mat --interp=trilinear --super --superlevel=4
echo "Creating ${subj}_MNIatlas..."
${FSLDIR}/bin/flirt -in ${FSLDIR}/data/atlases/MNI/MNI-maxprob-thr0-1mm -ref ${srcout}/${subj}/${subj}_MZeroScan_brain -out ${srcout}/${subj}/${subj}_MNIatlas -init ${srcout}/${subj}/${subj}_MNI2MZeroScan.mat -interp nearestneighbour -applyxfm -v
  
echo "Applying Transforms "

MRIParam='default paramA paramB'
for MRIParam in ${MRIParam}
do 
${FSLDIR}/bin/flirt -in ${srcout}/StatsOutput/CVRmaps_percentage/CVRmaps_concat_${subj}_${MRIParam}_percentage.nii.gz -ref ${srcout}/${subj}/${subj}_MZeroScan_brain.nii.gz -out ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_${MRIParam}_cvr_MZeroScan -init ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_${MRIParam}_2MZeroScan.mat -applyxfm -v
echo "CVRmaps_concat_${subj}_${MRIParam}_percentage.nii.gz + ${subj}_${MRIParam}_2MZeroScan ==> ${subj}_${MRIParam}_cvr_MZeroScan..."
done 

echo "Make ROIs"

## Make ROIs ##
${FSLDIR}/bin/fslmaths ${srcout}/${subj}/${subj}_MNIatlas -thr 3 -uthr 6 ${srcout}/${subj}/${subj}_MNIatlas_cortex
${FSLDIR}/bin/fslmaths ${srcout}/${subj}/${subj}_MNIatlas -thr 8 -uthr 8 -add ${srcout}/${subj}/${subj}_MNIatlas_cortex ${srcout}/${subj}/${subj}_MNIatlas_cortex
echo "Binarise MNIatlas ..."
${FSLDIR}/bin/fslmaths ${srcout}/${subj}/${subj}_MNIatlas -bin ${srcout}/${subj}/${subj}_MNIatlas_cortex
echo "gmcortex ${subj}_GM_MZeroScan ${subj}_MNIatlas_cortex"
${FSLDIR}/bin/fslmaths ${srcout}/${subj}/${subj}_GM_MZeroScan -thr ${PVthresh} -mas ${srcout}/${subj}/${subj}_MNIatlas_cortex -bin ${srcout}/${subj}/${subj}_gmcortex

## ApplyxfmD4 ##
MRIParam='default paramA paramB'
for MRIParam in ${MRIParam}
do 
echo "Running applyxfm4D for: ${subj} ${MRIParam} ..."
applyxfm4d ${srcout}/StatsOutput/CVRmaps_percentage/CVRmaps_concact_${subj}_${MRIParam}_percentage.nii.gz ${srcout}/${subj}/${subj}_MZeroScan_brain.nii.gz ${srcout}/${subj}/${subj}_${MRIParam}_analysis/CVRmaps_concact_${subj}_${MRIParam}_MZeroScan.nii.gz ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_${MRIParam}_2MZeroScan.mat -singlematrix
#echo "applyxfm4d ${srcout}/StatsOutput/CVRmaps_percentage/CVRmaps_concact_${subj}_${MRIParam}_percentage.nii.gz ${srcout}/${subj}/${subj}_MZeroScan_brain.nii.gz ${srcout}/${subj}/${subj}_${MRIParam}_analysis/CVRmaps_concact_${subj}_${MRIParam}_MZeroScan.nii.gz ${srcout}/${subj}/${subj}_${MRIParam}_analysis/${subj}_${MRIParam}_2MZeroScan.mat -singlematrix
done

echo "~~~ End of do_analysis for $2 ~~~"