%% README asl_analysis
 
%% Line to make it into a function
% function T=asl_analysis(src,subj)
%% Housekeeping
clear all % clear workspace
close all % close any open figures
clc %clear command window


% %% Set Variables
% src = '/Users/colette/sourcedata/' %set source as the folder to sourcedata
% subj = 'sub-01' %set subject to the subject you wish to run in format sub-01
% srcout = [src '/derivatives/' subj]
% MRIParam = 'default' %set MRI Parameter to default paramA paramB

% addpath = '/usr/local/fsl/etc/fslconf/'
% addpath = '/usr/local/fsl/etc/fslconf/fsl.sh'
% '${FSLDIR}' = '/usr/local/'
%% set FSL environment 
fprintf 'Setting up FSL environment...'
setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is 
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be 

%% *read_avw* 
 
cvr_default = '/Users/colette/sourcedata/derivatives/StatsOutput/CVRmaps_all/CVRmaps_concat_sub-01_default.nii.gz'
[cvr_default, dims, scales, bpp, endian] = read_avw([cvr_default]);


%
%/Users/colette/sourcedata/derivatives/StatsOutput/sub-01/sub-01_default/stats_sub-01_default_060.nii.gz


cvr_paramA = '/Users/colette/sourcedata/derivatives/StatsOutput/CVRmaps_all/CVRmaps_concat_sub-01_paramA.nii.gz'
[cvr_paramA, dims, scales, bpp, endian] = read_avw([cvr_paramA]);
cvr_paramB = '/Users/colette/sourcedata/derivatives/StatsOutput/CVRmaps_all/CVRmaps_concat_sub-01_paramB.nii.gz'
[cvr_paramB, dims, scales, bpp, endian] = read_avw([cvr_paramB]);

% * comes with fsl 
% Builds up the directory for source subj number %%
% need to update these to derivatives
%fills in the variables - 
%*Dataset* (CVR_MRI parameter)
%*dims* = dimension of the matrix: 64 rows, 64 cols, 12 slices, 1 volumes 
%*scales* = list of dimesons e.g. 3.5x3.5.4 x time
%*bpp* = byte per pixel - so that if you can save out with same bytes
% *endian* = how numbers stored - determine the number, some read from big 
% end some read second then first big endian.
% Path read in = sourcedata folder/subjectnumber/func/subject_cvr_MRIparameter
% *Inputs*: sourcedata folder/subjectnumber/func/subject_cvr_MRIparameter
% *Outputs*: Loads in default, paramA & paramB cvr maps and assigns variables 

% gm_cortex %
gm_cortex=read_avw(['/Users/colette/sourcedata/derivatives/sub-01/sub-01_GM_MZeroScan.nii.gz']);
%% Reads in the cortex map that we created GM_MZeroScan and assigns it to variable 
% called gm_cortex
% *Input*: sourcedata folder/subjectnumber/func/subject_ GM_MZeroScan
% *Output*: variable gm_cortex
%%figure 

figure; %creates a figure
imagesc(cvr_default(:,:,6),[0 2]); %displays cvr_default as a sanity check
colorbar; %colour bar displayed in figure

%% reshape
default_cvr_r=reshape(cvr_default,64*64*12,31);
paramA_cvr_r=reshape(cvr_paramA,64*64*12,31);
paramB_cvr_r=reshape(cvr_paramB,64*64*12,31);
gm_cortexr=reshape(gm_cortex,64*64*12,1);

% *Cvr_defaultr=reshape*  - takes cvr default map  put an r on the end of the 
% variable for a reshaped map'  reshape pass through matrix and tell it how 
% to look like instead of 4D it goes down to 2 each column is a time point 
% The last variable on the thing
% *Input*: cvr_MRIparameter, dimensions of the file
% *Output*: reshaped file called ‘cvr_MRI parameterr

%% Make Mask for gm_cortex
default_cvr_vals=default_cvr_r(gm_cortexr>0,:);
paramA_cvr_vals=paramA_cvr_r(gm_cortexr>0,:);
paramB_cvr_vals=paramB_cvr_r(gm_cortexr>0,:); 
% Gone from matrix with a stack  but also in time (not in this case) to each 
% time point you have a column
% Time is subsequent columns so that if you have  a cvr map and another column 
% w. GM mask 
% For the GM column you can say GM >0  ===> another col - > has a 1 or 0 for 
% voxel w. any datapoint >>0 so it creates a mask 
% Then takes cvr map and put the mask into the argument it will give values 
% for each point that is picked out and then from there calc the median of the 
% cvrmask 
% 
% If you have more than one column then change it to go to CVR(mask, ;) so it 
% takes all the things instead of just one column of numbers have a column for 
% even time point and median down the column instead of across their rows
% 
% When you do reshape the quickest way is to take the param which has multiple 
% dimensions – have   param=pparam(:) 
% If you have more than one column then you need to change the number on the 
% end to have the number of time points 
%  you could pull these from dimes - > 64 dims 1  12  dims3, 
%
% Could end up with 30 time points - > have that in the dimensions 
% *gm_cortex* will only have one time point 
%% Set up table inputs

rows={'default'; 'paramA'; 'paramB'};

cvrmeans(1,:)=mean(default_cvr_vals,2); %default as row, added comma 2  
%% 
% Cvr_default _vals ==> so gm_cx >>0 col of 1 and 0s  so those w. > 0 
% > default_cvr_vals is acol of numbers picked out pf cvrmap and from there 
% do stuff to is
% Cvrmeans - numbers into row 1
% default_cvr_vals will only be one number – number that have been picked 
% out to give column of one number
% so it could just say row 1 is the default values
% Row 2 is the parameterA values
% row 3 is he parameter B values 
% - made it more confusing by saying that calc the mean and medians in difference 
% ways
% Cvrmeans takes means of num you get 

cvrmeansnoZ(1,:)=mean(default_cvr_vals(:,default_cvr_vals~=0)); %not =/= to 0 - find for 0 andonly applying to one time point
% 


%% Create Mean and Median variables for table
% MeansnoZ removes any numbers that are 0 -> so only actually calculating 
% stuff w.in the mask 
% so look at cvr_default values - matrix will have different size columns
%check any are not == to 0 - not able to have empty emelemnts 
% loop through the columns -> calc the median for each column hek
% do for including and not including 0's
% do same for zstat to see if a peak in the 
%see if it is tag and ctrl diffr
% work out peak in zstat and cvr maps 

cvrmedians(1,:)=median(default_cvr_vals,2);
cvrmediansnoZ(1,:)=median(default_cvr_vals(default_cvr_vals~=0),2);
cvrmeans(2,:)=mean(paramA_cvr_vals,2);
cvrmeansnoZ (2,:)=mean(paramA_cvr_vals(paramA_cvr_vals~=0),2);
cvrmedians(2,:)=median(paramA_cvr_vals,2);
cvrmediansnoZ(2,:)=median(paramA_cvr_vals(paramA_cvr_vals~=0),2);
cvrmeans(3,:)=mean(paramB_cvr_vals,2);
cvrmeansnoZ(3,:)=mean(paramB_cvr_vals(paramB_cvr_vals~=0),2);
cvrmedians(3,:)=median(paramB_cvr_vals,2);
cvrmediansnoZ(3,:)=median(paramB_cvr_vals(paramB_cvr_vals~=0),2);
%1x31 or 31x1 or by 3 - one for each MRIparam??

% for SD - weighted 
% Does same for medians And mediansnoZ 
%For all 3 parameters
%% Make Table
T=table(rows,cvrmeans,cvrmeansnoZ,cvrmedians,cvrmediansnoZ); 
% T=table names of rows variables of the columns names - > in command window 
% you have a table w. headings of rows, cvrm, and numbers listed down the columns 
% - that table T is produced by the funcn 
% should have 32 columns - could look at cvrmediansnoZ - plot cvrmedianz
%% Make Histograms %%
figure;
hist(default_cvr_vals(:,16),1000) %distribution for each time delay - 
%hist(default_cvr_vals,1000) %1000 bins - can specify what bins, could
%specify range of valsues e..g 1:1:10 to do in steps of 10 seconds - inside
%bins - then outside bins with everythign else - look
%
title('CVR default');
fprintf 'default'

figure;
hist(paramA_cvr_vals,1000)
fprintf 'paramA'

figure;
hist(paramB_cvr_vals,1000)
fprintf 'paramB'
%% 
% *hist ->* Put it in that histogram 
% - do one for each 
% - so you can look at it and see if mean is better representation to the peak 
% or the median 0 
% Often the median is better, in some cases – when less noisy data then mean 
% and median is sim 
% If noisy then mean and median will be far apart but median more likely to 
% be close to the peak

%% run through index of skip number time point
% -60+(index-1)*4