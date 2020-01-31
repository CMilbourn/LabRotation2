%% README asl_analysis

% function T=asl_analysis(src,subj)

[cvr_default, dims, scales, bpp, endian]= read_avw([src subj '/func/' subj '_cvr_default']);
[cvr_paramA, dims, scales, bpp, endian] = read_avw([src subj '/func/' subj '_cvr_paramA']);
[cvr_paramB, dims, scales, bpp, endian] = read_avw([src subj '/func/' subj '_cvr_paramB']);

%% *read_avw* 
%% 
% * comes with fsl 
%% Builds up the directory for source subj number
% need to update these to derivatives
%fills in the variables - 
%*Dataset* (CVR_MRI parameter)
%*dims* = dimension of the matrix – 64 rows, cols and so many slices , and 
% so many volumes 
%*scales* = list of dimesons e.g. 3.5x3.5.4 x time
%*bpp* = byte per pixel – so that if you can save out with same ‘byte’ 
% *endian* = how numbers stored - determine the number, some read from big 
% end some read second then first big endian.
% Path read in = sourcedata folder/subjectnumber/func/subject_cvr_MRIparameter
% *Inputs*: sourcedata folder/subjectnumber/func/subject_cvr_MRIparameter
% *Outputs*: Loads in default, paramA & paramB cvr maps and assigns variables 

% gm_cortex %
gm_cortex=read_avw([src subj '/func/' subj '_GM_MZeroScan']);
%% Reads in the cortex map that we created GM_MZeroScan and assigns it to variable 
% called gm_cortex
% *Input*: sourcedata folder/subjectnumber/func/subject_ GM_MZeroScan
% *Output*: variable gm_cortex
%%figure 

figure; %creates a figure
imagesc(cvr_default(:,:,6),[0 2]); %displays cvr_default as a sanity check
colorbar; %colour bar displayed in figure
%% reshape
cvr_defaultr=reshape(cvr_default,64*64*12,1);
cvr_paramAr=reshape(cvr_paramA,64*64*12,1);
cvr_paramBr=reshape(cvr_paramB,64*64*12,1);
gm_cortexr=reshape(gm_cortex,64*64*12,1);;

% *Cvr_defaultr=reshape*  - takes cvr default map  put an r on the end of the 
% variable for a reshaped map'  reshape pass through matrix and tell it how 
% to look like instead of 4D it goes down to 2 each column is a time point 
% The last variable on the thing
% *Input*: cvr_MRIparameter, dimensions of the file
% *Output*: reshaped file called ‘cvr_MRI parameterr

%% Make Mask for gm_cortex
cvr_default_vals=cvr_defaultr(gm_cortexr>0,:);
cvr_paramA_vals=cvr_paramAr(gm_cortexr>0,:);
cvr_paramB_vals=cvr_paramBr(gm_cortexr>0,:); 
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

cvrmeans(1,:)=mean(cvr_default_vals);
%% 
% Cvr_default _vals ==> so gm_cx >>0 col of 1 and 0s  so those w. > 0 
% •	> cvr_default_vals is acol of numbers picked out pf cvrmap and from there 
% do stuff to is
% Cvrmeans - numbers into row 1
% Cvr_default_vals will only be one number – number that have been picked 
% out to give column of one number
% so it could just say row 1 is the default values
% Row 2 is the parameterA values
% row 3 is he parameter B values 
% - made it more confusing by saying that calc the mean and medians in difference 
% ways
% Cvrmeans takes means of num you get 

cvrmeansnoZ(1,:)=mean(cvr_default_vals(cvr_default_vals~=0));
%% Create Mean and Median variables for table
% MeansnoZ removes any numbers that are 0 -> so only actually calculating 
% stuff w.in the mask 
% so look at cvr_default vales – check any are not == to 0 
cvrmedians(1,:)=median(cvr_default_vals);
cvrmediansnoZ(1,:)=median(cvr_default_vals(cvr_default_vals~=0));
cvrmeans(2,:)=mean(cvr_paramA_vals);
cvrmeansnoZ (2,:)=mean(cvr_paramA_vals(cvr_paramA_vals~=0));
cvrmedians(2,:)=median(cvr_paramA_vals);
cvrmediansnoZ(2,:)=median(cvr_paramA_vals(cvr_paramA_vals~=0));
cvrmeans(3,:)=mean(cvr_paramB_vals);
cvrmeansnoZ(3,:)=mean(cvr_paramB_vals(cvr_paramB_vals~=0));
cvrmedians(3,:)=median(cvr_paramB_vals);
cvrmediansnoZ(3,:)=median(cvr_paramB_vals(cvr_paramB_vals~=0));

% Does same for medians And mediansnoZ 
%For all 3 parameters
%% Make Table
T=table(rows,cvrmeans,cvrmeansnoZ,cvrmedians,cvrmediansnoZ);
% T=table names of rows variables of the columns names - > in command window 
% you have a table w. headings of rows, cvrm, … and numbers listed down the columns 
% - that table T is produced by the funcn 

%% Make Histograms %%
figure;
hist(cvr_default_vals,1000)
title('CVR default');

figure;
hist(cvr_paramA_vals,1000)
figure;
hist(cvr_paramB_vals,1000)
%% 
% *hist ->* Put it in that histogram 
% - do one for each 
% - so you can look at it and see if mean is better representation to the peak 
% or the median 0 
% Often the median is better, in some cases – when less noisy data then mean 
% and median is sim 
% If noisy then mean and median will be far apart but median more likely to 
% be close to the peak