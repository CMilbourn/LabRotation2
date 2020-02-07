%% README asl_analysis
 
%% Line to make it into a function
% function T=asl_analysis(src,subj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all % clear workspace
close all % close any open figures
clc %clear command window
diary on %turn on diary
diary my_data.out %ouput the information into a file called 'my_data.out'
echo on %echo is turned on - this will print most things into the command window to turn of use 'echo off'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set FSL environment 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf 'Setting up FSL environment...' %prints to command window 'Setting up FSL environment...'
setenv('FSLDIR','/usr/local/fsl');  % path to FSL folder
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % Set up output type 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source In - Path to datafolder
src = '/Users/colette/sourcedata/'
%Source Out - Path to previously processed data
srcout = '/Users/colette/sourcedata/derivatives/StatsOutput/'

%Check if the folder 'DataOutput' exists, if it does not then creates it
    if ~exist(sprintf('%s/derivatives/DataOutput',src), 'dir')
       mkdir(sprintf('%s/derivatives/DataOutput',src))
    end
DataOutput = sprintf('%s/derivatives/DataOutput',src) %Assign variable DataOutput to the folder


%Check if the subfolder 'PerSubj' exists, if it does not then creates it
    if ~exist(sprintf('%s/derivatives/DataOutput/PerSubj',src), 'dir')
       mkdir(sprintf('%s/derivatives/DataOutput/PerSubj',src))
    end
DataOutputPerSubj = sprintf('%s/derivatives/DataOutput/PerSubj',src) %Assign variable DataOutput to the folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For Loop to create files of mean, median Zstats, both with & without 0's %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(i.e. without 0's excludes areas outside the Grey matter mask) 
%for k=1:3 %to run through subjects 1 to 3
numlist = {1,2,3,5,7,8,9}; %Creates a list with the subject numbers to run through called 'numlist'
for k = 1:length(numlist) %start of for loop - runs through one number in 'numlist' per iteration
    %% Variables
    subj= sprintf('%02d',numlist{k}) %assigns the variable 'subj' the number from numlist as a string in 01 format
    src = '/Users/colette/sourcedata/' %source in
    srcout2 = sprintf('%sderivatives/sub-%s/sub-%s_', src, subj, subj) %alternative source out if it is not the same as srcout

    %% *read_avw*
    %Reads in the zstat files created by 'do_analysis_zstat.sh'
    %e.g. file in: /Users/colette/sourcedata/derivatives/sub-01/sub-01_default_analysis/Zstat2_concat_sub-01_default_MZeroScan.nii.gz 
    MRIParam1= 'default'
    [default_zstat, dims, scales, bpp, endian] = read_avw([srcout2 MRIParam1 '_analysis/Zstat2_concat_sub-' subj '_' MRIParam1 '_MZeroScan.nii.gz']);
    
    MRIParam2 = 'paramA'
    [paramA_zstat, dims, scales, bpp, endian] = read_avw([srcout2 MRIParam2 '_analysis/Zstat2_concat_sub-' subj '_' MRIParam2 '_MZeroScan.nii.gz']);
    
    MRIParam3 = 'paramB'
    [paramB_zstat, dims, scales, bpp, endian] = read_avw([srcout2 MRIParam3 '_analysis/Zstat2_concat_sub-' subj '_' MRIParam3 '_MZeroScan.nii.gz']);
   
    %% gm_cortex %
    %Reads in the grey matter cortex mask created by 'do_analysis_edit2.sh'
    %e.g. file: /Users/colette/sourcedata/derivatives/sub-01/sub-01_gmcortex.nii.gz
    gm_cortex = read_avw([sprintf('%sderivatives/sub-%s/sub-%s_gmcortex.nii.gz', src, subj, subj)]);
    % *Output*: variable gm_cortex
    
    %% figure for checking image input
%     figure; %creates a figure for one slice/volume 
%     imagesc(default_zstat(:,:,6),[0 2]); %displays default_zstat as a sanity check
%     colorbar; %colour bar displayed in figure
    
    %% reshape
    %Reshapes the 4D singles into 419152x31 single - appends '_r' to new
    %variable
    default_zstat_r=reshape(default_zstat,64*64*12,31); %64x64x12 matrix, 31 timepoints(skip numbers)
    paramA_zstat_r=reshape(paramA_zstat,64*64*12,31);
    paramB_zstat_r=reshape(paramB_zstat,64*64*12,31);
    gm_cortex_r=reshape(gm_cortex,64*64*12,1);
    
    %% Make Mask for gm_cortex
    default_zstat_vals=default_zstat_r(gm_cortex_r>0,:);
    paramA_zstat_vals=paramA_zstat_r(gm_cortex_r>0,:);
    paramB_zstat_vals=paramB_zstat_r(gm_cortex_r>0,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    %% Create Mean and Median variables for table %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Make Medians without 0's %%
    default_zstat_valsNOZ=default_zstat_vals;
    default_zstat_valsNOZ(default_zstat_valsNOZ(:,1)==0,:)=NaN;
    zstatmedians_NOZ(1,:)= nanmedian(default_zstat_valsNOZ)
    
    zstatmedians_NOZ(1,:)= nanmedian(default_zstat_valsNOZ)
    
    paramA_zstat_valsNOZ=paramA_zstat_vals;
    paramA_zstat_valsNOZ(paramA_zstat_valsNOZ(:,1)==0,:)=NaN;
    zstatmedians_NOZ(2,:)=nanmedian(paramA_zstat_valsNOZ)
    
    paramB_zstat_valsNOZ=paramB_zstat_vals;
    paramB_zstat_valsNOZ(paramB_zstat_valsNOZ(:,1)==0,:)=NaN;
    zstatmedians_NOZ(3,:)=nanmedian(paramB_zstat_valsNOZ)
    
    %save zstatmedians to one variable column per iteration
    zstatmediansNOZ_final{k}=(zstatmedians_NOZ)'
    dlmwrite(sprintf('%s/zstatmediansNOZ_final_sub-%s_%s.tsv', DataOutputPerSubj, subj, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmediansNOZ_final{k},'\t')
    
    %% Medians WITH 0's  %%
    zstatmedians(1,:)=median(default_zstat_vals,1);
    zstatmedians(2,:)=median(paramA_zstat_vals,1);
    zstatmedians(3,:)=median(paramB_zstat_vals,1);
    
    %save to one part per iteration
    zstatmedians_final{k}=(zstatmedians)'
    dlmwrite(sprintf('%s/zstatmedians_final_sub-%s_%s.tsv', DataOutputPerSubj, subj, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmedians_final{k},'\t')
    
    %% Make Means without 0's %%
    default_zstat_valsNOZ=default_zstat_vals; %Create variable for zstat_vals with NO Zeros
    default_zstat_valsNOZ(default_zstat_valsNOZ(:,1)==0,:)=NaN; %replace any 0' with NAN 
    zstatmeans_NOZ(1,:)= nanmean(default_zstat_valsNOZ) %Mean of zstat_vals without any 0's(areas outside the mask)
    
    paramA_zstat_valsNOZ=paramA_zstat_vals; %Create variable for zstat_vals with NO Zeros
    paramA_zstat_valsNOZ(paramA_zstat_valsNOZ(:,1)==0,:)=NaN;%replace any 0' with NAN
    zstatmeans_NOZ(2,:)=nanmean(paramA_zstat_valsNOZ) %Mean of zstat_vals without any 0's(areas outside the mask)
    
    paramB_zstat_valsNOZ=paramB_zstat_vals; %Create variable for zstat_vals with NO Zeros
    paramB_zstat_valsNOZ(paramB_zstat_valsNOZ(:,1)==0,:)=NaN; %replace any 0' with NAN
    zstatmeans_NOZ(3,:)=nanmean(paramB_zstat_valsNOZ) %Mean of zstat_vals without any 0's(areas outside the mask)
    
    %save to zstatmeansNOZ for each iteration
    zstatmeansNOZ_final{k}=(zstatmeans_NOZ)' %transpose 
    dlmwrite(sprintf('%s/zstatmeansNOZ_final_sub-%s_%s.tsv', DataOutputPerSubj, subj, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmeansNOZ_final{k},'\t')  
    
    %% Means WITH 0's %%
    zstatmeans(1,:)=mean(default_zstat_vals,1);
    zstatmeans(2,:)=mean(paramA_zstat_vals,1);
    zstatmeans(3,:)=mean(paramB_zstat_vals,1);
    
    zstatmeans_final{k}=(zstatmeans)'
    dlmwrite(sprintf('%s/zstatmeans_final_sub-%s_%s.tsv', DataOutputPerSubj, subj, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmeans_final{k},'\t')
    
    %% Save zstat_final files %%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Save total zstat files as .tsv files %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%zstatmedians
dlmwrite(sprintf('%s/zstatmedians_final_all_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmedians_final,'\t')
%zstatmediansNOZ
dlmwrite(sprintf('%s/zstatmediansNOZ_final_all_%s.tsv', DataOutput,datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmediansNOZ_final,'\t')
%zstatmeans
dlmwrite(sprintf('%s/zstatmeans_final_all_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmeans_final,'\t')
%zstatmeansNOZ
dlmwrite(sprintf('%s/zstatmeansNOZ_final_all_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmeansNOZ_final,'\t')  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up table inputs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%rows - skip numbers -60 secs to +60 secs in intervals of 4secs
%row={-60:4:60}';
%row=(-60:4:60)';
% Headers %
%Each parameter is entered mannually cos it does not like the loop
%Headers = {'SubjectNo_MRIParam','sub01default','sub01paramA','sub01paramB',...
Headers = {'sub01default','sub01paramA','sub01paramB',...
'sub02default','sub02paramA','sub02paramB',...
'sub03default','sub03paramA','sub03paramB',...
'sub05default','sub05paramA','sub05paramB',...
'sub07default','sub07paramA','sub07paramB',...
'sub08default','sub08paramA','sub08paramB'...
'sub09default','sub09paramA','sub09paramB'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TABLES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% zstatmeans %
% Tablezstatmeans = cell2mat(zstatmeans_final) %converts from cell array to matrix
% %Tablezstatmeans2 = [row; num2cell(Tablezstatmeans)]
% Tablezstatmeans2 = [Headers; num2cell(Tablezstatmeans)] %adds headers to columns not sure why it has to be cells 
% Tablezstatmeans3 = cell2table(Tablezstatmeans2) %converts from cell to table format
% Tablezstatmeans4 = Tablezstatmeans2(2:end,:); %assigns from row 2 until end and all columns 
% Table_zstatmeans5 = cell2table(Tablezstatmeans4) %converts from cell to table format
% Table_zstatmeans6 = table(rows,Table_zstatmeans5)
% Table_zstatmeans5.Properties.VariableNames = Tablezstatmeans2(1,:) %adds correct headers in table format
% zstatmeans_final_table = Table_zstatmeans5 %saves as '_final_table'

%sumzstatmeans=summary(zstatmeans_final_table)
%zstatmeans_final_table.Overallmean = mean(zstatmeans_final_table{:,2:end},2)
%zstatmeans_final_table.eachtimepointmean = mean(zstatmeans_final_table{2:end,:})

% zstatmeans %
Tablezstatmeans = cell2mat(zstatmeans_final)
Tablezstatmeans2 = [Headers; num2cell(Tablezstatmeans)]
Tablezstatmeans3 = cell2table(Tablezstatmeans2)
Tablezstatmeans4 = Tablezstatmeans2(2:end,:);
Table_zstatmeans5 = cell2table(Tablezstatmeans4)
Table_zstatmeans5.Properties.VariableNames = Tablezstatmeans2(1,:)
zstatmeans_final_table = Table_zstatmeans5

% zstatmeansNOZ %
TablezstatmeansNOZ = cell2mat(zstatmeansNOZ_final)
TablezstatmeansNOZ2 = [Headers; num2cell(TablezstatmeansNOZ)]
TablezstatmeansNOZ3 = cell2table(TablezstatmeansNOZ2)
TablezstatmeansNOZ4 = TablezstatmeansNOZ2(2:end,:);
Table_zstatmeansNOZ5 = cell2table(TablezstatmeansNOZ4)
Table_zstatmeansNOZ5.Properties.VariableNames = TablezstatmeansNOZ2(1,:)
zstatmeansNOZ_final_table = Table_zstatmeansNOZ5




% zstatmedians %
Tablezstatmedians = cell2mat(zstatmedians_final)
Tablezstatmedians2 = [Headers; num2cell(Tablezstatmedians)]
Tablezstatmedians3 = cell2table(Tablezstatmedians2)
Tablezstatmedians4 = Tablezstatmedians2(2:end,:);
Table_zstatmedians5 = cell2table(Tablezstatmedians4)
Table_zstatmedians5.Properties.VariableNames = Tablezstatmedians2(1,:)
zstatmedians_final_table = Table_zstatmedians5

%zstatmediansNOZ
TablezstatmediansNOZ = cell2mat(zstatmediansNOZ_final)
TablezstatmediansNOZ2 = [Headers; num2cell(TablezstatmediansNOZ)]
TablezstatmediansNOZ3 = cell2table(TablezstatmediansNOZ2)
TablezstatmediansNOZ4 = TablezstatmediansNOZ2(2:end,:);
Table_zstatmediansNOZ5 = cell2table(TablezstatmediansNOZ4)
Table_zstatmediansNOZ5.Properties.VariableNames = TablezstatmediansNOZ2(1,:)
zstatmediansNOZ_final_table = Table_zstatmediansNOZ5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT means and medians with & without 0's %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots Tablezstat's for mean, median both w. & w.o 0's
%adds a title
figure;
plot(Tablezstatmeans);
title('Tablezstatmeans')
figure;
plot(Tablezstatmedians);
title('Tablezstatmedians')
figure;
plot(TablezstatmeansNOZ);
title('TablezstatmeansNOZ')
figure;
plot(TablezstatmediansNOZ);
title('TablezstatmediansNOZ')

% Plots mean of the default, paramA & paramB medians %
figure
plot((-60:4:60),mean(Tablezstatmedians(:,1:3:end),2),'r')
hold on %plots the following two onto the same graph
plot((-60:4:60),mean(Tablezstatmedians(:,2:3:end),2),'g')
plot((-60:4:60),mean(Tablezstatmedians(:,3:3:end),2),'b')
title('Zstat Means of Medians');
grid
legend('Default','paramA','paramB') %default = red, paramA = green, paramB = blue
%this is not RGB friendly so might change the colours later. 
hold off

% Plots mean of the default, paramA & paramB means %
figure
plot((-60:4:60),mean(Tablezstatmeans(:,1:3:end),2),'r')
hold on
plot((-60:4:60),mean(Tablezstatmeans(:,2:3:end),2),'g')
plot((-60:4:60),mean(Tablezstatmeans(:,3:3:end),2),'b')
title('Zstat Means of Means');
grid
legend('Default','paramA','paramB')
hold off

%% Need to get a mean of each column of the new table - excluding the header %%

%dlmwrite(sprintf('%s/zstatmeans_final_all_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),Table_new,'\t')

% save('savefile.mat', 'T')
% save(sprintf('%s/zstatmeans_table_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),Table_new,'\t')
% save('zstatmeans_table.tsv', 'Table_new')
% dlmwrite(sprintf('%s/zstatmeans_table_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),Table_new,'\t') 
% 
% dlmwrite('zstatmeans_table',Table_new,'\t') 


%% Make loop for tables here: 
% zstatlist = [zstatmeans_final zstatmeansNOZ_final zstatmedians_final zstatmedians_final];
% for c= 1:4
%  %final = sprintf('%s_final',c)
%   final= zstatlist({c})
% keyboard;
%     TableTest = cell2mat(final)
% end
%     Headers = {'sub01default','sub01paramA','sub01paramB',...
%         'sub02default','sub02paramA','sub02paramB',...
%         'sub03default','sub03paramA','sub03paramB',...
%         'sub05default','sub05paramA','sub05paramB',...
%         'sub07default','sub07paramA','sub07paramB',...
%         'sub08default','sub08paramA','sub08paramB'...
%         'sub09default','sub09paramA','sub09paramB'};
%     TableTest2 = [Headers; num2cell(TableTest)]
%     TableTest3 = cell2table(TableTest2)
%     TableTest4 = TableTest2(2:end,:);
%     Table_new = cell2table(TableTest4)
%     Table_new.Properties.VariableNames = TableTest2(1,:)
%     sprintf('%s_table',final) = Table_new
% end    
%     
%     
%% run through index of skip number time point
% -60+(index-1)*4
%% End of Script %%
fprintf('~~~ End of Script ~~~'); %prints to command window '~~~ End of Script ~~~'