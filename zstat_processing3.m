%% README asl_analysis
%What does this script do?
%set up FSL environment
%Read in zstatmap percentage data
%reads in GM cortex mask
%creates output folders
%loops through subjects {1,2,3,5,7,8,9,10,11,12,13,14}
%removes 0's
%outputs mean, medians both w. & w.o 0's
%Creates tables with this data in format sub1 default, sub1 paramA, sub1
%Paramb, sub2 default etc
%creates graphs for all subjs, also with mean of mean and medians 
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

%Check if folder exists
    if ~exist(sprintf('%s/derivatives/FiguresOut',src), 'dir')
       mkdir(sprintf('%s/derivatives/FiguresOut',src))
    end
FiguresOut = sprintf('%s/derivatives/FiguresOut',src) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For Loop to create files of mean, median Zstats, both with & without 0's %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(i.e. without 0's excludes areas outside the Grey matter mask) 
%for k=1:3 %to run through subjects 1 to 3
numlist = {1,2,3,5,7,8,9,10,11,12,13,14}; %Creates a list with the subject numbers to run through called 'numlist'
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
%row={-60:4:60}'; %row=(-60:4:60)';
% Headers %
%Each parameter is entered mannually cos it does not like the loop
%Headers = {'SubjectNo_MRIParam','sub01default','sub01paramA','sub01paramB',...
% Headers = {'sub01default','sub01paramA','sub01paramB',...
% 'sub02default','sub02paramA','sub02paramB',...
% 'sub03default','sub03paramA','sub03paramB',...
% 'sub05default','sub05paramA','sub05paramB',...
% 'sub07default','sub07paramA','sub07paramB',...
% 'sub08default','sub08paramA','sub08paramB'...
% 'sub09default','sub09paramA','sub09paramB'}; 

Headers = {'sub01default','sub01paramA','sub01paramB',...
'sub02default','sub02paramA','sub02paramB',...
'sub03default','sub03paramA','sub03paramB',...
'sub05default','sub05paramA','sub05paramB',...
'sub07default','sub07paramA','sub07paramB',...
'sub08default','sub08paramA','sub08paramB'...
'sub09default','sub09paramA','sub09paramB'...
'sub10default','sub10paramA','sub10paramB'...
'sub11default','sub11paramA','sub11paramB'...
'sub12default','sub12paramA','sub12paramB'...
'sub13default','sub13paramA','sub13paramB'...
'sub14default','sub14paramA','sub14paramB'...
};

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
saveas(gcf,fullfile(FiguresOut,['zstatmeans_all_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save
figure;
plot(Tablezstatmedians);
title('Tablezstatmedians')
saveas(gcf,fullfile(FiguresOut,['zstatmedians_all_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save
figure;
plot(TablezstatmeansNOZ);
title('TablezstatmeansNOZ')
saveas(gcf,fullfile(FiguresOut,['zstatmeansNOZ_all_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save
figure;
plot(TablezstatmediansNOZ);
title('TablezstatmediansNOZ')
saveas(gcf,fullfile(FiguresOut,['zstatmediansNOZ_all_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save


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
saveas(gcf,fullfile(FiguresOut,['zstatmeans_means_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save


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
saveas(gcf,fullfile(FiguresOut,['zstatmedians_means_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save

% Plots mean of the default, paramA & paramB means NOZ %
figure
plot((-60:4:60),mean(TablezstatmeansNOZ(:,1:3:end),2),'r')
hold on
plot((-60:4:60),mean(TablezstatmeansNOZ(:,2:3:end),2),'g')
plot((-60:4:60),mean(TablezstatmeansNOZ(:,3:3:end),2),'b')
title('Zstat Means of Means NOZ');
grid
legend('Default','paramA','paramB')
hold off
saveas(gcf,fullfile(FiguresOut,['zstatmeans_meansNOZ_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save

% Plots mean of the default, paramA & paramB medians NOZ %
figure
plot((-60:4:60),mean(TablezstatmediansNOZ(:,1:3:end),2),'r')
hold on %plots the following two onto the same graph
plot((-60:4:60),mean(TablezstatmediansNOZ(:,2:3:end),2),'g')
plot((-60:4:60),mean(TablezstatmediansNOZ(:,3:3:end),2),'b')
title('Zstat Means of Medians NOZ');
grid
legend('Default','paramA','paramB') %default = red, paramA = green, paramB = blue
%this is not RGB friendly so might change the colours later. 
hold off
saveas(gcf,fullfile(FiguresOut,['zstatmediansNOZ_means_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save


%% BoxPlot %%


%Medians - all subjects
tmp_Tablezstatmedians{1} = Tablezstatmedians(:,1:3:end); 
tmp_Tablezstatmedians{2} = Tablezstatmedians(:,2:3:end); 
tmp_Tablezstatmedians{3} = Tablezstatmedians(:,3:3:end); 

y = num2cell(1:numel(tmp_Tablezstatmedians)); 
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], tmp_Tablezstatmedians, y, 'UniformOutput', 0); % adding labels to the cells 
X = vertcat(x{:}); 
figure;
boxplot(X(:,1), X(:,2), 'Labels',{'Default','paramA', 'paramB'})
title('zstatmedians all subjects')
saveas(gcf,fullfile(FiguresOut,['zstatmedians_boxplot_allsubjs' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save

%MediansNOZ - all subjects
tmp_TablezstatmediansNOZ{1} = TablezstatmediansNOZ(:,1:3:end); 
tmp_TablezstatmediansNOZ{2} = TablezstatmediansNOZ(:,2:3:end); 
tmp_TablezstatmediansNOZ{3} = TablezstatmediansNOZ(:,3:3:end); 

y = num2cell(1:numel(tmp_TablezstatmediansNOZ)); 
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], tmp_TablezstatmediansNOZ, y, 'UniformOutput', 0); % adding labels to the cells 
X = vertcat(x{:}); 
figure;
boxplot(X(:,1), X(:,2), 'Labels',{'Default','paramA', 'paramB'})
title('zstatmediansNOZ all subjects')
saveas(gcf,fullfile(FiguresOut,['zstatmediansNOZ_boxplot_allsubjs' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save

%% sandpit
% Data1=tmp_Tablezstatmedians{1}
% Data2=tmp_Tablezstatmedians{2}
% Data3=tmp_Tablezstatmedians{3}
% allData = {Data1; Data2; Data3}; 
% h = boxplot(cell2mat(allData),group); % old version: h = boxplot([allData{:}],group);
% set(h, 'linewidth' ,2)
% set(gca,'XTickLabel', {'Data1'; 'Data2'; 'Data3'})
% hold on
% xCenter = 1:numel(allData); 
% spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
% for i = 1:numel(allData)
%     plot(rand(size(allData{i}))*spread -(spread/2) + xCenter(i), allData{i}, 'mo','linewidth', 2)
% end


%MediansNOZ

tmp_rows=(-60:4:60)'
single(tmp_rows)
tmp_TablezstatmediansNOZ{1} = TablezstatmediansNOZ(:,1:3:end); 
tmp_TablezstatmediansNOZ{2} = TablezstatmediansNOZ(:,2:3:end); 
tmp_TablezstatmediansNOZ{3} = TablezstatmediansNOZ(:,3:3:end);
%tmp_TablezstatmediansNOZ{4}=tmp_rows

y = num2cell(1:numel(tmp_TablezstatmediansNOZ)); 
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], tmp_TablezstatmediansNOZ, y, 'UniformOutput', 0); % adding labels to the cells 
X = vertcat(x{:}); 
figure;
boxplot(X(:,1), X(:,2), 'Labels',{'Default','paramA', 'paramB'})
title('zstatmediansNOZ')


%test1{1}=table(tmp_rows, tmp_TablezstatmediansNOZ{1}(1,:))
test1{1}=table(tmp_rows, tmp_TablezstatmediansNOZ{1}(1,:))
test2{1}=tmp_TablezstatmediansNOZ{1}(1,:)
y = num2cell(1:numel(tmp_TablezstatmediansNOZ)); 
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], tmp_TablezstatmediansNOZ, y, 'UniformOutput', 0); % adding labels to the cells 
X = vertcat(x{:}); 
figure;
boxplot(X(:,1), X(:,2))
%boxplot(X(:,1), X(:,2), 'Labels',{'Default','paramA', 'paramB'})
title('zstatmediansNOZ')

figure;
boxplot(tmp_TablezstatmediansNOZ{1}(1,:), tmp_TablezstatmediansNOZ{1}(2,:),tmp_TablezstatmediansNOZ{1}(3,:),tmp_TablezstatmediansNOZ{1}(4,:),tmp_TablezstatmediansNOZ{1}(5,:),tmp_TablezstatmediansNOZ{1}(6,:),tmp_TablezstatmediansNOZ{1}(7,:))
%%

% figure
% boxplot((-60:4:60),mean(Tablezstatmedians(:,1:3:end),2),'r')
% hold on %plots the following two onto the same graph
% boxplot((-60:4:60),mean(Tablezstatmedians(:,2:3:end),2),'g')
% boxplot((-60:4:60),mean(Tablezstatmedians(:,3:3:end),2),'b')
% title('Zstat Means of Medians');
% grid
% legend('Default','paramA','paramB') %default = red, paramA = green, paramB = blue
% %this is not RGB friendly so might change the colours later. 
% hold off
% 
% 
% allthedefaults=reshape(Tablezstatmedians(:,1:3:end),64*64*12,31);
% allthedefauts=Tablezstatmedians(:,1:3:end)
% alltheparamAs=Tablezstatmedians(:,2:3:end)
% alltheparamBs=Tablezstatmedians(:,3:3:end)
% 
% 
% group = [ ones(size(Data1));
% 
% 2 * ones(size(Data2))
% 
% 3 * ones(size(Data3))];
% 
% boxplot([Data1; Data2; Data3],group)
% 
% h = boxplot([Data1; Data2; Data3],group)
% 
% set(h,{'linew'},{2})
% 
% set(gca,'XTickLabel', {'Data1'; 'Data2'; 'Data3'})
% 
% %dlmwrite(sprintf('%s/zstatmeans_final_all_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),Table_new,'\t')
% 
% % save('savefile.mat', 'T')
% % save(sprintf('%s/zstatmeans_table_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),Table_new,'\t')
% % save('zstatmeans_table.tsv', 'Table_new')
% % dlmwrite(sprintf('%s/zstatmeans_table_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),Table_new,'\t') 
% % 
% % dlmwrite('zstatmeans_table',Table_new,'\t') 
% 
% %% BOXPLOTS %%
% % Tablezstatmedians(:,3)=Tablezstatmedians(:,1).*3; 
% % 
% % boxplot(Tablezstatmedians,'Notch','off','Labels',{'mu = 5','mu = 6','mu = 6'},'Whisker',1) 
% % 
% % lines = findobj(gcf, 'type', 'line', 'Tag', 'Median'); 
% % 
% % set(lines, 'Color', 'g'); 
% % 
% % % Change the boxplot color from blue to green 
% % 
% % a = get(get(gca,'children'),'children'); % Get the handles of all the objects 
% % 
% % %t = get(a,'tag'); % List the names of all the objects  
% % 
% % %box1 = a(7); % The 7th object is the first box 
% % 
% % set(a, 'Color', 'r'); % Set the color of the first box to green 
% % 
% % hold on 
% % 
% % x=ones(length(Tablezstatmedians)).*(1+(rand(length(Tablezstatmedians))-0.5)/5); 
% % 
% % x1=ones(length(Tablezstatmedians)).*(1+(rand(length(Tablezstatmedians))-0.5)/10); 
% % 
% % x2=ones(length(Tablezstatmedians)).*(1+(rand(length(Tablezstatmedians))-0.5)/15); 
% 
% Default = (Tablezstatmeans(:,1:3:end),2),'r')
% hold on
% plot((-60:4:60),mean(Tablezstatmeans(:,2:3:end),2),'g')
% plot((-60:4:60),mean(Tablezstatmeans(:,3:3:end),2),'b')
% 

% 
% f1=scatter(x(:,1),Tablezstatmedians(:,1),'k','filled');f1.MarkerFaceAlpha = 0.4;hold on  
% 
% f2=scatter(x1(:,2).*2,Tablezstatmedians(:,2),'k','filled');f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on 
% 
% f3=scatter(x2(:,3).*3,Tablezstatmedians(:,3),'k','filled');f3.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on 

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