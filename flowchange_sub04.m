%%% Flow Change %%%
%part 1 is cvr_processing7.m from 20200218 - CVRmediansNOZ {1,2,3,5,7,8,9,10,11,12,13,14}
%part 2 is CO2_processing2_3.m from 20200218 numlist = {1,2,3,5,7,8,9,10,13,14};
%part 3 combine and extract the flow change

%part 1:
%What does this script do?
%set up FSL environment
%Read in cvrmap percentage data
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
% clear all % clear workspace
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
%%subjects to run %%
%subjlist = {1,2,3,5,6,7,8,9,10,11,12,13,14,4}; %Creates a list with the subject numbers to run through called 'numlist'
subjlist = {4}; 
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
%% PART 1 - For Loop to create files of median cvrs, without 0's %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(i.e. without 0's excludes areas outside the Grey matter mask) 
%for k=1:3 %to run through subjects 1 to 3

for Subj = 1:length(subjlist) %start of for loop - runs through one number in 'numlist' per iteration
    %% Variables
    subj= sprintf('%02d',subjlist{Subj}) %assigns the variable 'subj' the number from numlist as a string in 01 format
    src = '/Users/colette/sourcedata/' %source in
    srcout2 = sprintf('%sderivatives/sub-%s/sub-%s_', src, subj, subj) %alternative source out if it is not the same as srcout

    srcout = '/Users/colette/sourcedata/derivatives/StatsOutput/'
    src = '/Users/colette/sourcedata/'
    %srcout = sprintf('/Users/colette/sourcedata/derivatives/sub-%s', subj)
    srcout2 = sprintf('%sderivatives/sub-%s/sub-%s_', src, subj, subj)
    
    %% *read_avw*
    
    %default_cvr
    MRIParam1= 'default'
    [default_cvr, dims, scales, bpp, endian] = read_avw([srcout2 MRIParam1 '_analysis/CVRmaps_concact_sub-' subj '_' MRIParam1 '_MZeroScan.nii.gz']);
    
    MRIParam2 = 'paramA'
    [paramA_cvr, dims, scales, bpp, endian] = read_avw([srcout2 MRIParam2 '_analysis/CVRmaps_concact_sub-' subj '_' MRIParam2 '_MZeroScan.nii.gz']);
%     
%     MRIParam3 = 'paramB'
%     [paramB_cvr, dims, scales, bpp, endian] = read_avw([srcout2 MRIParam3 '_analysis/CVRmaps_concact_sub-' subj '_' MRIParam3 '_MZeroScan.nii.gz']);
%     
    %% gm_cortex %
    %Reads in the grey matter cortex mask created by 'do_analysis_edit2.sh'
    %e.g. file: /Users/colette/sourcedata/derivatives/sub-01/sub-01_gmcortex.nii.gz
    gm_cortex = read_avw([sprintf('%sderivatives/sub-%s/sub-%s_gmcortex.nii.gz', src, subj, subj)]);
    % *Output*: variable gm_cortex
    
    %% figure for checking image input
%     figure; %creates a figure for one slice/volume 
%     imagesc(default_cvr(:,:,6),[0 2]); %displays default_cvr as a sanity check
%     colorbar; %colour bar displayed in figure
    
    %% reshape
    %Reshapes the 4D singles into 419152x31 single - appends '_r' to new
    %variable
    default_cvr_r=reshape(default_cvr,64*64*12,31); %64x64x12 matrix, 31 timepoints(skip numbers)
    paramA_cvr_r=reshape(paramA_cvr,64*64*12,31);
%     paramB_cvr_r=reshape(paramB_cvr,64*64*12,31);
    gm_cortex_r=reshape(gm_cortex,64*64*12,1);
    
    %% Make Mask for gm_cortex
    default_cvr_vals=default_cvr_r(gm_cortex_r>0,:);
    paramA_cvr_vals=paramA_cvr_r(gm_cortex_r>0,:);
%     paramB_cvr_vals=paramB_cvr_r(gm_cortex_r>0,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    %% Create Mean and Median variables for table %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Make Medians without 0's %%
    default_cvr_valsNOZ=default_cvr_vals;
    default_cvr_valsNOZ(default_cvr_valsNOZ(:,1)==0,:)=NaN;
    cvrmedians_NOZ(1,:)= nanmedian(default_cvr_valsNOZ)
    
    paramA_cvr_valsNOZ=paramA_cvr_vals;
    paramA_cvr_valsNOZ(paramA_cvr_valsNOZ(:,1)==0,:)=NaN;
    cvrmedians_NOZ(2,:)=nanmedian(paramA_cvr_valsNOZ)
    
%     paramB_cvr_valsNOZ=paramB_cvr_vals;
%     paramB_cvr_valsNOZ(paramB_cvr_valsNOZ(:,1)==0,:)=NaN;
%     cvrmedians_NOZ(3,:)=nanmedian(paramB_cvr_valsNOZ)
%     
    %save cvrmedians to one variable column per iteration
    cvrmediansNOZ_final{Subj}=(cvrmedians_NOZ)'
    dlmwrite(sprintf('%s/cvrmediansNOZ_final_sub-%s_%s.tsv', DataOutputPerSubj, subj, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),cvrmediansNOZ_final{Subj},'\t')
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Save total cvr files as .tsv files %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%cvrmediansNOZ
%dlmwrite(sprintf('%s/cvrmediansNOZ_final_all_%s.tsv', DataOutput,datestr(clock,'yyyy-mm-dd_HH-MM-SS')),cvrmediansNOZ_final,'\t')
dlmwrite(sprintf('%s/cvrmediansNOZ_final_sub04_%s.tsv', DataOutput,datestr(clock,'yyyy-mm-dd_HH-MM-SS')),cvrmediansNOZ_final,'\t')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up table inputs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%rows - skip numbers -60 secs to +60 secs in intervals of 4secs
%SkipNumbers = {(-60:4:60)'}; %row=(-60:4:60)';
rowNames = {'-60','-56','-52','-48','-44','-40','-36','-32','-28','-24','-20','-16','-12','-8','-4','0','4','8','12','16','20',...
'24','28','32','36','40','44','48','52','56','60'}
rowNames=rowNames'
% 
% Headers = {'sub01default','sub01paramA','sub01paramB',...
% 'sub02default','sub02paramA','sub02paramB',...
% 'sub03default','sub03paramA','sub03paramB',...
% 'sub05default','sub05paramA','sub05paramB',...
% 'sub06default','sub06paramA','sub06paramB',...
% 'sub07default','sub07paramA','sub07paramB',...
% 'sub08default','sub08paramA','sub08paramB'...
% 'sub09default','sub09paramA','sub09paramB'...
% 'sub10default','sub10paramA','sub10paramB'...
% 'sub11default','sub11paramA','sub11paramB'...
% 'sub12default','sub12paramA','sub12paramB'...
% 'sub13default','sub13paramA','sub13paramB'...
% 'sub14default','sub14paramA','sub14paramB'...
% 'sub4default','sub4paramA',...
% };

Headers = {
'sub4default','sub4paramA',...
};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TABLES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timestamp=datestr(clock,'yyyy-mm-dd_HH-MM-SS') 

%cvrmediansNOZ
TablecvrmediansNOZ = cell2mat(cvrmediansNOZ_final)
cvrmediansNOZ_final_table= array2table(TablecvrmediansNOZ,'RowNames',rowNames,'VariableNames',Headers) 
writetable(cvrmediansNOZ_final_table, sprintf('%s/cvrmediansNOZ_final_table_%s.csv',FiguresOut,(datestr(clock,'yyyy-mm-dd_HH-MM-SS'))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT medians with & without 0's %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots Tablecvr's for mean, median both w. & w.o 0's
%adds a title
figure;
plot(TablecvrmediansNOZ);
title('TablecvrmediansNOZ')
%saveas(gcf,fullfile(FiguresOut,['cvrmediansNOZ_all_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save
saveas(gcf,fullfile(FiguresOut,['cvrmediansNOZ_sub04_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save

% Plots mean of the default, paramA & paramB medians NOZ %
figure
plot((-60:4:60),mean(TablecvrmediansNOZ(:,1:3:end),2),'r')
hold on %plots the following two onto the same graph
plot((-60:4:60),mean(TablecvrmediansNOZ(:,2:3:end),2),'g')
%plot((-60:4:60),mean(TablecvrmediansNOZ(:,3:3:end),2),'b')
title('cvr Means of Medians NOZ');
grid
%legend('Default','paramA','paramB') %default = red, paramA = green, paramB = blue
legend('Default','paramA') %default = red, paramA = green, paramB = blue
%this is not RGB friendly so might change the colours later. 
hold off
% saveas(gcf,fullfile(FiguresOut,['cvrmediansNOZ_means_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save
saveas(gcf,fullfile(FiguresOut,['cvrmediansNOZ_means_sub04' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save

%% BoxPlots %%

%% All Subjects, all skip numbers by MRI Param
%MediansNOZ - all subjects
tmp_TablecvrmediansNOZ{1} = TablecvrmediansNOZ(:,1:3:end); 
tmp_TablecvrmediansNOZ{2} = TablecvrmediansNOZ(:,2:3:end); 
%tmp_TablecvrmediansNOZ{3} = TablecvrmediansNOZ(:,3:3:end); 

y = num2cell(1:numel(tmp_TablecvrmediansNOZ)); 
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], tmp_TablecvrmediansNOZ, y, 'UniformOutput', 0); % adding labels to the cells 
X = vertcat(x{:}); 
figure;
%boxplot(X(:,1), X(:,2), 'Labels',{'Default','paramA', 'paramB'})
boxplot(X(:,1), X(:,2), 'Labels',{'Default','paramA'})
title('cvrmediansNOZ subjects')
% saveas(gcf,fullfile(FiguresOut,['cvrmediansNOZ_boxplot_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save
saveas(gcf,fullfile(FiguresOut,['cvrmediansNOZ_boxplot_sub04' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save
%%All subjects, only skip number 19 (+12 seconds), by MRI param
%% All Subjects, by MRI Param
%MediansNOZ - all subjects for time point +12seconds from center point
%(skip number col 19)
tmp_TablecvrmediansNOZ_12{1} = TablecvrmediansNOZ(19,1:3:end); 
tmp_TablecvrmediansNOZ_12{2} = TablecvrmediansNOZ(19,2:3:end); 
%tmp_TablecvrmediansNOZ_12{3} = TablecvrmediansNOZ(19,3:3:end); 

y = num2cell(1:numel(tmp_TablecvrmediansNOZ_12)); 
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], tmp_TablecvrmediansNOZ_12, y, 'UniformOutput', 0); % adding labels to the cells 
X = vertcat(x{:}); 
figure;
%boxplot(X(:,1), X(:,2), 'Labels',{'Default','paramA', 'paramB'})
boxplot(X(:,1), X(:,2), 'Labels',{'Default','paramA'})
title('cvrmediansNOZ all subjects, skip number +12 seconds')
%saveas(gcf,fullfile(FiguresOut,['cvrmediansNOZ_boxplot__12secs' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save
saveas(gcf,fullfile(FiguresOut,['cvrmediansNOZ_boxplot_syub04_12secs' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save

fprintf('~ ~ ~ End of part 1 cvr median NOZ processing ~ ~ ~ ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2 - CO2 change calculation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Subj= 1:length(subjlist) %start of for loop - runs through one number in 'numlist' per iteration
%% Variables
%subj=numlist{a} %Subj=in number form subj = in string form
subj= sprintf('%02d',subjlist{Subj})
sub_folder= sprintf('%s/sub-%s/CVR/', src, subj) %subject CVR sub_folder

%Load breath by breath file%
bbb_mod= sprintf('%s/sub-%s_bbb_mod.txt', sub_folder,subj) %subj bbb_mod.txt file - breath by breath file path
bbb_table=readtable(bbb_mod,'Delimiter','	'); %imports as a table
%check bbb_mod file exists
if ~exist(bbb_mod, 'file');
    fprintf('bbb_mod file not found: %s\n', bbb_mod);
else
    fprintf('bbb_mod file found!\n');
end

%bbb file from table to matrix
bbb_array=table2array(bbb_table)
%col 1 = > PCO2Time_min
%col2 = PCO2mmHg
%col3= PO2mmHg

%%Load eventsRobustCVR file%
eventsRobustCVR= sprintf('%s/sub-%s_events_RobustCVR.txt', sub_folder,subj) %subj eventsRobustCVR.txt file - timings of default, paramA, parmaB events
events_table=readtable(eventsRobustCVR,'Delimiter','	'); %imports as a table
%check eventsRobust CVR file exists
if ~exist(eventsRobustCVR, 'file')
     fprintf('eventsRobustCVR file not found: %s\n', eventsRobustCVR);
 else
     fprintf('eventsRobustCVR file found!\n');
end

%% Load Events file
eventsfile= sprintf('%s/sub-%s_events_start.csv', sub_folder,subj);
events_start=load(eventsfile); %load eventsfile
% 
% figure; %open figure
% hold on; %hold figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through paramaters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for mriparam=1:3 %loop through 3 of the MRI parameters row 1 = default, row2 = paramA row 3 = paramB 
     for mriparam=1:2
        start_ind(mriparam,Subj)=find(bbb_array(:,1)>events_start(mriparam),1); %find start index number from bbb_array using the event from eventsfile 'row mriparam' 
        end_ind(mriparam,Subj)=find(bbb_array(:,1)>(events_start(mriparam)+7),1); %find the end index number by adding 7 to the start_inx - assume the full time is 7mins, add 7mins to the  
       %figure();hold on 
%       plot(bbb_array(start_ind(mriparam,:):end_ind(mriparam,:),1)-bbb_array(start_ind(mriparam,:),1),bbb_array(start_ind(mriparam,:):end_ind(mriparam,:),2)) %plot the three windows on a graph overlaying each other  
% 
%       title('BBB array Windows'); 
%       ylabel('PCO2'); 
%       xlabel('Time (Minutes)'); 
                 %before hypercapnia
                start_bef_ind= start_ind(mriparam,Subj);%time point  %full window - 0 to 2mins
                end_bef_ind= find(bbb_array(:,1)>(events_start(mriparam)+1.5),1); 
                fprintf('start_bef_ind is %d ',start_bef_ind); 
                fprintf('end_bef_ind is %d ',end_bef_ind);            
    %           start_bef_co2=bbb_array(start_bef_ind,2); 
    %           fprintf(start_bef_co2) 
    %           end_bef_co2= 
    %           bbb_array(end_bef_ind,2);(fprintf(end_bef_co2) 
                window_bef=bbb_array(start_bef_ind:end_bef_ind,2) %select range of CO2 column from start to end of bef indices  
                bef_mean(mriparam,Subj)=mean(window_bef) % find mean of bef window 
%                 
%                 start_ind.bef= start_ind(mriparam,a);%time point  %full window - 0 to 2mins
%                 end_ind.bef= find(bbb_array(:,1)>(events_start(mriparam)+1.5),1); 
%                 fprintf('start_bef_ind is %d ',start_ind.bef); 
%                 fprintf('end_bef_ind is %d ',end_ind.bef);            
%     %           start_bef_co2=bbb_array(start_bef_ind,2); 
%     %           fprintf(start_bef_co2) 
%     %           end_bef_co2= 
%     %           bbb_array(end_bef_ind,2);(fprintf(end_bef_co2) 
%                 window.bef=bbb_array(start_ind.bef:end_ind.bef,2) %select range of CO2 column from start to end of bef indices  
%                mean.bef(mriparam,a)=mean(window.bef) % find mean of bef window 
%                 
%                 %during hypercapnia window
                start_dur_ind=find(bbb_array(:,1)>(events_start(mriparam)+2.5),1); %full window = 2 to 5mins 
                end_dur_ind=find(bbb_array(:,1)>(events_start(mriparam)+4.5),1); 
                fprintf('start_dur_ind is %d ', start_dur_ind); 
                fprintf('end_dur_ind is %d ', end_dur_ind); 

                window_dur=bbb_array(start_dur_ind:end_dur_ind,2); %select range of CO2 column from start to end of default indices  
                dur_mean(mriparam,Subj)=mean(window_dur) % find mean of default window         
                %after hypercapnia window 
                start_aft_ind=find(bbb_array(:,1)>(events_start(mriparam)+5.5),1); %full window = 5 to 7 mins 
                end_aft_ind= find(bbb_array(:,1)>(events_start(mriparam)+7),1); 
                fprintf('start_aft_ind is %d ', start_aft_ind);
                fprintf('end_aft_ind is %d ', end_aft_ind); 
 
                window_aft=bbb_array(start_aft_ind:end_aft_ind,2); %select range of CO2 column from start to end of default indices  
                aft_mean(mriparam,Subj)=mean(window_aft) % find mean of default window 
                
                %Average CO2 Change == average hypercapnia CO2 - Avereage
                %Baseline CO2
                baseline(mriparam,Subj)= (bef_mean(mriparam,Subj) + aft_mean(mriparam,Subj))/2 %before + after baseline divided by 2 == baseline
                avgco2change(mriparam,Subj)= dur_mean(mriparam,Subj) - baseline(mriparam,Subj)
                
                %rownames(:,a)=
%                 
%                     if ismember(mriparam,[1])
%                         rownames(mriparam,:)='default'
%                     elseif ismember(mriparam,[2])
%                         rownames(mriparam,:)='paramA'
%                     else ismember(mriparam,[3])
%                         rownames(mriparam,:)='paramB'
%                     end
                % Headers(:,a)=(sprintf('SubjNo_%s',a))
                %Table(mriparam,a . ) = avgco2change(mriparam,a), 'RowNames',rownames(:,a), 'VariableNames',Headers)
                table_co2_changes(mriparam,Subj) = avgco2change(mriparam,Subj)
                
                %flow change per mmHg
                %could pick one with a slight delay shift ~ could pick one
                %w. optimal corrolation 
                
                
                
                
    end 
end 

%% extract from loop CO2 changes and put into a table
base=baseline
avgco2=avgco2change

%mriparameters = {'default', 'paramA', 'paramB'}
mriparameters = {'default', 'paramA'}
mriparameters = (mriparameters)'
% Headers2={'sub01',...
% 'sub02',...
% 'sub03',...
% 'sub05',...
% 'sub06',...
% 'sub07',...
% 'sub08'...
% 'sub09'...
% 'sub10'...
% 'sub11'...
% 'sub12'...
% 'sub13'...
% 'sub14'...
% 'sub4'...
% };
Headers2={
'sub4'...
};
T=table_co2_changes
co2_changes_final= array2table(table_co2_changes,'RowNames',mriparameters,'VariableNames',Headers2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 3 - calculate flow change %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate flow change

% for refence
% tmp_TablecvrmediansNOZ_12{1} = TablecvrmediansNOZ(19,1:3:end); %default,
% skip no. 119
% tmp_TablecvrmediansNOZ_12{2} = TablecvrmediansNOZ(19,2:3:end); %paramA
% tmp_TablecvrmediansNOZ_12{3} = TablecvrmediansNOZ(19,3:3:end);  %paramB
%dlmwrite(sprintf('%s/tmp_TablecvrmediansNOZ_12%s.tsv', DataOutput,datestr(clock,'yyyy-mm-dd_HH-MM-SS')),tmp_TablecvrmediansNOZ_12,'\t')


%re-format tmp_TablecvrmediansNOZ_12

%cells_cvrmediansNOZ_12={tmp_TablecvrmediansNOZ_12{1};tmp_TablecvrmediansNOZ_12{2};tmp_TablecvrmediansNOZ_12{3}} %format transverse, col = subjno, row = mriparam
cells_cvrmediansNOZ_12={tmp_TablecvrmediansNOZ_12{1};tmp_TablecvrmediansNOZ_12{2}} %format transverse, col = subjno, row = mriparam

array_cvrmediansNOZ_12=cell2mat(cells_cvrmediansNOZ_12) %format to array col = subjno, row = mriparam

Table_cvrmediansNOZ_12= array2table(array_cvrmediansNOZ_12,'RowNames',mriparameters,'VariableNames',Headers2)


%flow change = cvr divided by average co2
array_flowchange = array_cvrmediansNOZ_12./avgco2
dlmwrite(sprintf('%s/array_flowchange_%s.tsv', DataOutput,datestr(clock,'yyyy-mm-dd_HH-MM-SS')),array_flowchange,'\t')
Table_flowchange= array2table(array_flowchange,'RowNames',mriparameters,'VariableNames',Headers2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boxplots %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flow_change_transpose=(array_flowchange)';
figure;
Data1=flow_change_transpose(:,1)
Data2=flow_change_transpose(:,2)
%Data3=flow_change_transpose(:,3)
%allData = {Data1; Data2; Data3}; 
allData = {Data1; Data2}; 


group = [    ones(size(Data1));
         2 * ones(size(Data2))
         %3 * ones(size(Data3))];
%boxplot([Data1; Data2; Data3],group)
boxplot([Data1; Data2],group)

h = boxplot(cell2mat(allData),group); % old version: h = boxplot([allData{:}],group);
set(h, 'linewidth' ,2)
%set(gca,'XTickLabel', {'Default';'paramA'; 'paramB'})
set(gca,'XTickLabel', {'Default';'paramA'})
hold on
xCenter = 1:numel(allData); 
spread = 0; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(allData)
    plot(rand(size(allData{i}))*spread -(spread/2) + xCenter(i), allData{i}, 'mo','linewidth', 2)
end
ylim([0 10])
get(h,'tag')
set(h(6,:),'color','k','linewidth',3.0)
colorchoice=['y','g','b'];
h=findobj(gca,'Tag','Box')'
for i = 1:length(h)
        patch(get(h(i),'XData'),get(h(i),'YData'),colorchoice(i),'FaceAlpha',0.5);
end


%% BOXPLOT - Young participants only

flow_change_transpose=(array_flowchange)';
figure;
Data1=flow_change_transpose(1:9,1)
Data2=flow_change_transpose(1:9,2)
%Data3=flow_change_transpose(1:9,3)
allDatayoung = {Data1; Data2; Data3}; 


group = [    ones(size(Data1));
         2 * ones(size(Data2))
         3 * ones(size(Data3))];
boxplot([Data1; Data2; Data3],group)

h = boxplot(cell2mat(allDatayoung),group); % old version: h = boxplot([allDatayoung{:}],group);
set(h, 'linewidth' ,2)
set(gca,'XTickLabel', {'Default';'paramA'; 'paramB'})
hold on
xCenter = 1:numel(allDatayoung); 
spread = 0; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(allDatayoung)
    plot(rand(size(allDatayoung{i}))*spread -(spread/2) + xCenter(i), allDatayoung{i}, 'mo','linewidth', 2)
end
title('Young participants only')
get(h,'tag')
ylim([0 10])
set(h(6,:),'color','k','linewidth',3.0)
colorchoice=['y','g','b'];
h=findobj(gca,'Tag','Box')'
for i = 1:length(h)
        patch(get(h(i),'XData'),get(h(i),'YData'),colorchoice(i),'FaceAlpha',0.5);
end



%% BOXPLOT - old participants only

flow_change_transpose=(array_flowchange)';
figure;
Data4=flow_change_transpose(10:end,1)
Data5=flow_change_transpose(10:end,2)
%Data6=flow_change_transpose(10:end,3)
allDataold = {Data4; Data5; Data6}; 


group = [    ones(size(Data4));
         2 * ones(size(Data5))
%         3 * ones(size(Data6))];
%boxplot([Data4; Data5; Data6],group)
boxplot([Data4; Data5],group)

h = boxplot(cell2mat(allDataold),group); % old version: h = boxplot([allDataold{:}],group);
set(h, 'linewidth' ,2)
%set(gca,'XTickLabel', {'Default';'paramA'; 'paramB'})
set(gca,'XTickLabel', {'Default';'paramA'})
hold on
xCenter = 1:numel(allDataold); 
spread = 0; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(allDataold)
    plot(rand(size(allDataold{i}))*spread -(spread/2) + xCenter(i), allDataold{i}, 'mo','linewidth', 2)
end
title('Old participants only')
get(h,'tag')
ylim([0 10])
set(h(6,:),'color','k','linewidth',3.0)
colorchoice=['y','g','b'];
h=findobj(gca,'Tag','Box')'
for i = 1:length(h)
        patch(get(h(i),'XData'),get(h(i),'YData'),colorchoice(i),'FaceAlpha',0.5);
end


%% Boxplot young vs. old
figure;

%allDatanew = {Data1;Data2; Data3; Data4; Data5; Data6}; 
allDatanew = {Data1;Data2; Data4; Data5;}; 

% group = [   ones(size(Data1));
%          2 * ones(size(Data2))
%          3 * ones(size(Data3))
%          4 * ones(size(Data4));
%          5 * ones(size(Data5))
%          6 * ones(size(Data6))];

group = [   ones(size(Data1))
         2 * ones(size(Data4))
         3 * ones(size(Data2))
         4 * ones(size(Data5))
%         5 * ones(size(Data3))
%         6 * ones(size(Data6))];
%boxplot([Data1; Data4; Data2; Data5; Data3; Data6],group)
boxplot([Data1; Data4; Data2; Data5;],group)

h = boxplot(cell2mat(allDatanew),group); % old version: h = boxplot([allDatanew{:}],group);
set(h, 'linewidth' ,2)
%set(gca,'XTickLabel', {'Defaultyoung';'Defaultold';'paramAyoung'; 'paramAold'; 'paramByoung'; 'paramBold'})
set(gca,'XTickLabel', {'Defaultyoung';'Defaultold';'paramAyoung'; 'paramAold'})
hold on
xCenter = 1:numel(allDatanew); 
spread = 0; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(allDatanew)
    plot(rand(size(allDatanew{i}))*spread -(spread/2) + xCenter(i), allDatanew{i}, 'mo','linewidth', 2)
end
title('Compare participants young v.s old')
get(h,'tag')
ylim([0 10])
set(h(6,:),'color','k','linewidth',3.0)
colorchoice=['y','g','b'];
h=findobj(gca,'Tag','Box')'
for i = 1:length(h)
        patch(get(h(i),'XData'),get(h(i),'YData'),colorchoice(i),'FaceAlpha',0.5);
end

%% WIP  from here %%%
% flow_change_transpose=(array_flowchange)'
% figure;
% H=boxplot(flow_change_transpose,'Labels',{'Default','paramA', 'paramB'})
% title('Flow Change')
% xlabel('MRI Parameter')
% ylabel(' Flow Change? (%)')
% hold on
% %plot data points
% xCenter = 1:numel(flow_change_transpose); 
% spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
% 
% for i = 1:numel(flow_change_transpose)
%     plot(rand(size(flow_change_transpose(:,1)))*spread -(spread/2) + xCenter(i), flow_change_transpose, 'yo','linewidth', 1)
% end
% 
% for i = 1:numel(flow_change_transpose)
%     plot(rand(size(flow_change_transpose(:,2)))*spread -(spread/2) + xCenter(i), flow_change_transpose, 'go','linewidth', 1)
% end
% 
% for i = 1:numel(flow_change_transpose)
%     plot(rand(size(flow_change_transpose(:,3)))*spread -(spread/2) + xCenter(i), flow_change_transpose, 'bo','linewidth', 1)
% end
% %set line width
% set(findobj(gca,'type','line'),'linew',3)
% set(gca,'linew',2)
% %set colours of boxes;
% get(H,'tag')
% set(H(6,:),'color','k','linewidth',3.0)
% colorchoice=['y','g','b'];
% h=findobj(gca,'Tag','Box')'
% for i = 1:length(h)
%         patch(get(h(i),'XData'),get(h(i),'YData'),colorchoice(i),'FaceAlpha',0.5);
% end
% %legend('Default','paramA','paramB')
% 
% %% 
% 
% 
% plot(flow_change_transpose)
% f1=scatter(x(:,1),flow_change_transpose(:,1),'k','filled');f1.MarkerFaceAlpha = 0.4;hold on  
% 
% f2=scatter(x1(:,2).*2,flow_change_transpose(:,2),'k','filled');f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on 
% 
% f3=scatter(x2(:,3).*3,flow_change_transpose(:,3),'k','filled');f3.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on 
% 
% %% Boxplot?
% %inputdata
% 
% mriparams = {'default', 'paramA','paramB'};
% mriparams= (mriparams)'
% 
% t = table(mriparams,flow_change_transpose(:,1),flow_change_transpose(:,2),flow_change_transpose(:,3),...
% 'VariableNames',{'meas1','meas2','meas3',});
% Meas = dataset([1; 2; 3]','VarNames',{'Measurements'});
% %Fit a repeated measures model, where the measurements are the responses and the species is the predictor variable.
% 
% rm = fitrm(t,'meas1-meas3~mriparams','WithinDesign',Meas);
% %Plot data grouped by the factor species.
% 
% plot(rm,'group','mriparams')
% 
% %%
% 
% 
% 
% mriparams = {'default', 'paramA','paramB'};
% mriparams= (mriparams)'
% 
% t = table(mriparams,array_flowchange(1,:),array_flowchange(2,:),array_flowchange(3,:),...
% 'VariableNames',{'meas1','meas2','meas3',});
% Meas = dataset([1 2 3]','VarNames',{'Measurements'});
% %Fit a repeated measures model, where the measurements are the responses and the species is the predictor variable.
% 
% rm = fitrm(t,'meas1-meas3~mriparams','WithinDesign',Meas);
% %Plot data grouped by the factor species.
% 
% plot(rm,'group','mriparams')
% 
% %%
% 
% t = table(Headers2,array_flowchange(1,:),array_flowchange(2,:),array_flowchange(3,:),...
% 'VariableNames',{'Headers2','meas1','meas2','meas3'});
% Meas = dataset([1 2 3]','VarNames',{'Measurements'});
% rm = fitrm(t,'meas1-meas3~Headers2','WithinDesign',array_flowchange);
% plot(rm,'group','Headers2')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistics: Two-way ANOVA %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANOVA 1 on array_flowchange transposed
% [p tbl sts]=anova1(array_flowchange') 
% 
% %Multi compare on stats from above
% multcompare(sts)
% 
% %% Young only
% 
% [p tbl sts]=anova1(array_flowchange(:,1:8)')
% multcompare(sts)
% 
% %% plot sub 06
% %(skip number col 19)
% tmp_TablecvrmediansNOZ_12{1} = TablecvrmediansNOZ(19,1:3:end); 
% tmp_TablecvrmediansNOZ_12{2} = TablecvrmediansNOZ(19,2:3:end); 
% tmp_TablecvrmediansNOZ_12{3} = TablecvrmediansNOZ(19,3:3:end); 
% 
% y = num2cell(1:numel(tmp_TablecvrmediansNOZ_12)); 
% x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], tmp_TablecvrmediansNOZ_12, y, 'UniformOutput', 0); % adding labels to the cells 
% X = vertcat(x{:}); 
% figure;
% boxplot(X(:,1), X(:,2), 'Labels',{'Default','paramA', 'paramB'})
% title('cvrmediansNOZ all subjects, skip number +12 seconds')
% saveas(gcf,fullfile(FiguresOut,['cvrmediansNOZ_boxplot__12secs' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg'])); %save
% 
% fprintf('~ ~ ~ End of part 1 cvr median NOZ processing ~ ~ ~ ')

%% End of Script %%
fprintf('~~~ End of Script ~~~'); %prints to command window '~~~ End of Script ~~~'