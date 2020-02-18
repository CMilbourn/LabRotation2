% CO2 Processing %
%process a single CVR CO2 mean

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
%% Find and load files %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

src= '/Users/colette/sourcedata/' %source data path 
%subj= '2' %subject number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subject List
numlist = {2,3,5,7,8,9,10,11,12,13,14}; %Creates a list with the subject numbers to run through called 'numlist'
%numlist = {2,3,5,7,8,9,10,13,14};
for a = 1:length(numlist) %start of for loop - runs through one number in 'numlist' per iteration
%% Variables
%subj=numlist{a}
subj= sprintf('%02d',numlist{a})
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
    for k=1:3 %loop through 3 of the MRI parameters row 1 = default, row2 = paramA row 3 = paramB 

        start_ind(k,a)=find(bbb_array(:,1)>events_start(k),1); %find start index number from bbb_array using the event from eventsfile 'row k' 
        end_ind(k,a)=find(bbb_array(:,1)>(events_start(k)+7),1); %find the end index number by adding 7 to the start_inx - assume the full time is 7mins, add 7mins to the  
       %figure();hold on 
%       plot(bbb_array(start_ind(k,:):end_ind(k,:),1)-bbb_array(start_ind(k,:),1),bbb_array(start_ind(k,:):end_ind(k,:),2)) %plot the three windows on a graph overlaying each other  
% 
%       title('BBB array Windows'); 
%       ylabel('PCO2'); 
%       xlabel('Time (Minutes)'); 
                 %before hypercapnia
                start_bef_ind= start_ind(k,a);%time point  %full window - 0 to 2mins
                end_bef_ind= find(bbb_array(:,1)>(events_start(k)+1.5),1); 
                fprintf('start_bef_ind is %d ',start_bef_ind); 
                fprintf('end_bef_ind is %d ',end_bef_ind);            
    %           start_bef_co2=bbb_array(start_bef_ind,2); 
    %           fprintf(start_bef_co2) 
    %           end_bef_co2= 
    %           bbb_array(end_bef_ind,2);(fprintf(end_bef_co2) 
                window_bef=bbb_array(start_bef_ind:end_bef_ind,2) %select range of CO2 column from start to end of bef indices  
                bef_mean(k,a)=mean(window_bef) % find mean of bef window 
%                 
%                 start_ind.bef= start_ind(k,a);%time point  %full window - 0 to 2mins
%                 end_ind.bef= find(bbb_array(:,1)>(events_start(k)+1.5),1); 
%                 fprintf('start_bef_ind is %d ',start_ind.bef); 
%                 fprintf('end_bef_ind is %d ',end_ind.bef);            
%     %           start_bef_co2=bbb_array(start_bef_ind,2); 
%     %           fprintf(start_bef_co2) 
%     %           end_bef_co2= 
%     %           bbb_array(end_bef_ind,2);(fprintf(end_bef_co2) 
%                 window.bef=bbb_array(start_ind.bef:end_ind.bef,2) %select range of CO2 column from start to end of bef indices  
%                mean.bef(k,a)=mean(window.bef) % find mean of bef window 
%                 
%                 %during hypercapnia window
                start_dur_ind=find(bbb_array(:,1)>(events_start(k)+2),1); %full window = 2 to 5mins 
                end_dur_ind=find(bbb_array(:,1)>(events_start(k)+5),1); 
                fprintf('start_dur_ind is %d ', start_dur_ind); 
                fprintf('end_dur_ind is %d ', end_dur_ind); 

                window_dur=bbb_array(start_dur_ind:end_dur_ind,2); %select range of CO2 column from start to end of default indices  
                dur_mean(k,a)=mean(window_dur) % find mean of default window         
                %after hypercapnia window 
                start_aft_ind=find(bbb_array(:,1)>(events_start(k)+5),1); %full window = 5 to 7 mins 
                end_aft_ind= find(bbb_array(:,1)>(events_start(k)+7),1); 
                fprintf('start_aft_ind is %d ', start_aft_ind);
                fprintf('end_aft_ind is %d ', end_aft_ind); 
 
                window_aft=bbb_array(start_aft_ind:end_aft_ind,2); %select range of CO2 column from start to end of default indices  
                aft_mean(k,a)=mean(window_aft) % find mean of default window 
                
                %Average CO2 Change == Avereage Baseline - average hypercapnia 
                baseline(k,a)= (bef_mean(k,a) + aft_mean(k,a))/2 %before + after baseline divided by 2 == baseline
                avgco2change(k,a)=baseline(k,a) - dur_mean(k,a) 
                
                %rownames(:,a)=
%                 
%                     if ismember(k,[1])
%                         rownames(k,:)='default'
%                     elseif ismember(k,[2])
%                         rownames(k,:)='paramA'
%                     else ismember(k,[3])
%                         rownames(k,:)='paramB'
%                     end
                 %Headers(:,a)=(sprintf('SubjNo_%s',a))
                %Table(k,a . ) = avgco2change(k,a), 'RowNames',rownames(:,a), 'VariableNames',Headers)
                table_co2_changes(k,a) = avgco2change(k,a)
                
                
    end 
end 

%extract from loop or maybe just save them?
base=baseline
avgco2change=avgco2change
T=table_co2_changes

% Calculate the mean and SD 
