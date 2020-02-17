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
subj= '2' %subject number
sub_folder= sprintf('%s/sub-0%s/CVR/', src, subj) %subject CVR sub_folder

%Load breath by breath file%
bbb_mod= sprintf('%s/sub-0%s_bbb_mod.txt', sub_folder,subj) %subj bbb_mod.txt file - breath by breath file path
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
eventsRobustCVR= sprintf('%s/sub-0%s_events_RobustCVR.txt', sub_folder,subj) %subj eventsRobustCVR.txt file - timings of default, paramA, parmaB events
events_table=readtable(eventsRobustCVR,'Delimiter','	'); %imports as a table
%check eventsRobust CVR file exists
if ~exist(eventsRobustCVR, 'file')
     fprintf('eventsRobustCVR file not found: %s\n', eventsRobustCVR);
 else
     fprintf('eventsRobustCVR file found!\n');
end

%% Load Events file
eventsfile= sprintf('%s/sub-0%s_events_start.csv', sub_folder,subj);
events_start=load(eventsfile); %load eventsfile

figure; %open figure
hold on; %hold figure

%% Loop through paramaters 
    for k=1:3 %loop through 3 of the MRI parameters row 1 = default, row2 = paramA row 3 = paramB 

        start_ind(k,:)=find(bbb_array(:,1)>events_start(k),1); %find start index number from bbb_array using the event from eventsfile 'row k' 

        end_ind(k,:)=find(bbb_array(:,1)>(events_start(k)+7),1); %find the end index number by adding 7 to the start_inx - assume the full time is 7mins, add 7mins to the  

        figure();hold on 

        plot(bbb_array(start_ind(k,:):end_ind(k,:),1)-bbb_array(start_ind(k,:),1),bbb_array(start_ind(k,:):end_ind(k,:),2)) %plot the three windows on a graph overlaying each other  

         

        title('BBB array Windows'); 

        ylabel('PCO2'); 

        xlabel('Time (Minutes)'); 

        %for a=subj 

            if ismember(k,[1]) 

                start_default_ind= start_ind(k,:);%time point  

                end_default_ind= find(bbb_array(:,1)>(events_start(k)+1.5),1); 

                fprintf('start_default_ind is %d ',start_default_ind); 

                fprintf('end_default_ind is %d ',end_default_ind);            

    %           start_default_co2=bbb_array(start_default_ind,2); 

    %           fprintf(start_default_co2) 

    %           end_default_co2= 

    %           bbb_array(end_default_ind,2);(fprintf(end_default_co2) 

                window_default=bbb_array(start_default_ind:end_default_ind,2) %select range of CO2 column from start to end of default indices  

                default_mean=mean(window_default) % find mean of default window 

  

            elseif ismember(k,[2]) 

                start_paramA_ind=start_ind(k,:); 

                end_paramA_ind=find(bbb_array(:,1)>(events_start(k)+2),1); 

                fprintf('start_paramA_ind is %d ', start_paramA_ind); 

                fprintf('end_paramA_ind is %d ', end_paramA_ind); 

  

                window_paramA=bbb_array(start_paramA_ind:end_paramA_ind,2) %select range of CO2 column from start to end of default indices  

                paramA_mean=mean(window_paramA) % find mean of default window         

  

            else ismember(k,[3]) 

                start_paramB_ind=start_ind(k,:); 

                end_paramB_ind= find(bbb_array(:,1)>(events_start(k)+1),1); 

                fprintf('start_paramB_ind is %d ', start_paramB_ind); 

                fprintf('end_paramB_ind is %d ', end_paramB_ind); 

  

                window_paramB=bbb_array(start_paramB_ind:end_paramB_ind,2) %select range of CO2 column from start to end of default indices  

                paramB_mean=mean(window_paramB) % find mean of default window 

  

  

  

            end 

        %end 

  

    end 

%end 
% 
% for k=1:3 %loop through 3 of the MRI parameters row 1 = default, row2 = paramA row 3 = paramB
%     start_ind(k,:)=find(bbb_array(:,1)>events_start(k),1); %find start index number from bbb_array using the event from eventsfile 'row k'
%     end_ind(k,:)=find(bbb_array(:,1)>(events_start(k)+7),1); %find the end index number by adding 7 to the start_inx - assume the full time is 7mins, add 7mins to the 
%     plot(bbb_array(start_ind(k,:):end_ind(k,:),1)-bbb_array(start_ind(k,:),1),bbb_array(start_ind(k,:):end_ind(k,:),2)) %plot the three windows on a graph overlaying each other 
%     
%     %for a=subj
%         if ismember(k,[1])
%             start_default_ind= find(bbb_array(:,1)>events_start(k),1);%time point 
%             end_default_ind= find(bbb_array(:,1)>(events_start(k)+1.5),1);
%             fprintf('start_default_ind is %d ',start_default_ind);
%             fprintf('end_default_ind is %d ',end_default_ind);
%             default_mean=mean(bbb_array(start_default:end_default)) % find mean of default window
%             
%           
%             indices = find(A(:,1) == 3);
%             A(indices,2)
%             
%             indices = find(A(:,1)>start_default_ind,1);
%             A= bbb_array(indices,2)
%             
%             %find bbb_array(start_default_ind:2)
%             start_default_co2=find(bbb_array(:,2)>start_default_ind,1)
%             
%             A(find(A(:,1)>2 & A(:,1)<6),2)
%             
%             start_default_co2=find(bbb_array(:,2)>start_default_ind,1)
%             
%             find(bbb_array(:,1)>events_start(k),1)
%             
%             end_default_co2=find(bbb_array(:,2)>start_end_ind,1)
%             
%         elseif ismember(k,[2])
%             start_paramA_ind=find(bbb_array(:,1)>events_start(k),1);
%             end_paramA_ind=find(bbb_array(:,1)>(events_start(k)+2),1);
%             fprintf('start_paramA_ind is %d ', start_paramA_ind);
%             fprintf('end_paramA_ind is %d ', end_paramA_ind);
%             paramA_mean=mean(bbb_array(start_paramA:end_paramA)) % find mean of paramA window
%         
%         else ismember(k,[3])
%             start_paramB_ind= find(bbb_array(:,1)>events_start(k)+0.5,1);
%             end_paramB_ind= find(bbb_array(:,1)>(events_start(k)+1),1);
%             fprintf('start_paramB_ind is %d ', start_paramB_ind);
%             fprintf('end_paramB_ind is %d ', end_paramB_ind);
%             paramB_mean=mean(bbb_array(start_paramB:end_paramB)) % find mean of paramB window
%         end
%     %end
%               
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract the times from events sequences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for k = 1:4
%   DefaultTime=EventsList(
    

%find the event_RobustCVR that has the correct times sub-01_eventsRobustCVR.txt

%c=textread(events_file,'%s','delimiter','\n');

% BBB table to cell array
%bbb_array = table2cell(bbb_table)
%[bbb_table.Properties.VariableNames;bbb_array]
%%events table to cell array
%events_array = table2cell(events_file)
%[events_file.Properties.VariableNames;events_array]

%% Events List %%



%% if clause so that if there is more than one start i.e. a restart that the one with the biggest line is used for the variable %%

% Has to take in the respdata mean and SD 
% 
% - take the different sections -> find start and finish point of each section for bef,dur,aft 

%find column
%find bef (~ 2mins)
%find dur (~3mins)
%find aft (~2mins)

%select the middle 30seconds ish to sample

% find the average of each 
% 
% Needs to be able to search both the time and the event column in unison  
% 
%  
% 
% Calculate the mean and SD 
% 
%  
% 
% See if possible to do this easily for one.  