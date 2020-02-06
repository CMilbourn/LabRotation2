%% README asl_analysis
 
%% Line to make it into a function
% function T=asl_analysis(src,subj)
%% Housekeeping
clear all % clear workspace
close all % close any open figures
clc %clear command window
diary on
diary my_data.out
echo on
%% set FSL environment 
fprintf 'Setting up FSL environment...'
setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is 
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be 

srcout = '/Users/colette/sourcedata/derivatives/StatsOutput/'
src = '/Users/colette/sourcedata/'

    if ~exist('/Users/colette/sourcedata/derivatives/DataOutput', 'dir')
       mkdir('/Users/colette/sourcedata/derivatives/DataOutput')
    end

DataOutput = '/Users/colette/sourcedata/derivatives/DataOutput'

    if ~exist('/Users/colette/sourcedata/derivatives/DataOutput/PerSubj', 'dir')
       mkdir('/Users/colette/sourcedata/derivatives/DataOutput/PerSubj')
    end

DataOutputPerSubj = '/Users/colette/sourcedata/derivatives/DataOutput/PerSubj'

%for k=1:3
%for k= [1, 2, 3, 5, 7, 8, 9]    
%for k= [1 2 3 5 7 8 9]
%for k= {'01' '02' '03' '05' '07' '08' '09'}
%for k= 1 2 3 5 7 8 9
numlist = {1,2,3,5,7,8,9};
for k = 1:length(numlist)
    %numlist = {1,2,3,5,7,8,9};
    %% Variables
    %subj = sprintf('%02d',k);
    subj= sprintf('%02d',numlist{k})
    %subj = {k};
    srcout = '/Users/colette/sourcedata/derivatives/'
    src = '/Users/colette/sourcedata/'
    srcout2 = sprintf('%sderivatives/sub-%s/sub-%s_', src, subj, subj)
    %srcout = sprintf('/Users/colette/sourcedata/derivatives/sub-%s', subj)
   % srcout2 = sprintf('%sderivatives/sub-%s/sub-%s_', src, subj, subj)
    
    %% *read_avw*
    %/Users/colette/sourcedata/derivatives/sub-01/sub-01_default_analysis/Zstat2_concat_sub-01_default_MZeroScan.nii.gz 
    %default_zstat
    MRIParam1= 'default'
    %[default_zstat, dims, scales, bpp, endian] = read_avw([srcout subj '/' 'zstat2maps_concat_sub-' subj '_' MRIParam1 '_percentage.nii.gz']);
 
    [default_zstat, dims, scales, bpp, endian] = read_avw([srcout2 MRIParam1 '_analysis/Zstat2_concat_sub-' subj '_' MRIParam1 '_MZeroScan.nii.gz']);
    
    MRIParam2 = 'paramA'
    
    [paramA_zstat, dims, scales, bpp, endian] = read_avw([srcout2 MRIParam2 '_analysis/Zstat2_concat_sub-' subj '_' MRIParam2 '_MZeroScan.nii.gz']);
    
    MRIParam3 = 'paramB'
    [paramB_zstat, dims, scales, bpp, endian] = read_avw([srcout2 MRIParam3 '_analysis/Zstat2_concat_sub-' subj '_' MRIParam3 '_MZeroScan.nii.gz']);
   
    % gm_cortex %
    %gm_cortex=read_avw(['/Users/colette/sourcedata/derivatives/sub-01/sub-01_gmcortex.nii.gz']);
    gm_cortex = read_avw([sprintf('%sderivatives/sub-%s/sub-%s_gmcortex.nii.gz', src, subj, subj)]);
    %% Reads in the cortex map that we created GM_MZeroScan and assigns it to variable
    % called gm_cortex
    % *Input*: sourcedata folder/subjectnumber/func/subject_ GM_MZeroScan
    % *Output*: variable gm_cortex
    %%figure
    
%     figure; %creates a figure
%     imagesc(default_zstat(:,:,6),[0 2]); %displays default_zstat as a sanity check
%     colorbar; %colour bar displayed in figure
    
    %% reshape
    default_zstat_r=reshape(default_zstat,64*64*12,31); %skip number 000
    %default_zstat_r=reshape(default_zstat,64*64*12,31);
    paramA_zstat_r=reshape(paramA_zstat,64*64*12,31);
    paramB_zstat_r=reshape(paramB_zstat,64*64*12,31);
    gm_cortex_r=reshape(gm_cortex,64*64*12,1);
    
    % *default_zstatr=reshape*  - takes zstat default map 
    %% Make Mask for gm_cortex
    default_zstat_vals=default_zstat_r(gm_cortex_r>0,:);
    paramA_zstat_vals=paramA_zstat_r(gm_cortex_r>0,:);
    paramB_zstat_vals=paramB_zstat_r(gm_cortex_r>0,:);
    
    %% Set up table inputs
    
    %rows={'default'; 'paramA'; 'paramB'};

    
    % zstatmeans takes means of num you get
    %zstatmeansnoZ(1,:)=mean(default_zstat_vals(:,default_zstat_vals~=0)); %not =/= to 0 - find for 0 andonly applying to one time point
    
    %% Create Mean and Median variables for table
   
%keyboard;
    %% Make Medians without 0's
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
%     plot(zstatmedians_final);
%     title(sprintf('zstatMedians No Zeros sub%s',subj))

    
    %% Medians WITH 0's 
    zstatmedians(1,:)=median(default_zstat_vals,1);
    zstatmedians(2,:)=median(paramA_zstat_vals,1);
    zstatmedians(3,:)=median(paramB_zstat_vals,1);
    
    %save to one part per iteration
    zstatmedians_final{k}=(zstatmedians)'
    dlmwrite(sprintf('%s/zstatmedians_final_sub-%s_%s.tsv', DataOutputPerSubj, subj, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmedians_final{k},'\t')
    
    %1x31 or 31x1 or by 3 - one for each MRIparam??
        %% Make Means without 0's
    default_zstat_valsNOZ=default_zstat_vals;
    default_zstat_valsNOZ(default_zstat_valsNOZ(:,1)==0,:)=NaN;
    zstatmeans_NOZ(1,:)= nanmean(default_zstat_valsNOZ)
    
    paramA_zstat_valsNOZ=paramA_zstat_vals;
    paramA_zstat_valsNOZ(paramA_zstat_valsNOZ(:,1)==0,:)=NaN;
    zstatmeans_NOZ(2,:)=nanmean(paramA_zstat_valsNOZ)
    
    paramB_zstat_valsNOZ=paramB_zstat_vals;
    paramB_zstat_valsNOZ(paramB_zstat_valsNOZ(:,1)==0,:)=NaN;
    zstatmeans_NOZ(3,:)=nanmean(paramB_zstat_valsNOZ)
    
    %save to zstatmeansNOZ for each iteration
    zstatmeansNOZ_final{k}=(zstatmeans_NOZ)'
    dlmwrite(sprintf('%s/zstatmeansNOZ_final_sub-%s_%s.tsv', DataOutputPerSubj, subj, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmeansNOZ_final{k},'\t')  
    %% Means WITH 0's
    %zstatmeans(1,:)=mean(default_zstat_vals,1); %default as row, added comma 2
%     zstatmeans_default(:,k)=mean(default_zstat_vals,1)';
%     zstatmeans_paramA(:,k)=mean(paramA_zstat_vals,1)';
%     zstatmeans_paramB(:,k)=mean(default_zstat_vals,1)';
    
    zstatmeans(1,:)=mean(default_zstat_vals,1);
    zstatmeans(2,:)=mean(paramA_zstat_vals,1);
    zstatmeans(3,:)=mean(paramB_zstat_vals,1);
    
    
    zstatmeans_final{k}=(zstatmeans)'
    dlmwrite(sprintf('%s/zstatmeans_final_sub-%s_%s.tsv', DataOutputPerSubj, subj, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmeans_final{k},'\t')
%     zstatmeans(2,:)=mean(paramA_zstat_vals,1);
%     zstatmeans(3,:)=mean(paramB_zstat_vals,1);
    
    %saveas(gcf,fullfile(srcout,['T' num2str(subj) 'table_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg']));
%keyboard;

end


%write out for total zstat files
%zstatmedians
dlmwrite(sprintf('%s/zstatmedians_final_all_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmedians_final,'\t')
%zstatmediansNOZ
dlmwrite(sprintf('%s/zstatmediansNOZ_final_all_%s.tsv', DataOutput,datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmediansNOZ_final,'\t')

%zstatmeans
dlmwrite(sprintf('%s/zstatmeans_final_all_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmeans_final,'\t')
%zstatmeansNOZ
dlmwrite(sprintf('%s/zstatmeansNOZ_final_all_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),zstatmeansNOZ_final,'\t')  

% for SD - weighted 
% Does same for medians And mediansnoZ 
%For all 3 parameters
%% Make Table
%T=table(rows,zstatmeans,zstatmeansnoZ,zstatmedians,zstatmediansnoZ); 
% T=table names of rows variables of the columns names - > in command window 
% you have a table w. headings of rows, zstatm, and numbers listed down the columns 
% - that table T is produced by the funcn 
% should have 32 columns - could look at zstatmediansnoZ - plot zstatmedianz
%% Make a new table here %%

%% TABLES %%
%zstatmeans
Tablezstatmeans = cell2mat(zstatmeans_final)
Headers = {'sub01default','sub01paramA','sub01paramB',...
'sub02default','sub02paramA','sub02paramB',...
'sub03default','sub03paramA','sub03paramB',...
'sub05default','sub05paramA','sub05paramB',...
'sub07default','sub07paramA','sub07paramB',...
'sub08default','sub08paramA','sub08paramB'...
'sub09default','sub09paramA','sub09paramB'};
Tablezstatmeans2 = [Headers; num2cell(Tablezstatmeans)]
Tablezstatmeans3 = cell2table(Tablezstatmeans2)
Tablezstatmeans4 = Tablezstatmeans2(2:end,:);
Table_zstatmeans5 = cell2table(Tablezstatmeans4)
Table_zstatmeans5.Properties.VariableNames = Tablezstatmeans2(1,:)
zstatmeans_final_table = Table_zstatmeans5


%sumzstatmeans=summary(zstatmeans_final_table)
zstatmeans_final_table.Overallmean = mean(zstatmeans_final_table{:,2:end},2)
%zstatmeans_final_table.eachtimepointmean = mean(zstatmeans_final_table{2:end,:})

%rows={-60:4:60}';

%zstatmeansNOZ
TablezstatmeansNOZ = cell2mat(zstatmeansNOZ_final)
Headers = {'sub01default','sub01paramA','sub01paramB',...
'sub02default','sub02paramA','sub02paramB',...
'sub03default','sub03paramA','sub03paramB',...
'sub05default','sub05paramA','sub05paramB',...
'sub07default','sub07paramA','sub07paramB',...
'sub08default','sub08paramA','sub08paramB'...
'sub09default','sub09paramA','sub09paramB'};
TablezstatmeansNOZ2 = [Headers; num2cell(TablezstatmeansNOZ)]
TablezstatmeansNOZ3 = cell2table(TablezstatmeansNOZ2)
TablezstatmeansNOZ4 = TablezstatmeansNOZ2(2:end,:);
Table_zstatmeansNOZ5 = cell2table(TablezstatmeansNOZ4)
Table_zstatmeansNOZ5.Properties.VariableNames = TablezstatmeansNOZ2(1,:)
zstatmeansNOZ_final_table = Table_zstatmeansNOZ5

%zstatmedians
Tablezstatmedians = cell2mat(zstatmedians_final)
Headers = {'sub01default','sub01paramA','sub01paramB',...
'sub02default','sub02paramA','sub02paramB',...
'sub03default','sub03paramA','sub03paramB',...
'sub05default','sub05paramA','sub05paramB',...
'sub07default','sub07paramA','sub07paramB',...
'sub08default','sub08paramA','sub08paramB'...
'sub09default','sub09paramA','sub09paramB'};
Tablezstatmedians2 = [Headers; num2cell(Tablezstatmedians)]
Tablezstatmedians3 = cell2table(Tablezstatmedians2)
Tablezstatmedians4 = Tablezstatmedians2(2:end,:);
Table_zstatmedians5 = cell2table(Tablezstatmedians4)
Table_zstatmedians5.Properties.VariableNames = Tablezstatmedians2(1,:)
zstatmedians_final_table = Table_zstatmedians5

%zstatmediansNOZ
TablezstatmediansNOZ = cell2mat(zstatmediansNOZ_final)
Headers = {'sub01default','sub01paramA','sub01paramB',...
'sub02default','sub02paramA','sub02paramB',...
'sub03default','sub03paramA','sub03paramB',...
'sub05default','sub05paramA','sub05paramB',...
'sub07default','sub07paramA','sub07paramB',...
'sub08default','sub08paramA','sub08paramB'...
'sub09default','sub09paramA','sub09paramB'};
TablezstatmediansNOZ2 = [Headers; num2cell(TablezstatmediansNOZ)]
TablezstatmediansNOZ3 = cell2table(TablezstatmediansNOZ2)
TablezstatmediansNOZ4 = TablezstatmediansNOZ2(2:end,:);
Table_zstatmediansNOZ5 = cell2table(TablezstatmediansNOZ4)
Table_zstatmediansNOZ5.Properties.VariableNames = TablezstatmediansNOZ2(1,:)
zstatmediansNOZ_final_table = Table_zstatmediansNOZ5

%% PLOT


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


figure
plot((-60:4:60),mean(Tablezstatmedians(:,1:3:end),2),'r')
hold on
plot((-60:4:60),mean(Tablezstatmedians(:,2:3:end),2),'g')
plot((-60:4:60),mean(Tablezstatmedians(:,3:3:end),2),'b')
title('Zstat Means of Medians');
grid
legend('Default','paramA','paramB')

figure
plot((-60:4:60),mean(Tablezstatmeans(:,1:3:end),2),'r')
hold on
plot((-60:4:60),mean(Tablezstatmeans(:,2:3:end),2),'g')
plot((-60:4:60),mean(Tablezstatmeans(:,3:3:end),2),'b')
title('Zstat Means of Means');
grid
legend('Default','paramA','paramB')
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


%% add new figures here for the means and medians here %%

% figure
% plot(zstatmeans_final)
% title('zstatmeans_final')
% figure;
% hist(default_zstat_vals(:,16),1000); %distribution for each time delay - 
%hist(default_zstat_vals,1000) %1000 bins - can specify what bins, could
%specify range of valsues e..g 1:1:10 to do in steps of 10 seconds - inside
%bins - then outside bins with everythign else - look
%
% title('zstat default');
% 
% figure;
% hist(paramA_zstat_vals,1000)
% title('zstat paramA');
% 
% figure;
% hist(paramB_zstat_vals,1000)
% title('zstat paramB');
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

fprintf('~~~ End of Script ~~~');