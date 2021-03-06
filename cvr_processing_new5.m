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
    srcout = '/Users/colette/sourcedata/derivatives/StatsOutput/'
    src = '/Users/colette/sourcedata/'
    %srcout = sprintf('/Users/colette/sourcedata/derivatives/sub-%s', subj)
    CVRmaps_dir = sprintf('%sCVRmaps_all', srcout)
    
    %% *read_avw*
    
    %default_cvr
    
    %default_cvr = '/Users/colette/sourcedata/derivatives/StatsOutput/CVRmaps_all/CVRmaps_concat_default.nii.gz'
    MRIParam1= 'default'
    [default_cvr, dims, scales, bpp, endian] = read_avw([CVRmaps_dir '/CVRmaps_concat_sub-' subj '_' MRIParam1 '.nii.gz']);
    
    MRIParam2 = 'paramA'
    [paramA_cvr, dims, scales, bpp, endian] = read_avw([CVRmaps_dir '/CVRmaps_concat_sub-' subj '_' MRIParam2 '.nii.gz']);
    
    MRIParam3 = 'paramB'
    [paramB_cvr, dims, scales, bpp, endian] = read_avw([CVRmaps_dir '/CVRmaps_concat_sub-' subj '_' MRIParam3 '.nii.gz']);
    
    % gm_cortex %
    %gm_cortex=read_avw(['/Users/colette/sourcedata/derivatives/sub-01/GM_MZeroScan.nii.gz']);
    gm_cortex = read_avw([sprintf('%sderivatives/sub-%s/sub-%s_GM_MZeroScan.nii.gz', src, subj, subj)]);
    %% Reads in the cortex map that we created GM_MZeroScan and assigns it to variable
    % called gm_cortex
    % *Input*: sourcedata folder/subjectnumber/func/subject_ GM_MZeroScan
    % *Output*: variable gm_cortex
    %%figure
    
%     figure; %creates a figure
%     imagesc(default_cvr(:,:,6),[0 2]); %displays default_cvr as a sanity check
%     colorbar; %colour bar displayed in figure
    
    %% reshape
    default_cvr_r=reshape(default_cvr,64*64*12,31); %skip number 000
    %default_cvr_r=reshape(default_cvr,64*64*12,31);
    paramA_cvr_r=reshape(paramA_cvr,64*64*12,31);
    paramB_cvr_r=reshape(paramB_cvr,64*64*12,31);
    gm_cortex_r=reshape(gm_cortex,64*64*12,1);
    
    % *default_cvrr=reshape*  - takes cvr default map 
    %% Make Mask for gm_cortex
    default_cvr_vals=default_cvr_r(gm_cortex_r>0,:);
    paramA_cvr_vals=paramA_cvr_r(gm_cortex_r>0,:);
    paramB_cvr_vals=paramB_cvr_r(gm_cortex_r>0,:);
    
    %% Set up table inputs
    
    %rows={'default'; 'paramA'; 'paramB'};

    
    % Cvrmeans takes means of num you get
    %cvrmeansnoZ(1,:)=mean(default_cvr_vals(:,default_cvr_vals~=0)); %not =/= to 0 - find for 0 andonly applying to one time point
    
    %% Create Mean and Median variables for table
   
%keyboard;
    %% Make Medians without 0's
    default_cvr_valsNOZ=default_cvr_vals;
    default_cvr_valsNOZ(default_cvr_valsNOZ(:,1)==0,:)=NaN;
    cvrmedians_NOZ(1,:)= nanmedian(default_cvr_valsNOZ)
    
    cvrmedians_NOZ(1,:)= nanmedian(default_cvr_valsNOZ)
    
    paramA_cvr_valsNOZ=paramA_cvr_vals;
    paramA_cvr_valsNOZ(paramA_cvr_valsNOZ(:,1)==0,:)=NaN;
    cvrmedians_NOZ(2,:)=nanmedian(paramA_cvr_valsNOZ)
    
    paramB_cvr_valsNOZ=paramB_cvr_vals;
    paramB_cvr_valsNOZ(paramB_cvr_valsNOZ(:,1)==0,:)=NaN;
    cvrmedians_NOZ(3,:)=nanmedian(paramB_cvr_valsNOZ)
    
    %save cvrmedians to one variable column per iteration
    cvrmediansNOZ_final{k}=(cvrmedians_NOZ)'
    dlmwrite(sprintf('%s/cvrmediansNOZ_final_sub-%s_%s.tsv', DataOutputPerSubj, subj, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),cvrmediansNOZ_final{k},'\t')
%     plot(cvrmedians_final);
%     title(sprintf('cvrMedians No Zeros sub%s',subj))

    
    %% Medians WITH 0's 
    cvrmedians(1,:)=median(default_cvr_vals,1);
    cvrmedians(2,:)=median(paramA_cvr_vals,1);
    cvrmedians(3,:)=median(paramB_cvr_vals,1);
    
    %save to one part per iteration
    cvrmedians_final{k}=(cvrmedians)'
    dlmwrite(sprintf('%s/cvrmedians_final_sub-%s_%s.tsv', DataOutputPerSubj, subj, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),cvrmedians_final{k},'\t')
    
    %1x31 or 31x1 or by 3 - one for each MRIparam??
        %% Make Means without 0's
    default_cvr_valsNOZ=default_cvr_vals;
    default_cvr_valsNOZ(default_cvr_valsNOZ(:,1)==0,:)=NaN;
    cvrmeans_NOZ(1,:)= nanmean(default_cvr_valsNOZ)
    
    paramA_cvr_valsNOZ=paramA_cvr_vals;
    paramA_cvr_valsNOZ(paramA_cvr_valsNOZ(:,1)==0,:)=NaN;
    cvrmeans_NOZ(2,:)=nanmean(paramA_cvr_valsNOZ)
    
    paramB_cvr_valsNOZ=paramB_cvr_vals;
    paramB_cvr_valsNOZ(paramB_cvr_valsNOZ(:,1)==0,:)=NaN;
    cvrmeans_NOZ(3,:)=nanmean(paramB_cvr_valsNOZ)
    
    %save to cvrmeansNOZ for each iteration
    cvrmeansNOZ_final{k}=(cvrmeans_NOZ)'
    dlmwrite(sprintf('%s/cvrmeansNOZ_final_sub-%s_%s.tsv', DataOutputPerSubj, subj, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),cvrmeansNOZ_final{k},'\t')  
    %% Means WITH 0's
    %cvrmeans(1,:)=mean(default_cvr_vals,1); %default as row, added comma 2
%     cvrmeans_default(:,k)=mean(default_cvr_vals,1)';
%     cvrmeans_paramA(:,k)=mean(paramA_cvr_vals,1)';
%     cvrmeans_paramB(:,k)=mean(default_cvr_vals,1)';
    
    cvrmeans(1,:)=mean(default_cvr_vals,1);
    cvrmeans(2,:)=mean(paramA_cvr_vals,1);
    cvrmeans(3,:)=mean(paramB_cvr_vals,1);
    
    
    cvrmeans_final{k}=(cvrmeans)'
    dlmwrite(sprintf('%s/cvrmeans_final_sub-%s_%s.tsv', DataOutputPerSubj, subj, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),cvrmeans_final{k},'\t')
%     cvrmeans(2,:)=mean(paramA_cvr_vals,1);
%     cvrmeans(3,:)=mean(paramB_cvr_vals,1);
    
    %saveas(gcf,fullfile(srcout,['T' num2str(subj) 'table_' datestr(clock,'yyyy-mm-dd_HH-MM-SS') '.jpg']));
%keyboard;

end


%write out for total cvr files
%cvrmedians
dlmwrite(sprintf('%s/cvrmedians_final_all_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),cvrmedians_final,'\t')
%cvrmediansNOZ
dlmwrite(sprintf('%s/cvrmediansNOZ_final_all_%s.tsv', DataOutput,datestr(clock,'yyyy-mm-dd_HH-MM-SS')),cvrmediansNOZ_final,'\t')

%cvrmeans
dlmwrite(sprintf('%s/cvrmeans_final_all_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),cvrmeans_final,'\t')
%cvrmeansNOZ
dlmwrite(sprintf('%s/cvrmeansNOZ_final_all_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),cvrmeansNOZ_final,'\t')  

% for SD - weighted 
% Does same for medians And mediansnoZ 
%For all 3 parameters
%% Make Table
%T=table(rows,cvrmeans,cvrmeansnoZ,cvrmedians,cvrmediansnoZ); 
% T=table names of rows variables of the columns names - > in command window 
% you have a table w. headings of rows, cvrm, and numbers listed down the columns 
% - that table T is produced by the funcn 
% should have 32 columns - could look at cvrmediansnoZ - plot cvrmedianz
%% Make a new table here %%

%% TABLES %%
%Cvrmeans
Tablecvrmeans = cell2mat(cvrmeans_final)
Headers = {'sub01default','sub01paramA','sub01paramB',...
'sub02default','sub02paramA','sub02paramB',...
'sub03default','sub03paramA','sub03paramB',...
'sub05default','sub05paramA','sub05paramB',...
'sub07default','sub07paramA','sub07paramB',...
'sub08default','sub08paramA','sub08paramB'...
'sub09default','sub09paramA','sub09paramB'};
Tablecvrmeans2 = [Headers; num2cell(Tablecvrmeans)]
Tablecvrmeans3 = cell2table(Tablecvrmeans2)
Tablecvrmeans4 = Tablecvrmeans2(2:end,:);
Table_cvrmeans5 = cell2table(Tablecvrmeans4)
Table_cvrmeans5.Properties.VariableNames = Tablecvrmeans2(1,:)
cvrmeans_final_table = Table_cvrmeans5

%cvrmeansNOZ
TablecvrmeansNOZ = cell2mat(cvrmeansNOZ_final)
Headers = {'sub01default','sub01paramA','sub01paramB',...
'sub02default','sub02paramA','sub02paramB',...
'sub03default','sub03paramA','sub03paramB',...
'sub05default','sub05paramA','sub05paramB',...
'sub07default','sub07paramA','sub07paramB',...
'sub08default','sub08paramA','sub08paramB'...
'sub09default','sub09paramA','sub09paramB'};
TablecvrmeansNOZ2 = [Headers; num2cell(TablecvrmeansNOZ)]
TablecvrmeansNOZ3 = cell2table(TablecvrmeansNOZ2)
TablecvrmeansNOZ4 = TablecvrmeansNOZ2(2:end,:);
Table_cvrmeansNOZ5 = cell2table(TablecvrmeansNOZ4)
Table_cvrmeansNOZ5.Properties.VariableNames = TablecvrmeansNOZ2(1,:)
cvrmeansNOZ_final_table = Table_cvrmeansNOZ5

%cvrmedians
Tablecvrmedians = cell2mat(cvrmedians_final)
Headers = {'sub01default','sub01paramA','sub01paramB',...
'sub02default','sub02paramA','sub02paramB',...
'sub03default','sub03paramA','sub03paramB',...
'sub05default','sub05paramA','sub05paramB',...
'sub07default','sub07paramA','sub07paramB',...
'sub08default','sub08paramA','sub08paramB'...
'sub09default','sub09paramA','sub09paramB'};
Tablecvrmedians2 = [Headers; num2cell(Tablecvrmedians)]
Tablecvrmedians3 = cell2table(Tablecvrmedians2)
Tablecvrmedians4 = Tablecvrmedians2(2:end,:);
Table_cvrmedians5 = cell2table(Tablecvrmedians4)
Table_cvrmedians5.Properties.VariableNames = Tablecvrmedians2(1,:)
cvrmedians_final_table = Table_cvrmedians5
%cvrmediansNOZ
TablecvrmediansNOZ = cell2mat(cvrmediansNOZ_final)
Headers = {'sub01default','sub01paramA','sub01paramB',...
'sub02default','sub02paramA','sub02paramB',...
'sub03default','sub03paramA','sub03paramB',...
'sub05default','sub05paramA','sub05paramB',...
'sub07default','sub07paramA','sub07paramB',...
'sub08default','sub08paramA','sub08paramB'...
'sub09default','sub09paramA','sub09paramB'};
TablecvrmediansNOZ2 = [Headers; num2cell(TablecvrmediansNOZ)]
TablecvrmediansNOZ3 = cell2table(TablecvrmediansNOZ2)
TablecvrmediansNOZ4 = TablecvrmediansNOZ2(2:end,:);
Table_cvrmediansNOZ5 = cell2table(TablecvrmediansNOZ4)
Table_cvrmediansNOZ5.Properties.VariableNames = TablecvrmediansNOZ2(1,:)
cvrmediansNOZ_final_table = Table_cvrmediansNOZ5

%dlmwrite(sprintf('%s/cvrmeans_final_all_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),Table_new,'\t')

% save('savefile.mat', 'T')
% save(sprintf('%s/cvrmeans_table_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),Table_new,'\t')
% save('cvrmeans_table.tsv', 'Table_new')
% dlmwrite(sprintf('%s/cvrmeans_table_%s.tsv', DataOutput, datestr(clock,'yyyy-mm-dd_HH-MM-SS')),Table_new,'\t') 
% 
% dlmwrite('cvrmeans_table',Table_new,'\t') 



%% Make loop for tables here: 
% cvrlist = [cvrmeans_final cvrmeansNOZ_final cvrmedians_final cvrmedians_final];
% for c= 1:4
%  %final = sprintf('%s_final',c)
%   final= cvrlist({c})
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
% plot(cvrmeans_final)
% title('cvrmeans_final')
% figure;
% hist(default_cvr_vals(:,16),1000); %distribution for each time delay - 
%hist(default_cvr_vals,1000) %1000 bins - can specify what bins, could
%specify range of valsues e..g 1:1:10 to do in steps of 10 seconds - inside
%bins - then outside bins with everythign else - look
%
% title('CVR default');
% 
% figure;
% hist(paramA_cvr_vals,1000)
% title('CVR paramA');
% 
% figure;
% hist(paramB_cvr_vals,1000)
% title('CVR paramB');
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