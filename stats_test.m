

% load table

array_flowchange_edit= '/Users/colette/sourcedata/derivatives/array_flowchange_edit.csv'

array_flowchange_edit=readmatrix(array_flowchange_edit);
array_flowchange_edit=load(array_flowchange_edit);

array_flowchange_edit=readcsv(/Users/colette/sourcedata/derivatives/array_flowchange_edit.csv);



x1 = [3;2;6;5];
x2 = [4;1;3;6;9];
x3 = [2;6;1];
%% 
x1=[6.0516;	6.241;	8.1903;	4.124;	5.6317;	8.0258;	7.6941;	4.3068;	5.1381;	5.7263;	4.2767;	4.6715;	7.6979;	1.3156]
x2= [2.7959; 4.7036; 2.0738; 2.5919; 4.2842; 5.8807; 6.4422; 5.4161; 1.895;	3.8221;	3.0623;	3.4254;	5.3338;	2.0903]
x3=[4.993;	4.3463;	2.9947;	5.4513;	5.056;	6.3264;	6.933;	4.2107;	6.5163;	5.5372;	3.0251;	4.0347;	4.4336]

data = [x1' x2' x3']; %// Create row vector with your data
group = {'Default', 'Default','Default','Default','Default','Default','Default','Default','Default','Default','Default','Default','Default','Default', ...
    'ParamA','ParamA','ParamA','ParamA','ParamA','ParamA','ParamA','ParamA','ParamA','ParamA','ParamA','ParamA','ParamA','ParamA' ...
    'ParamB','ParamB','ParamB','ParamB','ParamB','ParamB','ParamB','ParamB','ParamB','ParamB','ParamB','ParamB','ParamB'}; %// set the groups according to the data above

[p1] = anova1(data, group,'off') %// Use the 'off' option to prevent the table/box plot from showing up.

p1 =
%%

[p tbl sts]=anova1(arrayflowchangeedit(:,1:8)')
multcompare(sts)

