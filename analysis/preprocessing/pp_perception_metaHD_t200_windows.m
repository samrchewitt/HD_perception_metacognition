%Sam Hewitt, UCL 2020
%script to preprocess raw task data for each group
%provides accuracy, response, confidence and rt data, and summary data for
%behavuiural analysis
%also preprocess raw data for hddm analysis:
%removes dd=0 trials and takes the last 160 (80%) of trials to allow for staircase stabilisation
%if savedata=1, saves the data and csv file for further analysis
%update paths if using a different device
clear all
close all

%specify jobs: 
preprocess=1; %1=preprocess raw data or 0=load previously preprocessed data
savedata=1;

if preprocess
%---HC group:
hcpath = ('D:\meta-hd\meta_dots-master\perceptData\perceptDataHC');
addpath(hcpath); %file folder
cd(hcpath);
addpath('D:\meta-hd\meta_dots-master\perceptData\perceptDataHC\datafiles'); %data folder
files = dir('perceptData*.mat') ;

%parameters:
N = length(files); %number of ppts
Nratings = 6; %number of confidence ratings available

%loop to create confidence, accuracy, response matrix (N x 8 blocks)
for n = 1:N
    load(files(n).name) %load the first file
    for b=1:8
    conf{n, b} = (DATA(b).results.responseConf.*5) + 1; %get confidence ratings       
    acc{n, b} = DATA(b).results.correct; %get accuracy
    resp{n, b} = DATA(b).results.response; %get the response chosen
    rt{n, b} = DATA(b).results.rt; %get the rt of each response
    acc_blocks_hc(n, b)=sum((acc{n, b})/25)*100;
    end
end

%calculate the dot difference for each trial:
for n = 1:N %for n participants
    load(files(n).name) %load the first file
    for b=1:8 %blocks 1 to 8
    dots_s{n, b} = DATA(b).results.dots; %get the dot value structure (N x blocks)
    for t=1:25 %add the trials 1-25
    dotscell{n, b}(t) = struct2cell(dots_s{n, b}(t));
    for i=1:2
        dots{n, b}(t, i) = length(dotscell{n, b}{t}{i}); %length of dot cell on each trial = n of dots, add this to an array 25 x 2
    end 
    diff{n, b} = (dots{n, b}(:,2) - dots{n, b}(:,1));
    diff_abs_hc{n, b}=abs(diff{n, b}); %get absolute difference values 
    difference{n, b}=diff_abs_hc{n, b}'; %transpose difference values to 1 row x 25 trials (to match: acc, resp, conf)
    difference_r{n, b}=diff{n, b}'; %transpose difference values to 1 row x 25 trials (to match: acc, resp, conf)
    end
    end
end

%convert conf, accuracy, response and rt to arrays
confidence_hc= cell2mat(conf);
accuracy_hc=cell2mat(acc);
response_hc=cell2mat(resp);
rts_hc=cell2mat(rt);
difference_abs_hc=cell2mat(difference);
difference_r_hc=cell2mat(difference_r);
correct_hc=mean(acc_blocks_hc, 2);

%get ID numbers
for n=1:N
ppt_ids(n)=extractBetween((files(n).name),["perceptData"], [".mat"]); 
end
id_hc=transpose(str2double(ppt_ids)); %transpose IDs

%save the data
datapath='D:\meta-hd\meta_dots-master\perceptData\perceptDataHC\datafiles';
save(fullfile(datapath, 'id_hc'), 'id_hc');
save(fullfile(datapath, 'confidence_hc'), 'confidence_hc');
save(fullfile(datapath, 'accuracy_hc'), 'accuracy_hc');
save(fullfile(datapath, 'response_hc'), 'response_hc');
save(fullfile(datapath, 'rts_hc'), 'rts_hc');
save(fullfile(datapath, 'diff_hc'), 'diff_abs_hc');
save(fullfile(datapath, 'difference_abs_hc'), 'difference_abs_hc');
save(fullfile(datapath, 'difference_r_hc'), 'difference_r_hc');
save(fullfile(datapath, 'acc_blocks_hc'), 'acc_blocks_hc');
save(fullfile(datapath, 'correct_hc'), 'correct_hc');

%clear workspace
clearvars -except preprocess savedata

%---pre-HD group:
prehdpath = ('D:\meta-hd\meta_dots-master\perceptData\perceptDataHD\prehd (copy)');
addpath(prehdpath); %file folder
cd(prehdpath);
datapath=('D:\meta-hd\meta_dots-master\perceptData\perceptDataHD\prehd (copy)\datafiles'); %data folder
addpath(datapath);
files = dir('perceptData*.mat') ;

%parameters:
N = length(files) ;

%loop to create confidence, accuracy, response matrix (N x 8 blocks)
for n = 1:N
    load(files(n).name) %load the first file
    for b=1:8
    conf{n, b} = (DATA(b).results.responseConf.*5) + 1; %get confidence ratings       
    acc{n, b} = DATA(b).results.correct; %get accuracy
    resp{n, b} = DATA(b).results.response; %get the response chosen
    rt{n, b} = DATA(b).results.rt; %get the rt of each response
    acc_blocks_pre(n, b)=sum((acc{n, b})/25)*100;
    end
end

%calculate the % of correct responses for each ppt
for n=1:N
    for b=1:8
    correct{n, b} = (sum(acc{n, b})/25)*100; %percentage correct (n x block)
    mean_rts_prehd(n, b) = mean(rt{n, b}); %get the rt of each response
    end
end

%calcuate dot difference on each trial:
for n = 1:N %for n participants
    load(files(n).name) %load the first file
    for b=1:8 %blocks 1 to 8
    dots_s{n, b} = DATA(b).results.dots; %get the dot value structure (N x blocks)
    for t=1:25 %add the trials 1-25
    dotscell{n, b}(t) = struct2cell(dots_s{n, b}(t));
    for i=1:2
        dots{n, b}(t, i) = length(dotscell{n, b}{t}{i}); %length of dot cell on each trial = n of dots, add this to an array 25 x 2

    end 
    diff{n, b} = (dots{n, b}(:,2) - dots{n, b}(:,1));
    diff_abs_pre{n, b}=abs(diff{n, b}); %get absolute difference
    difference{n, b}=diff_abs_pre{n, b}'; %transpose difference values to 1 row x 25 trials (to match: acc, resp, conf)
    difference_r{n, b}=diff{n, b}'; %transpose difference values to 1 row x 25 trials (to match: acc, resp, conf)
    end
    end
end

%convert conf, accuracy, response and rt to arrays to save
confidence_pre= cell2mat(conf);
accuracy_pre=cell2mat(acc);
response_pre=cell2mat(resp);
rts_pre=cell2mat(rt);
difference_abs_pre=cell2mat(difference);
difference_r_pre=cell2mat(difference_r);
correct_pre=mean(acc_blocks_pre, 2); %get mean correct per patient

%get ID numbers
for n=1:N
ppt_ids(n)=extractBetween((files(n).name),["perceptData"], [".mat"]); 
end
id_pre=transpose(str2double(ppt_ids)); %transpose IDs

%save the data
save(fullfile(datapath, 'id_prehd'), 'id_pre');
save(fullfile(datapath, 'confidence_pre'), 'confidence_pre');
save(fullfile(datapath, 'accuracy_pre'), 'accuracy_pre');
save(fullfile(datapath, 'response_pre'), 'response_pre');
save(fullfile(datapath, 'rts_pre'), 'rts_pre');
save(fullfile(datapath, 'diff_pre'), 'diff_abs_pre');
save(fullfile(datapath, 'difference_abs_pre'), 'difference_abs_pre');
save(fullfile(datapath, 'difference_r_pre'), 'difference_r_pre');
save(fullfile(datapath, 'acc_blocks_pre'), 'acc_blocks_pre');
save(fullfile(datapath, 'correct_pre'), 'correct_pre');

%clear workspace
clearvars -except preprocess savedata

%-----early-HD group:
earlyhdpath = ('D:\meta-hd\meta_dots-master\perceptData\perceptDataHD\earlyhd (copy)');
addpath(earlyhdpath); %file folder
cd(earlyhdpath);
datapath=('D:\meta-hd\meta_dots-master\perceptData\perceptDataHD\earlyhd (copy)\datafiles'); %data folder
addpath(datapath);
files = dir('perceptData*.mat') ;

%parameters:
N = length(files) ;

%loop to create confidence, accuracy, response matrix (N x 8 blocks)
for n = 1:N
    load(files(n).name) %load the first file
    for b=1:8
    conf{n, b} = (DATA(b).results.responseConf.*5) + 1; %get confidence ratings       
    acc{n, b} = DATA(b).results.correct; %get accuracy
    resp{n, b} = DATA(b).results.response; %get the response chosen
    rt{n, b} = DATA(b).results.rt; %get the rt of each response
    acc_blocks_early(n, b)=sum((acc{n, b})/25)*100;
    end
end

%calculate the % of correct responses for each ppt
for n=1:N
    for b=1:8
    correct{n, b} = (sum(acc{n, b})/25)*100; %percentage correct (n x block)
    mean_rts_earlyhd(n, b) = mean(rt{n, b}); %get the rt of each response
    end
end

for n = 1:N %for n participants
    load(files(n).name) %load the first file
    for b=1:8 %blocks 1 to 8
    dots_s{n, b} = DATA(b).results.dots; %get the dot value structure (N x blocks)
    for t=1:25 %add the trials 1-25
    dotscell{n, b}(t) = struct2cell(dots_s{n, b}(t));
    for i=1:2
        dots{n, b}(t, i) = length(dotscell{n, b}{t}{i}); %length of dot cell on each trial = n of dots, add this to an array 25 x 2

    end 
    diff{n, b} = (dots{n, b}(:,2) - dots{n, b}(:,1));
    diff_abs_early{n, b}=abs(diff{n, b}); %get absolute difference
    %n.b. making all diff values positive allows for mean calculation
    difference{n, b}=diff_abs_early{n, b}'; %transpose difference values to 1 row x 25 trials (to match: acc, resp, conf)
    difference_r{n, b}=diff{n, b}'; %transpose difference values to 1 row x 25 trials (to match: acc, resp, conf)
    end
    end
end

%convert conf, accuracy, response and rt to arrays to save
confidence_early= cell2mat(conf);
accuracy_early=cell2mat(acc);
response_early=cell2mat(resp);
rts_early=cell2mat(rt);
difference_abs_early=cell2mat(difference);
difference_r_early=cell2mat(difference_r);
correct_early=mean(acc_blocks_early, 2); %get mean correct per patient

%get ID numbers
for n=1:N
ppt_ids(n)=extractBetween((files(n).name),["perceptData"], [".mat"]); 
end
id_early=transpose(str2double(ppt_ids)); %transpose IDs

%save the data
save(fullfile(datapath, 'id_earlyhd'), 'id_early');
save(fullfile(datapath, 'confidence_early'), 'confidence_early');
save(fullfile(datapath, 'accuracy_early'), 'accuracy_early');
save(fullfile(datapath, 'response_early'), 'response_early');
save(fullfile(datapath, 'rts_early'), 'rts_early');
save(fullfile(datapath, 'diff_early'), 'diff_abs_early');
save(fullfile(datapath, 'difference_abs_early'), 'difference_abs_early');
save(fullfile(datapath, 'difference_r_early'), 'difference_r_early');
save(fullfile(datapath, 'acc_blocks_early'), 'acc_blocks_early');
save(fullfile(datapath, 'correct_early'), 'correct_early');

end

%clear workspace
clearvars -except preprocess savedata

%reload necessary files below:
%add datapaths:
addpath('D:\meta-hd\meta_dots-master\perceptData\perceptDataHC\datafiles'); %data folder
addpath('D:\meta-hd\meta_dots-master\perceptData\perceptDataHD\prehd (copy)\datafiles'); %data folder
addpath('D:\meta-hd\meta_dots-master\perceptData\perceptDataHD\earlyhd (copy)\datafiles'); %data folder

%and load data: 
%hc:
load('accuracy_hc.mat');
load('confidence_hc.mat');
load('response_hc.mat');
load('rts_hc.mat');
load('diff_hc.mat');
load('difference_abs_hc.mat');
load('difference_r_hc.mat');
load('id_hc.mat')
%prehd:
load('accuracy_pre.mat');
load('confidence_pre.mat');
load('response_pre.mat');
load('rts_pre.mat');
load('diff_pre.mat');
load('difference_abs_pre.mat');
load('difference_r_pre.mat');
load('id_prehd.mat')
%earlyhd:
load('accuracy_early.mat');
load('confidence_early.mat');
load('response_early.mat');
load('rts_early.mat');
load('diff_abs_early.mat');
load('difference_abs_early.mat');
load('difference_r_early.mat');
load('id_earlyhd.mat')

%concatenate the arrays together [hc;pre;early]:
diff_abs=[diff_abs_hc; diff_abs_pre; diff_abs_early];
difference_abs=[difference_abs_hc; difference_abs_pre; difference_abs_early];
difference_r=[difference_r_hc; difference_r_pre; difference_r_early];
response=[response_hc; response_pre; response_early];
accuracy=[accuracy_hc; accuracy_pre; accuracy_early];
confidence=[confidence_hc; confidence_pre; confidence_early];
rts=[rts_hc; rts_pre; rts_early];
%loop to create dd based on accuracy (neg=incorrect, pos=correct) 
N=58; %ppts
t=200; %trials
for n = 1:N %for n participants
    for i=1:t %trials 
        if accuracy(n, i)==0; %if accuracy==0:
            rdd_acc(n, i)=-difference_abs(n, i); %rdd_acc = negative dd
        else %aka if accuracy=1:
            rdd_acc(n, i)=difference_abs(n, i);
        end
    end
end
%zscore dd:
%z score within subjects (mean by row using DIM=2, n.b. 0 means Z score will use the sample sdev):
rdd_accZws=zscore(rdd_acc, 0, 2);
dd_absZws=zscore(difference_abs, 0, 2);
%z score between subjects (by all):
rdd_accZbs=zscore(rdd_acc, 0, 'all');
dd_absZbs=zscore(difference_abs, 0, 'all');

%-----create HDDM datafile csv 

%convert trials to long vector:
dd_long = reshape(difference_abs',[],1);
response_long = reshape(response',[],1)-1; %-1 to make 0(L) and 1(R)
accuracy_long = reshape(accuracy',[],1);
dd_r_long = reshape(difference_r',[],1);
rt_long = reshape(rts',[],1);
rdd_acc_long = reshape(rdd_acc',[],1);
rdd_accZws_long = reshape(rdd_accZws',[],1);
rdd_accZbs_long = reshape(rdd_accZbs',[],1);
dd_absZws_long = reshape(dd_absZws',[],1);
dd_absZbs_long = reshape(dd_absZbs',[],1);

%visualise to understand which adjusted DD is best to use:
figure;
histogram(rdd_acc_long)
title('RDD (accuracy)')
figure;
histogram(rdd_accZws_long)
title('RDD acc Z-score WS')
figure;
histogram(rdd_accZbs_long)
title('RDD acc Z-score BS')
figure;
histogram(dd_absZws_long)
title('DD (abs) Z-score WS')
figure;
histogram(dd_absZbs_long)
title('DD (abs) Z-score BS')
figure;
histogram(dd_long)
title('DD (abs) long');

%loop to create confidence, accuracy, response matrix (N x 8 blocks)
tablelength=t*N; %table length: 200 trials x n
ids=[id_hc; id_pre; id_early];
stimid_long=NaN(tablelength, 1); %create stimID array
%loop to assign stimID (correct stimulus left or right), based on accuracy and response
for r=1:tablelength %iterate for each row in the table
if accuracy_long(r, 1)==1; %if accuracy==1:
    stimid_long(r, 1)=response_long(r, 1); %then stimid==response, %-1 to make 0(left) and 1(right)
else %aka if accuracy=0:
    if response_long(r, 1)==0 %given acc=0, if resp=1, stimid=2
        stimid_long(r, 1)=1;
    else %else, given acc=0, resp must be 2, therefore stimid==1
        stimid_long(r, 1)=0;
    end
end
end

%make ID list
pptSindex(1:N)=NaN;
pptEindex(1:N)=NaN;
for s=1:N-1
pptSindex(1)=1; %first index = 1
pptSindex(s+1)=(t*s)+1; %get the starting index of each additional ppt
end
for e=1:N
pptEindex(e)=(t*e); %get the ending index of each ppt
end
for i=1:n
id_long((pptSindex(i):pptEindex(i)), 1)=ids(i); %create id_long for table
end
%check the ids are right:
xx = unique(id_long);       % temp vector of vals
x = sort(id_long);          % sorted input aligns with temp (lowest to highest)
f = zeros(size(xx)); % vector for freqs
% frequency for each value
for i = 1:length(xx)
    f(i) = sum(x == xx(i));
end
idcheck=all(f==160); %1=yes;

%define group variable:
group_hc=zeros(length(id_hc)*t, 1); %add group var HC=0
group_pre=ones(length(id_pre)*t,1); %pre=1
group_early=ones(length(id_early)*t,1)+1; %early=2
group=[group_hc; group_pre; group_early];
%make table:
hddmdata= table(id_long, accuracy_long, dd_long, rdd_acc_long, rdd_accZws_long, rdd_accZbs_long, dd_absZws_long, rt_long, group);
hddmdata.Properties.VariableNames={'subj_idx', 'response', 'dd', 'rdd', 'rdd_Zws', 'rdd_Zbs', 'dd_Zws', 'rt', 'group'}; %simplify column names

%---Remove outliers:
median_rt=median(hddmdata.rt); %get overall mean
%calculate median absolute deviation:
mediandev = mad(hddmdata.rt,1,'all')
%remove low rt (<0.1s) & high (> mean + 3*std) rts
mean_rt=mean(hddmdata.rt); %get overall mean
sdev_rt=std(hddmdata.rt);
outlier_rts=find((hddmdata.rt<0.1) | (hddmdata.rt > mean_rt+(3*std(hddmdata.rt))));
outlier_rts_MAD=find((hddmdata.rt<0.1) | (hddmdata.rt > median_rt+(3*mediandev)));

%percentage of trials are removed? 
outlier_pct=(length(outlier_rts_MAD)./height(hddmdata))*100;
hddm_outliers=hddmdata(outlier_rts,:); %save the outlier trials as separate table
hddm_outliers_MAD=hddmdata(outlier_rts_MAD,:); %save the outlier trials as separate table

%remove outliers from original table:
hddmdata(outlier_rts_MAD,:)=[];

%calculate significant differences in groups for excluded data?:
%observed proportion of excluded rts per group:
po_hc=sum(hddm_outliers_MAD.group==0)/height(hddm_outliers_MAD); 
po_pre=sum(hddm_outliers_MAD.group==1)/height(hddm_outliers_MAD);  
po_early=sum(hddm_outliers_MAD.group==2)/height(hddm_outliers_MAD); 
%expected proportion given group n
n=length(ids); %total n
pe_hc=length(id_hc)/n; 
pe_pre=length(id_pre)/n;
pe_early=length(id_early)/n;
observed=[po_hc po_pre po_early];
expected = [pe_hc pe_pre pe_early];
%calculate chi2stat:
chi2stat = sum((observed-expected).^2 ./ expected);
deg_freedom =2; %deg of freedom = n groups-1;
pchisq = 1-chi2cdf(chi2stat, deg_freedom);

%calculate summary statistic variables for behavioural analysis:
%correct (%), mean rts, mean dotdiff, median dotdiff, confidence 
%number of trials:
t=200;
%hc:
for n=1:length(id_hc)
    correct_hc_mean(:, n) = (sum(accuracy_hc(n, :))/t)*100; %percentage correct (1 x n) transposed
    rts_hc_mean(:, n) = mean(rts_hc(n, :)); %mean rt for each pt (transposed)
    dd_hc_mean(:, n)= mean(difference_abs_hc(n, :));
    dd_hc_median(:, n)= median(difference_abs_hc(n, :));
    conf_hc_mean(:, n) = mean(confidence_hc(n, :), 'OmitNaN'); %mean rt for each pt (transposed)
end
%pre:
for n=1:length(id_pre)
    correct_pre_mean(:, n) = (sum(accuracy_pre(n, :))/t)*100; %percentage correct (1 x n) transposed
    rts_pre_mean(:, n) = mean(rts_pre(n, :)); %mean rt for each pt (transposed)
    dd_pre_mean(:, n)= mean(difference_abs_pre(n, :));
    dd_pre_median(:, n)= median(difference_abs_pre(n, :));
    conf_pre_mean(:, n) = mean(confidence_pre(n, :), 'OmitNaN'); %mean rt for each pt (transposed)
end
%early:
for n=1:length(id_early)
    correct_early_mean(:, n) = (sum(accuracy_early(n, :))/t)*100; %percentage correct (1 x n) transposed
    rts_early_mean(:, n) = mean(rts_early(n, :)); %mean rt for each pt (transposed)
    dd_early_mean(:, n)= mean(difference_abs_early(n, :));
    dd_early_median(:, n)= median(difference_abs_early(n, :));
    conf_early_mean(:, n) = mean(confidence_early(n, :), 'OmitNaN'); %mean rt for each pt (transposed)

end

if savedata
%save the data:
%hc:
datapath_hc=('D:\meta-hd\meta_dots-master\perceptData\perceptDataHC\datafiles');
save(fullfile(datapath_hc, 'dd_hc_mean'), 'dd_hc_mean');
save(fullfile(datapath_hc, 'dd_hc_median'), 'dd_hc_median');
save(fullfile(datapath_hc, 'rts_hc_mean'), 'rts_hc_mean');
save(fullfile(datapath_hc, 'correct_hc_mean'), 'correct_hc_mean');
save(fullfile(datapath_hc, 'conf_hc_mean'), 'conf_hc_mean');

%pre
datapath_pre=('D:\meta-hd\meta_dots-master\perceptData\perceptDataHD\prehd (copy)\datafiles');
save(fullfile(datapath_pre, 'dd_pre_mean'), 'dd_pre_mean');
save(fullfile(datapath_pre, 'dd_pre_median'), 'dd_pre_median');
save(fullfile(datapath_pre, 'rts_pre_mean'), 'rts_pre_mean');
save(fullfile(datapath_pre, 'correct_pre_mean'), 'correct_pre_mean');
save(fullfile(datapath_pre, 'conf_pre_mean'), 'conf_pre_mean');

%early
datapath_early=('D:\meta-hd\meta_dots-master\perceptData\perceptDataHD\earlyhd (copy)\datafiles');
save(fullfile(datapath_early, 'dd_early_mean'), 'dd_early_mean');
save(fullfile(datapath_early, 'dd_early_median'), 'dd_early_median');
save(fullfile(datapath_early, 'rts_early_mean'), 'rts_early_mean');
save(fullfile(datapath_early, 'correct_early_mean'), 'correct_early_mean');
save(fullfile(datapath_early, 'conf_early_mean'), 'conf_early_mean');

%write hddm data to CSV
%writetable(hddmdata, '/Volumes/TOSHIBA EXT/meta-hd/meta_dots-master/perceptData/hddm csv/hddmdata_hlex_n58_MAD_Acoding.csv');

end

%check the id numbers again in table:
xx = unique(hddm_outliers_MAD.subj_idx);       % temp vector of vals
x = sort(hddm_outliers_MAD.subj_idx);          % sorted input aligns with temp (lowest to highest)
f = zeros(size(xx)); % vector for freqs
% frequency for each value
for i = 1:length(xx)
    f(i) = sum(x == xx(i));
end
f(:,2)=sort(xx);

%idcheck=all(f==160); %1=yes;