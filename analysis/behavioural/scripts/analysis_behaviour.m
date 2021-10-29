%behavioural data analysis script
%Sam Hewitt, UCL 2020
clear all
close all
%set default figure background to white:
set(0, 'defaultfigurecolor', 'w');
%define directories:
addpath('D:\matlab\Tools-master\plotting') %plotting functions
addpath('D:\matlab'); addpath('D:\matlab\Tools-master\');
datapath=['C:\Users\samrc\OneDrive\Documents\GitHub\samrchewitt\HD_perception_metacognition\analysis\behavioural\summary_data\']; %data folder
addpath(genpath(datapath));

%load required:
load('dd_hc_mean.mat')
load('dd_pre_mean.mat')
load('dd_early_mean.mat')
load('correct_pre_mean.mat')
load('correct_early_mean.mat')
load('correct_hc_mean.mat')

%get number of ppts:
n=length(dd_early_mean)+length(dd_pre_mean)+length(dd_hc_mean);

%--Figure 2a:
% 1: plot bar (mean + sem) and scatter (subject mean) values of accuracy (% correct)  (premanifest vs. early vs. hc)
% 2: conduct parametric assumptions checks
% 3: conduct appropriate statistical test and plot significance 

means(:,1)=mean(correct_hc_mean);%1 = HC
means(:,2)=mean(correct_pre_mean);%2 = preHD
means(:,3)=mean(correct_early_mean); %3 = earlyHD
sdev(:,1) = std(correct_hc_mean, 0, 2); %stdev HC
sem(:,1) = sdev(:,1) / sqrt(length(correct_hc_mean)) %sem HC
sdev(:,2) = std(correct_pre_mean, 0, 2); %stdev pre
sem(:,2) = sdev(:,2) / sqrt(length(correct_pre_mean)) %sem pre
sdev(:,3) = std(correct_early_mean, 0, 2); %stdev early HD
sem(:,3) = sdev(:,3) / sqrt(length(correct_early_mean)) %sem early HD

%plot the bars, errors=sem
t=tiledlayout(2,2) % use tiled layour to plot together
figure('PaperType', 'a4');
nexttile
b = bar([1, 2, 3], means, 0.75)
b.FaceColor = 'flat';
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [1 1 1];
b.CData(3,:) = [0.85 0.85 0.85];
hold on
errorbar([1, 2, 3], means, sem,'.', 'color', 'black', 'lineWidth', 1.2, 'CapSize', 30)
set(gca, 'YLim', [50 100], 'FontSize', 14); %larger font

%plot the individual data points (scatter)
ax1=gca;
hold(ax1, 'all')
x1(1:length(correct_hc_mean))=1;
x2(1:length(correct_pre_mean))=2;
x3(1:length(correct_early_mean))=3;
scatter(x1,correct_hc_mean, 30, 'black')
hold on
scatter(x2, correct_pre_mean, 30, 'black');
hold on
scatter(x3, correct_early_mean, 30, 'black');

%test equality of variance for KW vs. ANOVA:
%normal distribution?:
norm{:,1}=normalitytest(correct_hc_mean);
norm{:,2}=normalitytest(correct_pre_mean);
norm{:,3}=normalitytest(correct_early_mean);
for g=1:3
    if mean(norm{:,g}(:,3))==1
        disp(sprintf('Data from group %d comes from normal distribution', g));
    else
        disp(sprintf('Data from group %d does NOT come from normal distribution', g));
    end
 end
%initialise matrix:
corr(1:length(correct_hc_mean),1:3)=NaN;
corr(:,1)=correct_hc_mean;
corr(1:length(correct_pre_mean),2)=correct_pre_mean;
corr(1:length(correct_early_mean),3)=correct_early_mean;
[bt_p_corr]=vartestn(corr, 'display', 'off') %bartlett's 

 %structure data for one-way comparison
data=[correct_hc_mean, correct_pre_mean, correct_early_mean]
[group{1:length(correct_hc_mean)}]=deal('hc');
[group{end+1:end+length(correct_pre_mean)}]=deal('pre');
[group{end+1:n}]=deal('early');

if bt_p_corr>0.05
    disp('Groups have homogeneity of variance, conducting ANOVA...')
    disp('Find results in pval, tbl, stats and c for multicompare. pval_exact contains the exact sig value for 3-way comparison');
    %one way ANOVA:
[pval_c,tbl_c,stats_c] = anova1(data, group, 'off')
else
     disp('NO homogeneity of variance, conducting Kruskall Wallis');
    %kruskall wallis group comparison
[pval_c,tbl_c,stats_c] = kruskalwallis(data, group, 'off')
end

%display the test stats
%calculate partial eta squared effect size:
n2_c = round(rdivide(tbl_c{2, 2},tbl_c{4,2}), 2, 'significant');

if bt_p_corr>0.05
str=['{\it η^2} = ', num2str(n2_c), ',{\it p} = ', num2str(round(pval_c, 2, 'significant'))];
T = text(min(get(gca, 'xlim')+0.02), max(get(gca, 'ylim')-0.2), str); 
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
else
    str=['{\it H}', '(', num2str(tbl_c{2, 3}), ') = ', num2str(round(tbl_c{2, 5}, 2)), newline '{\it p} = ' num2str(round(vpa(pval_c), 8))];
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
end
xticklabels({'Control', 'Pre-HD', 'Early-HD'})
ylabel('% correct')
title('(A) Accuracy')
box off

%--Figure 2b:
% 1: plot bar (mean + sem) and scatter (subject mean) values of dot difference  (premanifest vs. early vs. hc)
% 2: conduct parametric assumptions checks 
% 3: conduct appropriate statistical test and plot significance 

means(:,1)=mean(dd_hc_mean); %1 = HC
means(:,2)=mean(dd_pre_mean);%1 = pre HD 
means(:,3)=mean(dd_early_mean);%2 = early HD
sdev(:,1) = std(dd_hc_mean, 0, 2); %stdev HC
sem(:,1) = sdev(:,1) / sqrt(length(dd_hc_mean)) %sem HC
sdev(:,2) = std(dd_pre_mean, 0, 2); %stdev pre
sem(:,2) = sdev(:,2) / sqrt(length(dd_pre_mean)) %sem pre
sdev(:,3) = std(dd_early_mean, 0, 2); %stdev early
sem(:,3) = sdev(:,3) / sqrt(length(dd_early_mean)) %sem early

%plot the bars, error bars=sem
nexttile
b = bar([1, 2, 3], means, 0.75)
b.FaceColor = 'flat';
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [1 1 1];
b.CData(3,:) = [0.85 0.85 0.85];
hold on
errorbar([1, 2, 3], means, sem,'.', 'color', 'black', 'lineWidth', 1.2, 'CapSize', 30)

%test equality of variance for KW vs. ANOVA:
%normal distribution?:
norm{:,1}=normalitytest(dd_hc_mean);
norm{:,2}=normalitytest(dd_pre_mean);
norm{:,3}=normalitytest(dd_early_mean);
for g=1:3
    if mean(norm{:,g}(:,3))==1
        disp(sprintf('Data from group %d comes from normal distribution', g));
    else
        disp(sprintf('Data from group %d does NOT come from normal distribution', g));
    end
end
%initialise matrix:
dd(1:length(dd_hc_mean),1:3)=NaN;
dd(:,1)=dd_hc_mean;
dd(1:length(dd_pre_mean),2)=dd_pre_mean;
dd(1:length(dd_early_mean),3)=dd_early_mean;
[bt_p]=vartestn(dd, 'display', 'off') %bartlett's test as data is normal
%structure data for one-way comparison
data=[dd_hc_mean, dd_pre_mean, dd_early_mean]

if bt_p>0.05
    disp('Groups have homogeneity of variance, conducting ANOVA...')
    disp('Find results in pval, tbl, stats and c for multicompare. pval_exact contains the exact sig value for 3-way comparison');
%one way ANOVA:
format long g
[pval_dd,tbl_dd,stats_dd] = anova1(data, group, 'off')

else
     disp('Groups have homogeneity of variance, conducting Kruskall Wallis');
    %kruskall wallis group comparison
[pval_dd,tbl_dd,stats_dd] = kruskalwallis(data, group, 'off')
end
if pval_dd < 0.001
    pval_dd_round='< 0.001'
end
%plot the individual data points (scatter)
ax1=gca;
hold(ax1, 'all')
x1(1:length(dd_hc_mean))=1;
x2(1:length(dd_pre_mean))=2;
x3(1:length(dd_early_mean))=3;
scatter(x1,dd_hc_mean, 30, 'black')
hold on
scatter(x2, dd_pre_mean, 30, 'black');
hold on
scatter(x3, dd_early_mean, 30, 'black');
set(gca, 'YLim', [2 13], 'FontSize', 14);
hold on
%display the test stats
%calculate eta squared effect size
n2_dd = round(rdivide(tbl_dd{2, 2},tbl_dd{4,2}), 2, 'significant');

if bt_p>0.05
str=['{\it η^2} =', num2str(n2_dd), ',{\it p} < 0.001'];
T = text(min(get(gca, 'xlim')+0.02), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
else
    str=['{\it H}', '(', num2str(tbl_dd{2, 3}), ') = ', num2str(round(tbl_dd{2, 5}, 2)), newline '{\it p} < 0.001'];
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
end
[c_dd, m_dd] = multcompare(stats_dd, 'display', 'off', 'CType', 'bonferroni')
%add sig stars
mysigstar(gca, [1,2], 10, c_dd(1, 6));
mysigstar(gca, [1,3], 11.5, c_dd(2, 6));
mysigstar(gca, [2,3], 10.8, c_dd(3, 6));

%x and y labels
ylabel('\Delta dots')
xticklabels({'Control', 'Pre-HD', 'Early-HD'})
title('(B) Stimulus strength')
box off
%c contains summary data, including corrected pvalues
%pval_exact is the precise pval for 3-way comparison

%Figure 2c:
%load mean values for rt 
%plot overall mean as bar chart with individual scatter
%%conduct parametric assumptions check and 3-way statistical comparison

load('rts_pre_mean.mat')
load('rts_early_mean.mat')
load('rts_hc_mean.mat')

means(:,1)=mean(rts_hc_mean); %1 = hc
means(:,2)=mean(rts_pre_mean); % 2=pre
means(:,3)=mean(rts_early_mean); %3 = early
sdev(:,1) = std(rts_hc_mean, 0, 2); %stdev HC
sem(:,1) = sdev(:,1) / sqrt(length(rts_hc_mean)) %sem HC
sdev(:,2) = std(rts_pre_mean, 0, 2); %stdev pre-HD
sem(:,2) = sdev(:,2) / sqrt(length(rts_pre_mean)) %sem pre-HD
sdev(:,3) = std(rts_early_mean, 0, 2); %stdev early
sem(:,3) = sdev(:,3) / sqrt(length(rts_early_mean)) %sem early

%plot the bars
nexttile
b = bar([1, 2, 3], means, 0.75)
b.FaceColor = 'flat';
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [1 1 1];
b.CData(3,:) = [0.85 0.85 0.85];
hold on
errorbar([1, 2, 3], means, sem,'.', 'color', 'black', 'lineWidth', 1.2, 'CapSize', 30)

%test equality of variance for KW vs. ANOVA:
%normal distribution?:
norm{:,1}=normalitytest(rts_hc_mean);
norm{:,2}=normalitytest(rts_pre_mean);
norm{:,3}=normalitytest(rts_early_mean);

%structure data for one-way comparison
data=[rts_hc_mean, rts_pre_mean, rts_early_mean];

for g=1:3
    if mean(norm{:,g}(:,3))==1
        disp(sprintf('Data from group %d comes from normal distribution', g));
        %initialise matrix:
rt(1:length(rts_hc_mean),1:3)=NaN;
rt(:,1)=rts_hc_mean;
rt(1:length(rts_pre_mean),2)=rts_pre_mean;
rt(1:length(rts_early_mean),3)=rts_early_mean;
[bt_p]=vartestn(rt, 'display', 'off'); %bartlett's test as data is normal
    else
        disp(sprintf('Data from group %d does NOT come from normal distribution', g));
        bt_p=NaN;
    end
end

if bt_p>0.05
    disp('Groups have homogeneity of variance, conducting ANOVA...');
    disp('Find results in stats_rt, tbl_rt, pval_rt, c_rt and m_rt (for post-hoc multicompare)');

    %one way ANOVA:
    [pval_rt,tbl_rt,stats_rt] = anova1(data, group, 'off');

else
     disp('Groups have homogeneity of variance, conducting Kruskall Wallis...');
    disp('Find results in stats_rt, tbl_rt, pval_rt, c_rt and m_rt (for post-hoc multicompare)');

     %kruskall wallis group comparison
[pval_rt,tbl_rt,stats_rt] = kruskalwallis(data, group, 'off');

end

%plot the individual data points (scatter)
ax1=gca;
hold(ax1, 'all')
x1(1:length(rts_hc_mean))=1;
x2(1:length(rts_pre_mean))=2;
x3(1:length(rts_early_mean))=3;
scatter(x1,rts_hc_mean, 30, 'black')
hold on
scatter(x2, rts_pre_mean, 30, 'black');
hold on
scatter(x3, rts_early_mean, 30, 'black');

%x and y labels
set(gca, 'YLim', [0 2.5], 'FontSize', 14);
ylabel('RT (s)');
xticklabels({'Control', 'Pre-HD', 'Early-HD'})
title('(C) Response time')
box off

%display the test stats:
%calculate effect size:
n2_rt = round(rdivide(tbl_rt{2, 2},tbl_rt{4,2}), 1, 'significant');

if bt_p>0.05
str=['{\it η^2} =', num2str(n2_rt), ',{\it p} = ', num2str(round(pval_rt, 2, 'significant'))];
T = text(min(get(gca, 'xlim')+0.02), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
else
    str=['{\it H}', '(', num2str(tbl_rt{2, 3}), ') = ', num2str(round(tbl_rt{2, 5}, 2)), newline '{\it p} = ' num2str(round(vpa(pval_rt), 8))];
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
end
[c_rt, m_rt] = multcompare(stats_rt, 'display', 'off', 'CType', 'bonferroni');

%Figure 2d:
%load mean values for confidence 
%plot overall mean as bar chart with individual scatter
%conduct parametric assumptions check then appropriate 3-way comparison

load('conf_hc_mean.mat')
load('conf_pre_mean.mat')
load('conf_early_mean.mat')

%calculate mean, sdev and sem
means(:,1)=mean(conf_hc_mean); %1 = hc
means(:,2)=mean(conf_pre_mean); % 2=pre
means(:,3)=mean(conf_early_mean); %3 = early
sdev(:,1) = std(conf_hc_mean, 0, 2); %sdev HC
sem(:,1) = sdev(:,1) / sqrt(length(conf_hc_mean)) %sem HC
sdev(:,2) = std(conf_pre_mean, 0, 2); %sdev pre
sem(:,2) = sdev(:,2) / sqrt(length(conf_pre_mean)) %sem pre
sdev(:,3) = std(conf_early_mean, 0, 2); %sdev early
sem(:,3) = sdev(:,3) / sqrt(length(conf_early_mean)) %sem early

%plot the bars, errors=sem
nexttile
b = bar([1, 2, 3], means, 0.75)
b.FaceColor = 'flat';
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [1 1 1];
b.CData(3,:) = [0.85 0.85 0.85];
hold on
errorbar([1, 2, 3], means, sem,'.', 'color', 'black', 'lineWidth', 1.2, 'CapSize', 30)
%plot the individual data points (scatter)
ax1=gca;
hold(ax1, 'all')
x1(1:length(conf_hc_mean))=1;
x2(1:length(conf_pre_mean))=2;
x3(1:length(conf_early_mean))=3;
scatter(x1,conf_hc_mean, 30, 'black')
hold on
scatter(x2, conf_pre_mean, 30, 'black');
hold on
scatter(x3, conf_early_mean, 30, 'black');

%test equality of variance for KW vs. ANOVA:
%normal distribution?:
norm{:,1}=normalitytest(conf_hc_mean);
norm{:,2}=normalitytest(conf_pre_mean);
norm{:,3}=normalitytest(conf_early_mean);
for g=1:3
    if mean(norm{:,g}(:,3))==1
        disp(sprintf('Data from group %d comes from normal distribution', g));
    else
        disp(sprintf('Data from group %d does NOT come from normal distribution', g));
    end
end
%initialise matrix:
conf(1:length(conf_hc_mean),1:3)=NaN;
conf(:,1)=conf_hc_mean;
conf(1:length(conf_pre_mean),2)=conf_pre_mean;
conf(1:length(conf_early_mean),3)=conf_early_mean;
[bt_p]=vartestn(conf, 'display', 'off') %bartlett's test as data is normal


%group comparison
data=[conf_hc_mean, conf_pre_mean, conf_early_mean]
if bt_p>0.05
    disp('Groups have homogeneity of variance, conducting ANOVA');
%one way ANOVA:
    [pval_conf,tbl_conf,stats_conf] = anova1(data, group, 'off')
else
     disp('NO homogeneity of variance, conducting Kruskall Wallis');
    %kruskall wallis group comparison
[pval_conf,tbl_conf,stats_conf] = kruskalwallis(data, group, 'off')
end
%display the test stats
n2_conf=round(tbl_conf{2,2}/tbl_conf{4,2}, 1, 'significant'); %calculate effect size

if bt_p>0.05
str=['{\it η^2} =', num2str(n2_conf) ',{\it p} = ', num2str(round(pval_conf, 2, 'significant'))];
T = text(min(get(gca, 'xlim')+0.02), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
else
    str=['{\it H}', '(', num2str(tbl_conf{2, 3}), ') = ', num2str(round(tbl_conf{2, 5}, 2)), newline '{\it p} = ' num2str(round(vpa(pval_conf), 8))];
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
end
%x and y labels
ylabel('Level')
xticklabels({'Control', 'Pre-HD', 'Early-HD'})
set(gca, 'YLim', [0 6], 'FontSize', 14);
title('(D) Confidence')
box off
%save the figure:
figuresdir = 'C:\Users\samrc\OneDrive\Documents\GitHub\samrchewitt\HD_perception_metacognition\analysis\behavioural\figs\';
saveas(b, fullfile(figuresdir, 'figure2_acc_dd_rt_conf_bw'), 'fig');
saveas(b, fullfile(figuresdir, 'figure2_acc_dd_rt_conf_bw'), 'jpeg');

%export the figure to pdf: 
%saveas(gcf, 'banalysis.png');

%figure 3 - check staircase stabilisation across blocks:
load('acc_blocks_hc.mat')
load('acc_blocks_pre.mat')
load('acc_blocks_early.mat')

load('acc_blocks_hc.mat')
load('acc_blocks_pre.mat')
load('acc_blocks_early.mat')

accuracy=[acc_blocks_hc; acc_blocks_pre; acc_blocks_early];
n_hc=length(acc_blocks_hc);
n_pre=length(acc_blocks_pre);
n_early=length(acc_blocks_early);
n=58;

%two way analysis of variance (blocks * group):
blocks(:,1)=1:8;
acc_a2 = reshape(accuracy',[],1);
block=repmat(blocks,[58, 1]); %f=blocks
[hc{1:8, 1}]=deal('hc');
[pre{1:8, 1}]=deal('pre');
[early{1:8, 1}]=deal('early');

grp1=repmat(hc, [n_hc, 1]);
grp2=repmat(pre, [n_pre, 1]);
grp3=repmat(early, [n_early, 1]);
group=[grp1; grp2; grp3];

[p2way, tbl2way] = anovan(acc_a2,{block,group}, 'model','interaction','varnames',{'block','group'});

mn(:,1)=mean(accuracy);%1 = early HD 
sd(:,1) = std(accuracy); %standard deviation HD
se(:,1) = sd(:,1) / sqrt(length(accuracy)); %standard error of mean HD


fig=figure;
er1=errorbar(mn(:,1), se(:,1), 'color', 'black', 'lineWidth', 2);
hold on;
a.FaceAlpha=0.2;
box off
hold on
set(gca, 'YLim', [50 100], 'FontSize', 16);
set(gca, 'XLim', [0.75 8.25], 'FontSize', 16);
xlabel('Block')
ylabel('% correct')
title('Accuracy')

%display statistics:
str=['{\it F}', '(', num2str(tbl2way{2, 3}), ',', num2str(tbl2way{5, 3}) ') = ', num2str(round(tbl2way{2, 6}, 2)), ',{\it p} = ', num2str(round(tbl2way{2, 7}, 2))];
T = text(min(get(gca, 'xlim')+0.02), max(get(gca, 'ylim')-1), str); 
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%SAVE
saveas(fig, fullfile(figuresdir, 'accuracy_blocks_bw'), 'jpeg');
