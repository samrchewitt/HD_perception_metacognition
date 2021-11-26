%Sam Hewitt 2021
% see https://github.com/metacoglab/HMeta-d/tree/master/Matlab for
% dependencies for model fitting 
% turn off model fitting in jobs 
% by default script loads previous fits, plots results and runs GLM
% regressions
clear all; close all; clc 
%% PATHS:
datapath=['C:\Users\samrc\OneDrive\Documents\GitHub\samrchewitt\HD_perception_metacognition\analysis\Hmetad_metacognition\summary_data\'];
behaviouraldata=['C:\Users\samrc\OneDrive\Documents\GitHub\samrchewitt\HD_perception_metacognition\analysis\behavioural\summary_data\'];
clinicaldata=['C:\Users\samrc\OneDrive\Documents\GitHub\samrchewitt\HD_perception_metacognition\raw_data\clinical'];
addpath(genpath(datapath)); addpath(genpath(behaviouraldata)); addpath(genpath(clinicaldata));
addpath(genpath('D:\meta-hd\analysis'));  
addpath(genpath('D:\matlab'));  
modelspath=['C:\Users\samrc\OneDrive\Documents\GitHub\samrchewitt\HD_perception_metacognition\analysis\Hmetad_metacognition\fits\final'];
figspath=['C:\Users\samrc\OneDrive\Documents\GitHub\samrchewitt\HD_perception_metacognition\analysis\Hmetad_metacognition\fits\final\figs'];

%% jobs:
runHmeta=0;
savemodels=0;
plotdata=1;
runGLM=1;
%% model fitting %%%
if runHmeta
% load the confidence rating vectors for each group
load(['nR_S1prehd']); load(['nR_S2prehd']);
load(['nR_S1earlyhd']); load(['nR_S2earlyhd']);
load(['nR_S1hc']); load(['nR_S2hc']);

%fit pre-HD group
fitpreHD = fit_meta_d_mcmc_group(nR_S1prehd, nR_S2prehd);
%fit early-HD group
fitearlyHD = fit_meta_d_mcmc_group(nR_S1earlyhd, nR_S2earlyhd);
%fit HC matched-group
fitHC = fit_meta_d_mcmc_group(nR_S1hc, nR_S2hc);

%%% otherwise just load the previous fits:
else 
    addpath(modelspath);
    load('fitHC.mat'); load('fitpreHD.mat'); load('fitearlyHD.mat');
end
%%%% also load accuracy data: 
load('correct_hc.mat'); load('correct_pre.mat'); load('correct_early.mat');
%% plot data %%%
if plotdata 
set(0,'defaultfigurecolor',[1 1 1]); n=58; %define participants 
%%% statistical differences between distributions in pair-wise comparison: 
%PRE-HD vs. CONTROL
sampleDiff_prehc = fitpreHD.mcmc.samples.mu_logMratio - fitHC.mcmc.samples.mu_logMratio;
hdi_prehc = calc_HDI(sampleDiff_prehc(:)); ci_prehc = calc_CI(sampleDiff_prehc(:));
p_prehc=mean(sampleDiff_prehc(:) < 0);

%EARLY-HD vs. CONTROL
sampleDiff_earlyhc = fitearlyHD.mcmc.samples.mu_logMratio - fitHC.mcmc.samples.mu_logMratio;
hdi_earlyhc = calc_HDI(sampleDiff_earlyhc(:)); ci_earlyhc = calc_CI(sampleDiff_earlyhc(:));
p_earlyhc=mean(sampleDiff_earlyhc(:) < 0);

%EARLY-HD vs. PRE-HD
sampleDiff_earlypre = fitearlyHD.mcmc.samples.mu_logMratio - fitpreHD.mcmc.samples.mu_logMratio;
hdi_earlypre = calc_HDI(sampleDiff_earlypre(:)); ci_earlypre = calc_CI(sampleDiff_earlypre(:));
p_earlypre=mean(sampleDiff_earlypre(:) < 0);

%calc HDIs of distributions alone: 
hdi_prehd=calc_HDI(exp(fitpreHD.mcmc.samples.mu_logMratio(:)))
hdi_earlyhd=calc_HDI(exp(fitearlyHD.mcmc.samples.mu_logMratio(:)))
hdi_hcmatched=calc_HDI(exp(fitHC.mcmc.samples.mu_logMratio(:)))

%%%% PLOT %%%%
plotSamples_hdi(sampleDiff_prehc, hdi_prehc, p_prehc, [0 0 0])
ylabel('Samples'); title('Difference pre-HD vs. HC');
plotSamples_hdi(sampleDiff_earlyhc, hdi_earlyhc, p_earlyhc, [0 0 0])
ylabel('Samples'); title('Difference early-HD vs. HC')

%%%% PLOT 3 GROUPS TOGETHER %%%%
plotSamples3grps(exp(fitpreHD.mcmc.samples.mu_logMratio), exp(fitearlyHD.mcmc.samples.mu_logMratio), exp(fitHC.mcmc.samples.mu_logMratio))
ylabel('Sample count', 'FontSize', 20); xlabel('M-ratio', 'FontSize', 20)
box off; legend boxoff; legend('Location', 'WestOutside');
%%% SAVE IT %%%:
exportgraphics(gcf,[figspath 'figure3_plotHmetaSamples.jpeg'],'Resolution',300)

%%% PRINT mean values %%%
display(['Control mean M-ratio = ' num2str(mean(fitHC.Mratio))]);
display(['PRE-HD mean M-ratio = ' num2str(mean(fitpreHD.Mratio))]);
display(['EARLY-HD mean M-ratio = ' num2str(mean(fitearlyHD.Mratio))]);

%%%% PLOT M-ratio mean values and accuracy together %%%%
mratio_hc=fitHC.Mratio; mratio_pre=fitpreHD.Mratio; mratio_early=fitearlyHD.Mratio;
%scale accuracy to 1:
correct_hc=correct_hc/100; correct_pre=correct_pre/100; correct_early=correct_early/100;
%vectorise both: 
correct=[correct_hc; correct_pre; correct_early];
mratio=[mratio_hc'; mratio_pre'; mratio_early'];

%%build artificial x-axis of participants
xhc=(1:length(correct_hc))';
xpre=(length(correct_hc)+1:(length(correct_hc)+length(correct_pre)))';
xearly=(xpre(end)+1:xpre(end)+length(correct_early))'; %early

%PUT the data together (col 1=mratio, col2=group, col3=accuracy); 
mratio(:,2)=zeros; 
mratio(xhc,2)=1; mratio(xpre,2)=2; mratio(xearly,2)=3; 
mratio(:,3)=correct;
% SORT based on m-ratio
smratio=sortrows(mratio, 1); %sort based on M-ratio

%figure parameters: 
x=(1:n)'; set(0,'defaultfigurecolor','w')
figure('Renderer', 'painters', 'Position', [0 0 1000 700]); 
labels={'Accuracy', 'Control M-ratio', 'Pre-HD M-ratio', 'Early-HD M-ratio'};
rgb=[0 0 0; 0 0 0; 0 0 0]; %marker colors
mrks1=['.'; '.'; 'o']; mrks2=['o'; '^'; 'p']; %markers types
sz1=[18, 18, 4]; sz2=[4, 6, 8]; %marker sizes
% and PLOT: 
plot(x,smratio(:,3), ':o', 'color', [0.5 0.5 0.5], 'MarkerSize', 2, 'LineWidth', 2) %plot accuracy
hold on; gscatter(x,smratio(:,1), smratio(:,2), rgb, mrks2, sz2) %plot mratios
hold on; 
%appearance
xlabel('Participants', 'FontSize', 20); lgnd = legend(labels, 'Location', 'WestOutside');
set(gca, 'XLim', [-1 n+1], 'FontSize', 20);
set(gca, 'YLim', [min(smratio(:,1))-0.1, max(smratio(:,1))+0.1], 'FontSize', 20);
set(gca,'xtick',[]); set(gca,'color','w');
box off; legend boxoff 
exportgraphics(gcf,[figspath 'Figure4_accuracy_mratio_individual.jpeg'],'Resolution',300)


end
%% run GLM %%%
if runGLM
%load clinical/behavioural data;
load([clinicaldata '\metaHD_clinicaldata.mat']); load('id_hmetaD');
load('tBehaviour.mat');

%%%MEANS FROM FITS %%%%:
%mratio (metacognitive efficiency) 
load('fitHC.mat'); mratio_hc=fitHC.Mratio; 
load('fitpreHD.mat'); mratio_pre=fitpreHD.Mratio;
load('fitearlyHD.mat'); mratio_early=fitearlyHD.Mratio;
mratio=[mratio_pre'; mratio_early'; mratio_hc']; %add to table
%dprime (perceptual sensitivity)
dphc=fitHC.d1; dppre=fitpreHD.d1; dpearly=fitearlyHD.d1;
dp=[dppre'; dpearly'; dphc']; %add to table
%metad-prime (metacognitive sensitivity):
meta_dhc=fitHC.meta_d; meta_dpre=fitpreHD.meta_d; meta_dearly=fitearlyHD.meta_d;
mdp=[meta_dpre'; meta_dearly'; meta_dhc']; %add to table

%%%TABULATE metacog parameters:
tMeta=table(id_hmeta, mratio, dp, mdp);

%%%% then merge clinical and Meta using ID as key
tM=innerjoin(tM, tMeta, 'Keys', 1);
tM=innerjoin(tM, tB, 'Keys', 1);

%calculate z scores for continuous predictors:
tM.Zage=zscore(tM.age_years);
tM.Znart=nanzscore(tM.nart_iq);
tM.Zmmse=nanzscore(tM.mmse_total);
tM.Zuhdrs=nanzscore(tM.uhdrs_motor_score);
tM.ZhadsA=nanzscore(tM.hads_a);
tM.ZhadsD=nanzscore(tM.hads_d);
%convert categoricals:
tM.sex = categorical(tM.sex);
tM.gene = categorical(tM.gene);

%%%% run regressions %%%%

dv = {'mratio', 'dp', 'mdp', 'conf'};
for d = 1:length(dv)
    modelspec = [dv{d} ' ~ gene + Zage + sex + Znart + Zmmse + ZhadsA + ZhadsD'];
    mdl{d} = fitglm(tM,modelspec);
    disp(dv{d})
    disp(mdl{d})
    fprintf('\n\n\n')
end

%plot the coefficients -MRATIO:
figure('Renderer', 'painters', 'Position', [0 0 1000 500]); 
bar(mdl{1, 1}.Coefficients.Estimate(2:end), 'w', 'LineWidth', 1.5)
hold on 
er=errorbar(mdl{1, 1}.Coefficients.Estimate(2:end), mdl{1, 1}.Coefficients.SE(2:end));
er.Color = [0 0 0]; er.LineStyle = 'none'; 
title('M-ratio', 'FontSize', 20)
x={'HD+', 'Gender (Male)', 'Age', 'IQ', 'MMSE', 'HADS-A', 'HADS-D'};
set(gca,'xticklabels', x, 'FontSize', 14); ylim([-0.1, 0.25]); 
ylabel('Regression coeffient', 'FontSize', 18);
%add signifance: 
mysigstar(gca, [1 1], 0.175, mdl{1, 1}.Coefficients.pValue(2), 'black');
mysigstar(gca, [2 2], 0.15, mdl{1, 1}.Coefficients.pValue(3), 'black');
mysigstar(gca, [3 3], 0.125, mdl{1, 1}.Coefficients.pValue(4), 'black');
mysigstar(gca, [4 7], 0.05, mdl{1, 1}.Coefficients.pValue(5), 'black');
%add stat to figure:
p=round(coefTest(mdl{1,1}), 3, 'significant');
str=['{\it R^2} = ', num2str(round(mdl{1, 1}.Rsquared.Ordinary, 2)), ', {\it p} = ', num2str(round(p, 1, 'significant'))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
box off
exportgraphics(gcf,[figspath 'Figure5_Mratio_GLM.jpeg'],'Resolution',300)


%plot the coefficients -DPRIME:
figure; bar(mdl{1, 2}.Coefficients.Estimate(2:end), 'w'); hold on 
er=errorbar(mdl{1, 2}.Coefficients.Estimate(2:end), mdl{1, 2}.Coefficients.SE(2:end));
er.Color = [0 0 0]; er.LineStyle = 'none'; 
title('Perceptual sensitivity (d-prime)'); ylim([-0.2, 0.2]);
xticklabels(x); set(gca,'xticklabels', x, 'FontSize', 18);
ylabel('Regression coeffient', 'FontSize', 18);
%add stat to figure:
p=round(coefTest(mdl{1,2}), 3, 'significant');
str=['{\it R-Sq} = ', num2str(mdl{1, 2}.Rsquared.Ordinary), newline, '{\it R-Sq adj} = ', num2str(mdl{1, 2}.Rsquared.Adjusted), newline, '{\it p} = ', num2str(p)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
box off

%plot the coefficients -METAD-PRIME:
figure; bar(mdl{1, 3}.Coefficients.Estimate(2:end), 'w'); hold on 
er=errorbar(mdl{1, 3}.Coefficients.Estimate(2:end), mdl{1, 3}.Coefficients.SE(2:end));
er.Color = [0 0 0]; er.LineStyle = 'none'; 
title('Metacognitive sensitivity');
set(gca,'xticklabels', x, 'FontSize', 18);
%add significance: 
mysigstar(gca, [1 1], 0.2, mdl{1, 3}.Coefficients.pValue(2), 'black');
mysigstar(gca, [2 2], 0.22, mdl{1, 3}.Coefficients.pValue(3), 'black');
mysigstar(gca, [3 7], 0.1, mdl{1, 3}.Coefficients.pValue(4), 'black');
ylim([-0.1, 0.3]);
ylabel('Regression coeffient', 'FontSize', 18);
%add stat to figure: 
p=round(coefTest(mdl{1,3}), 3, 'significant');
str=['{\it R-Sq} = ', num2str(mdl{1, 3}.Rsquared.Ordinary), newline, '{\it R-Sq adj} = ', num2str(mdl{1, 3}.Rsquared.Adjusted), newline, '{\it p} = ', num2str(p)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
box off

%plot the coefficients - CONFIDENCE : 
figure; bar(mdl{1, 4}.Coefficients.Estimate(2:end), 'w'); hold on 
er=errorbar(mdl{1, 4}.Coefficients.Estimate(2:end), mdl{1, 4}.Coefficients.SE(2:end));
er.Color = [0 0 0]; er.LineStyle = 'none'; 
title('Confidence');
ylabel('Regression coeffient', 'FontSize', 18);
x={'HD+', 'Gender (male)', 'Age', 'IQ', 'MMSE', 'HADS-A', 'HADS-D'};
set(gca,'xticklabels', x, 'FontSize', 18);
ylim([-0.6, 0.6]);
%add stat to figure:
p=round(coefTest(mdl{1,4}), 3, 'significant');
str=['{\it R-Sq} = ', num2str(mdl{1, 4}.Rsquared.Ordinary), newline, '{\it R-Sq adj} = ', num2str(mdl{1, 4}.Rsquared.Adjusted), newline, '{\it p} = ', num2str(p)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
box off
end

%%%% SAVE %%%
if savemodels
mkdir(modelspath);
save(fullfile(modelspath, 'fitHC'), 'fitHC');
save(fullfile(modelspath, 'fitpreHD'), 'fitpreHD');
save(fullfile(modelspath, 'fitearlyHD'), 'fitearlyHD');
end