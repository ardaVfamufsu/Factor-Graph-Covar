%% 
% Boxplots of variance estimates with biased and unbiased methods - linear model
% This chunk runs main_GpsLocalize_MeasAndTwoProcessCov.m in batch mode.
% Comment the initial (user inputs) section of the above code to run it in batch mode
clear all;
close all;
flag_robust=0; 
flag_outliers=0; 
flag_scale_cov_change=1;
n_MC=1000; 
n_rep_cycle=1;
% unbiased
flag_unbiased=1;
rng(0);
main_GpsLocalize_MeasAndTwoProcessCov;
var_vect_unbiased=var_vect;
% biased 
flag_unbiased=0;
rng(0);
main_GpsLocalize_MeasAndTwoProcessCov;
var_vect_biased=var_vect;


figure;
data = {var_vect_unbiased',  var_vect_biased'};
boxplotGroup(data, 'PrimaryLabels', {'Prop' 'Exis'}, ...
  'SecondaryLabels',{'\sigma^2_r', '\sigma^2_{q1}', '\sigma^2_{q2}'}, 'InterGroupSpace', 1,'symbol','')
hold on;
plot([.7 2.3],[sig_e sig_e],'g--','linewidth',2)
plot([3.8 5.4],[sig_p1 sig_p1],'g--','linewidth',2)
plot([6.8 8.4],[sig_p2 sig_p2],'g--','linewidth',2)
saveas(gcf,'VarEstBoxPlotLinModel.jpg');
%% 
% Boxplots of unbiased variance estimates with outliers and linear model, with and without M estimator 
% This chunk runs main_GpsLocalize_MeasAndTwoProcessCov.m. in batch mode.
% Comment the initial (user inputs) section of the above code to run it in batch mode
% Use w1=0.90 and w1=0.75 in the above code to consider the two different outlier cases
clear all;
close all;
flag_outliers=1;
flag_unbiased=1;
flag_scale_cov_change=1;
n_MC=1000; 
n_rep_cycle=1;
% M estimator
flag_robust=1; 
rng(0);
disp('Outliers, using M estimator');
main_GpsLocalize_MeasAndTwoProcessCov;
var_vect_robust=var_vect;
% no M estimator 
flag_robust=0; 
rng(0);
disp('Outliers, not using M estimator');
main_GpsLocalize_MeasAndTwoProcessCov;
var_vect_NonRobust=var_vect;

figure;
data = {var_vect_robust',  var_vect_NonRobust'};
boxplotGroup(data, 'PrimaryLabels', {'Prop-M' 'Prop'}, ...
  'SecondaryLabels',{'\sigma^2_r', '\sigma^2_{q1}', '\sigma^2_{q2}'}, 'InterGroupSpace', 1,'symbol','')
hold on;
plot([.7 2.3],[sig_e sig_e],'g--','linewidth',2)
plot([3.8 5.4],[sig_p1 sig_p1],'g--','linewidth',2)
plot([6.8 8.4],[sig_p2 sig_p2],'g--','linewidth',2)
%saveas(gcf,'VarEstBoxPlotLinModelOutliers.jpg');
%% 
% boxplot variance estimates with biased and unbiased methods - nonlinear
% model
% run in batch mode, main_GpsLocalize_MeasAndTwoProcessCov.m
% comment lines for user inputs
clear all;
close all;
flag_robust=0; 
flag_outliers=0; 
flag_scale_cov_change=1;
n_MC=1000; 
n_rep_cycle=1;
% unbiased
flag_unbiased=1;
rng(2);
main_GpsLocalize_MeasAndTwoProcessCov_NonlinearModel;
var_vect_unbiased=var_vect;
% biased 
flag_unbiased=0;
rng(2);
main_GpsLocalize_MeasAndTwoProcessCov_NonlinearModel;
var_vect_biased=var_vect;

figure;
data = {var_vect_unbiased',  var_vect_biased'};
boxplotGroup(data, 'PrimaryLabels', {'Prop' 'Exis'}, ...
  'SecondaryLabels',{'\sigma^2_r', '\sigma^2_{q1}', '\sigma^2_{q2}'}, 'InterGroupSpace', 1,'symbol','')
hold on;
plot([.7 2.3],[sig_e sig_e],'g--','linewidth',2)
plot([3.8 5.4],[sig_p1 sig_p1],'g--','linewidth',2)
plot([6.8 8.4],[sig_p2 sig_p2],'g--','linewidth',2)
%saveas(gcf,'VarEstBoxPlotNonLinModel.jpg');

%% 
% boxplot for data size: variance estimates with unbiased methods linear model
% run in batch mode, main_GpsLocalize_MeasAndTwoProcessCov.m
% comment lines for user inputs
clear all;
close all;
flag_robust=0; 
flag_outliers=0; 
flag_scale_cov_change=1;
n_MC=1000; 
% nrep =1, unbiased
flag_unbiased=1;
n_rep_cycle=1;
rng(0);
main_GpsLocalize_MeasAndTwoProcessCov;
var_vect_unbiasedRep1=var_vect;
clear Yn var_vect;
% nrep =5, unbiased 
flag_unbiased=1;
n_rep_cycle=5;
rng(0);
main_GpsLocalize_MeasAndTwoProcessCov;
var_vect_unbiasedRep5=var_vect;


figure;
data = {var_vect_unbiasedRep1',  var_vect_unbiasedRep5'};
boxplotGroup(data, 'PrimaryLabels', {'n=20', 'n=100'}, ...
  'SecondaryLabels',{'\sigma^2_r', '\sigma^2_{q1}', '\sigma^2_{q2}'}, 'InterGroupSpace', 1,'symbol','')
hold on;
plot([.7 2.3],[sig_e sig_e],'g--','linewidth',2)
plot([3.8 5.4],[sig_p1 sig_p1],'g--','linewidth',2)
plot([6.8 8.4],[sig_p2 sig_p2],'g--','linewidth',2)
%saveas(gcf,'VarEstBoxPlotLinModelDataSize.jpg');


