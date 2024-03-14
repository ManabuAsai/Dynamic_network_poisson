% Comparison of forecasts
clear
clc
load('confirmed1.mat');
mY = mY0;
mLY = log(1+mY);

[iT0,iN] = size(mY);
iP=6;
iT = iT0-iP;
[iT iN]
format longG

load('mR_conf_tdnp_p1q1_for1.mat'); %'mFFm',"mFFu","mFFl"
mFF1 = mFFm;
load('mR_conf_ednp_p1q1_for1.mat'); %'mFFm',"mFFu","mFFl"
mFF2 = mFFm;
load('mR_conf_cnp_p1q1_for1.mat'); %'mFFm',"mFFu","mFFl"
mFF3 = mFFm;
load('mR_conf_tdnpi_p1q1_for1.mat'); %'mFFm',"mFFu","mFFl"
mFF4 = mFFm;
load('mR_conf_ednpi_p1q1_for1.mat'); %'mFFm',"mFFu","mFFl"
mFF5 = mFFm;
load('mR_conf_cnpi_p1q1_for1.mat'); %'mFFm',"mFFu","mFFl"
mFF6 = mFFm;

mY_true = mY(end-iF+1:end,:);
vMSFE1 = sum((mFF1-mY_true).^2,2);
vMSFE2 = sum((mFF2-mY_true).^2,2);
vMSFE3 = sum((mFF3-mY_true).^2,2);
vMSFE4 = sum((mFF4-mY_true).^2,2);
vMSFE5 = sum((mFF5-mY_true).^2,2);
vMSFE6 = sum((mFF6-mY_true).^2,2);
mMSFE = [vMSFE1 vMSFE2 vMSFE3 vMSFE4 vMSFE5 vMSFE6];
vRMSFE = sqrt(mean(mMSFE))';
vRMSFE


