clear
clc
load('confirmed1.mat');
mY = mY0;
load('mR_conf_tdnp_p1q1_for1.mat') %'mFFm',"mFFu","mFFl"

% Pennsylvania (39), Texas (44), California (5)
iT = length(mY)-3;
vT = (1:iT)';
iS1 = iT-30+1; vT2 = (iS1:iT)';vT3=(iS1:95)';
plot(vT(1:end),mY(4:end,iI));xlim([0 110]);
hold on
mF = [mFFl(:,iI) mFFm(:,iI) mFFu(:,iI)];
plot(vT2,mFFm(:,iI))
hold off
