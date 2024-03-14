clear
clc
load('confirmed1.mat');
mY = mY0;
mLY = log(1+mY);

[iT0,iN] = size(mY);
iP=6;
iT = iT0-iP;
[iT iN]


iP=1;
iT=81-iP;iF=30;
mY0 = mY(end-iT-iF+1:end,:); 

% Pennsylvania (39), Texas (44), California (5)

[iT1,~]=size(mY0);[iT1 iN]
vT =(1:iT1)';
plot(vT,[mY0(:,5) mY0(:,44) mY0(:,39)]);xlim([0 110])
legend('California', 'Texas','Pennsylvania','Location','northwest');

mX = [mY0(:,5) mY0(:,44) mY0(:,39)];
ts1 = timeseries(mX(:,1),vT);
ts1.Name = 'California';
ts1.TimeInfo.Units = 'weeks';
ts1.TimeInfo.StartDate = '29-Jan-2020';
ts1.TimeInfo.Format = 'mm/dd/yyyy';
ts1.Time = ts1.Time - ts1.Time(1);
plot(ts1)
title(' ')
ylabel('Count')

ts2 = timeseries(mX(:,2),vT);
ts2.Name = 'Texas';
ts2.TimeInfo.Units = 'weeks';
ts2.TimeInfo.StartDate = '27-Jan-2020';
ts2.TimeInfo.Format = 'mm/dd/yyyy';
ts2.Time = ts2.Time - ts2.Time(1);
hold on
plot(ts2)

ts3 = timeseries(mX(:,3),vT);
ts3.Name = 'Texas';
ts3.TimeInfo.Units = 'weeks';
ts3.TimeInfo.StartDate = '27-Jan-2020';
ts3.TimeInfo.Format = 'mm/dd/yyyy';
ts3.Time = ts3.Time - ts3.Time(1);
hold on
plot(ts3)

legend('California', 'Texas','Pennsylvania','Location','northwest');

